import unittest
import parameterized

import jax
import jax.numpy as jnp
import meep as mp
import meep.adjoint as mpa
import numpy as onp

# The calculation of finite difference gradients requires that JAX be operated with double precision
jax.config.update('jax_enable_x64', True)

# The step size for the finite difference gradient calculation
_FD_STEP = 1e-4

# The relative tolerance for the adjoing and finitie difference gradient comparison
_RTOL = 5e-2

mp.verbosity(0)


def build_straight_wg_simulation(
  wg_width=0.7,
  wg_padding=1.0,
  wg_length=1.0,
  pml_width=1.0,
  source_to_pml=0.5,
  source_to_monitor=0.1,
  frequencies=[1 / 1.55],
  gaussian_rel_width=0.4,
  sim_resolution=20,
  design_region_resolution=20,
):
  """Builds a simulation of a straight waveguide with a design region segment."""
  design_region_shape = (1.0, wg_width)

  # Simulation domain size
  sx = 2 * pml_width + 2 * wg_length + design_region_shape[0]
  sy = 2 * pml_width + 2 * wg_padding + max(
    wg_width,
    design_region_shape[1],
  )

  # Mean / center frequency
  fmean = onp.mean(frequencies)

  si = mp.Medium(index=3.4)
  sio2 = mp.Medium(index=1.44)

  sources = [
    mp.EigenModeSource(
      mp.GaussianSource(frequency=fmean, fwidth=fmean * gaussian_rel_width),
      eig_band=1,
      direction=mp.NO_DIRECTION,
      eig_kpoint=mp.Vector3(1, 0, 0),
      size=mp.Vector3(0, wg_width + 2 * wg_padding, 0),
      center=[-sx / 2 + pml_width + source_to_pml, 0, 0],
    )
  ]

  nx = int(design_region_resolution * design_region_shape[0])
  ny = int(design_region_resolution * design_region_shape[1])
  mat_grid = mp.MaterialGrid(
    mp.Vector3(nx, ny),
    sio2,
    si,
    grid_type='U_DEFAULT',
  )

  design_regions = [
    mpa.DesignRegion(
      mat_grid,
      volume=mp.Volume(
        center=mp.Vector3(),
        size=mp.Vector3(
          design_region_shape[0],
          design_region_shape[1],
          0,
        ),
      ),
    )
  ]

  geometry = [
    mp.Block(
      center=mp.Vector3(x=-design_region_shape[0] / 2 - wg_length / 2 - pml_width / 2),
      material=si,
      size=mp.Vector3(wg_length + pml_width, wg_width, 0)),  # left wg
    mp.Block(
      center=mp.Vector3(x=+design_region_shape[0] / 2 + wg_length / 2 + pml_width / 2),
      material=si,
      size=mp.Vector3(wg_length + pml_width, wg_width, 0)),  # right wg
    mp.Block(
      center=design_regions[0].center,
      size=design_regions[0].size,
      material=mat_grid),  # design region
  ]

  simulation = mp.Simulation(
    cell_size=mp.Vector3(sx, sy),
    boundary_layers=[mp.PML(pml_width)],
    geometry=geometry,
    sources=sources,
    resolution=sim_resolution,
  )

  monitor_centers = [
    mp.Vector3(-sx / 2 + pml_width + source_to_pml + source_to_monitor),
    mp.Vector3(sx / 2 - pml_width - source_to_pml - source_to_monitor),
  ]
  monitor_size = mp.Vector3(y=wg_width + 2 * wg_padding)

  monitors = [
    mpa.EigenmodeCoefficient(
      simulation,
      mp.Volume(center=center, size=monitor_size),
      mode=1,
      forward=True) for center in monitor_centers
  ]
  return simulation, sources, monitors, design_regions, frequencies


class UtilsTest(unittest.TestCase):
  def setUp(self):
    super().setUp()
    (
      self.simulation,
      self.sources,
      self.monitors,
      self.design_regions,
      self.frequencies,
    ) = build_straight_wg_simulation()

  def test_mode_monitor_helpers(self):
    mpa.jax.utils.register_monitors(self.monitors, self.frequencies)
    self.simulation.run(until=100)
    monitor_values = mpa.jax.utils.gather_monitor_values(self.monitors)
    self.assertEqual(monitor_values.dtype, onp.complex128)
    self.assertEqual(monitor_values.shape,
                     (len(self.monitors), len(self.frequencies)))

  def test_design_region_monitor_helpers(self):
    design_region_monitors = mpa.jax.utils.install_design_region_monitors(
      self.simulation,
      self.design_regions,
      self.frequencies,
    )
    self.simulation.run(until=100)
    design_region_fields = mpa.jax.utils.gather_design_region_fields(
      self.simulation,
      design_region_monitors,
      self.frequencies,
    )

    self.assertIsInstance(design_region_fields, list)
    self.assertEqual(len(design_region_fields), len(self.design_regions))

    self.assertIsInstance(design_region_fields[0], list)
    self.assertEqual(len(design_region_fields[0]), len(mpa.jax.utils._ADJOINT_FIELD_COMPONENTS))

    for value in design_region_fields[0]:
      self.assertIsInstance(value, onp.ndarray)
      self.assertEqual(value.ndim, 4)  # dims: freq, x, y, pad
      self.assertEqual(value.dtype, onp.complex128)


class WrapperTest(unittest.TestCase):

  @parameterized.parameterized.expand([
    ('1550_1600bw_01relative_gaussian_10fill', onp.linspace(1 / 1.55, 1 / 1.60, 3).tolist(), 0.1, 1.0),
    ('1500_1600bw_02relative_gaussian_10fill', onp.linspace(1 / 1.55, 1 / 1.60, 4).tolist(), 0.2, 1.0),
    ('1600_1680bw_03relative_gaussian_10fill', onp.linspace(1 / 1.60, 1 / 1.68, 5).tolist(), 0.3, 1.0),
    ('1550_1600bw_01relative_gaussian_05fill', onp.linspace(1 / 1.55, 1 / 1.60, 3).tolist(), 0.1, 0.5),
    ('1500_1600bw_02relative_gaussian_05fill', onp.linspace(1 / 1.55, 1 / 1.60, 4).tolist(), 0.2, 0.5),
    ('1600_1680bw_03relative_gaussian_05fill', onp.linspace(1 / 1.60, 1 / 1.68, 5).tolist(), 0.3, 0.5),
  ])
  def test_wrapper_gradients(self, _, frequencies, gaussian_rel_width, design_variable_fill_value):
    """Tests gradient from the JAX-Meep wrapper against finite differences."""
    (
      simulation,
      sources,
      monitors,
      design_regions,
      frequencies,
    ) = build_straight_wg_simulation(frequencies=frequencies, gaussian_rel_width=gaussian_rel_width)

    wrapped_meep = mpa.jax.MeepJaxWrapper(
      simulation,
      sources,
      monitors,
      design_regions,
      frequencies,
      measurement_interval=50.0,
      dft_field_components=(mp.Ez,),
      dft_threshold=1e-6,
      minimum_run_time=0,
      maximum_run_time=onp.inf,
      until_after_sources=True
    )

    design_shape = tuple(int(i) for i in design_regions[0].design_parameters.grid_size)[:2]
    x = onp.ones(design_shape) * design_variable_fill_value

    # Define a loss function
    def loss_fn(x):
      monitor_values = wrapped_meep([x])
      t = monitor_values[1, :] / monitor_values[0, :]
      # Mean transmission vs wavelength
      return jnp.mean(jnp.square(jnp.abs(t)))

    value, adjoint_grad = jax.value_and_grad(loss_fn)(x)

    # Create dp
    random_perturbation_vector = _FD_STEP * onp.random.random(x.shape)

    # Calculate p + dp
    x_perturbed = x + random_perturbation_vector

    # Calculate T(p + dp)
    value_perturbed = loss_fn(x_perturbed)

    # Check that dp . ∇T ~ T(p + dp) - T(p)
    onp.testing.assert_allclose(
      onp.dot(
        random_perturbation_vector.ravel(),
        adjoint_grad.ravel(),
      ),
      value_perturbed - value,
      rtol=_RTOL,
      err_msg='The projection of the gradient computed via the adjoint simulation into a '
              'random direction of the design parameter space does not agree with a finite '
              'difference approximation of that same projection to the specified tolerance.'
    )


if __name__ == '__main__':
  unittest.main()

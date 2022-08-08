"""Wrapper for converting a Meep simulation into a differentiable JAX callable function.

Usage example:
```
import jax.numpy as jnp
import meep as mp
import meep.adjoint as mpa

sources = [
  mp.EigenModeSource(...)
]

monitors = [
  mpa.EigenmodeCoefficient(...),
  mpa.EigenmodeCoefficient(...),
]

design_regions = [
  mpa.DesignRegion(...)
]

frequencies = [1/1.55, 1/1.60, 1/1.65, ...]

simulation = mp.Simulation(...)

wrapped_meep = MeepJaxWrapper(
    simulation,
    sources,
    monitors,
    design_regions,
    frequencies,
    measurement_interval = 50.0,
    dft_field_components = (mp.Ez,),
    dft_threshold = 1e-6,
    minimum_run_time = 0,
    maximum_run_time = onp.inf,
    until_after_sources = True
)

def loss(x):
    monitor_values = wrapped_meep([x])
    t = monitor_values[0,:] / monitor_values[1,:]
    # Mean transmission vs wavelength
    return jnp.mean(jnp.abs(t))

value, grad = jax.value_and_grad(loss)(x)
```
"""
from typing import Callable, Iterable, List, Tuple

import jax
import jax.numpy as jnp
import numpy as onp

import meep as mp

from . import DesignRegion, EigenmodeCoefficient, utils

_norm_fn = onp.linalg.norm
_reduce_fn = onp.max


class MeepJaxWrapper:
    """Wraps a Meep simulation object into a JAX-differentiable callable.

    Attributes:
        simulation: the pre-configured Meep simulation object.
        sources: a list of Meep sources for the forward simulation.
        monitors: a list of eigenmode coefficient monitors from the `meep.adjoint`
          module.
        design_regions: a list of design regions from the `meep.adjoint` module.
        frequencies: the list of frequencies, in normalized Meep units.
        measurement_interval: the time interval between DFT field convergence
          measurements, in Meep time units. The default value is 50.
        dft_field_components: a list of Meep field components, such as `mp.Ex`,
          `mp.Hy`, etc, whose DFT will be monitored for convergence to stop the
          simulation. The default is `mp.Ez`.
        dft_threshold: the threshold for DFT field convergence. Once the norm of the
          change in the fields (the maximum over all design regions and field
          components) is less than this value, the simulation will be stopped. The
          default value is 1e-6.
        minimum_run_time: the minimum run time of the simulation, in Meep time
          units. The default value is 0.
        maximum_run_time: the maximum run time of the simulation, in Meep time
          units. The default value is infinity.
        until_after_sources: whether `maximum_run_time` should be ignored until the
          sources have turned off. This parameter specifies whether `until` or
          `until_after_sources` is used. See
          https://meep.readthedocs.io/en/latest/Python_User_Interface/#Simulation
          for more information. The default is true.
    """

    _log_fn = print

    def __init__(
        self,
        simulation: mp.Simulation,
        sources: List[mp.Source],
        monitors: List[EigenmodeCoefficient],
        design_regions: List[DesignRegion],
        frequencies: List[float],
        dft_threshold: float = 1e-11,
        minimum_run_time: float = 0,
        maximum_run_time: float = onp.inf,
        until_after_sources: bool = True,
        finite_difference_step: float = utils.FD_DEFAULT,
    ):
        self.simulation = simulation
        self.sources = sources
        self.monitors = monitors
        self.design_regions = design_regions
        self.frequencies = frequencies
        self.dft_threshold = dft_threshold
        self.minimum_run_time = minimum_run_time
        self.maximum_run_time = maximum_run_time
        self.until_after_sources = until_after_sources
        self.finite_difference_step = finite_difference_step

        self._simulate_fn = self._initialize_callable()

    def __call__(self, designs: List[jnp.ndarray]) -> jnp.ndarray:
        """Performs a Meep simulation, taking a list of designs and returning mode overlaps.

        Args:
          designs: a list of design variables as 1D, 2D, or 3D JAX arrays. Valid shapes for
          design variables are (Nx, Ny, Nz) where Nx{y,z} match the elements of the
          `grid_size` constructor argument of Meep's `MaterialGrid` used for the
          corresponding design region. Singleton dimensions of the `grid_size` may be
          omitted from the corresponding design variable. For example, a design variable
          with a shape of either (10, 20) or (10, 20, 1) would be compatible with a
          `grid_size` of (10, 20, 1). Similarly, a design variable with shapes of (25,),
          (25, 1), or (25, 1, 1) would be compatible with a `grid_size` of (25, 1, 1).

        Returns:
          a complex-valued JAX ndarray of differentiable mode monitor overlaps values with
          a shape of (num monitors, num frequencies).
        """
        return self._simulate_fn(designs)

    def _run_fwd_simulation(
        self, design_variables: Iterable[onp.ndarray]
    ) -> (jnp.ndarray, List[List[mp.DftFields]]):
        """Runs forward simulation, returning monitor values and design region fields."""
        utils.validate_and_update_design(self.design_regions, design_variables)
        self.simulation.reset_meep()
        self.simulation.change_sources(self.sources)
        utils.register_monitors(self.monitors, self.frequencies)
        fwd_design_region_monitors = utils.install_design_region_monitors(
            self.simulation,
            self.design_regions,
            self.frequencies,
        )
        self.simulation.init_sim()
        sim_run_args = {
            "until_after_sources"
            if self.until_after_sources
            else "until": mp.stop_when_dft_decayed(
                self.dft_threshold, self.minimum_run_time, self.maximum_run_time
            )
        }
        self.simulation.run(**sim_run_args)

        monitor_values = utils.gather_monitor_values(self.monitors)
        return (jnp.asarray(monitor_values), fwd_design_region_monitors)

    def _run_adjoint_simulation(
        self, monitor_values_grad: onp.ndarray
    ) -> List[List[mp.DftFields]]:
        """Runs adjoint simulation, returning design region fields."""
        if not self.design_regions:
            raise RuntimeError(
                "An adjoint simulation was attempted when no design "
                "regions are present."
            )
        adjoint_sources = utils.create_adjoint_sources(
            self.monitors, monitor_values_grad
        )
        # TODO refactor with optimization_problem.py #
        self.simulation.restart_fields()
        self.simulation.clear_dft_monitors()
        self.simulation.change_sources(adjoint_sources)
        #                                            #
        adj_design_region_monitors = utils.install_design_region_monitors(
            self.simulation,
            self.design_regions,
            self.frequencies,
        )
        self.simulation.init_sim()
        sim_run_args = {
            "until_after_sources"
            if self.until_after_sources
            else "until": mp.stop_when_dft_decayed(
                self.dft_threshold, self.minimum_run_time, self.maximum_run_time
            )
        }
        self.simulation.run(**sim_run_args)

        return adj_design_region_monitors

    def _calculate_vjps(
        self,
        fwd_fields: List[List[mp.DftFields]],
        adj_fields: List[List[mp.DftFields]],
        design_variable_shapes: List[Tuple[int, ...]],
        sum_freq_partials: bool = True,
    ) -> List[onp.ndarray]:
        """Calculates the VJP for a given set of forward and adjoint fields."""
        return utils.calculate_vjps(
            self.simulation,
            self.design_regions,
            self.frequencies,
            fwd_fields,
            adj_fields,
            design_variable_shapes,
            sum_freq_partials=sum_freq_partials,
            finite_difference_step=self.finite_difference_step,
        )

    def _initialize_callable(self) -> Callable[[List[jnp.ndarray]], jnp.ndarray]:
        """Initializes the callable JAX function and registers its VJP."""

        @jax.custom_vjp
        def simulate(design_variables: List[jnp.ndarray]) -> jnp.ndarray:
            monitor_values, _ = self._run_fwd_simulation(design_variables)
            return monitor_values

        def _simulate_fwd(design_variables):
            """Runs forward simulation, returning monitor values and fields."""
            monitor_values, self.fwd_design_region_monitors = self._run_fwd_simulation(
                design_variables
            )
            design_variable_shapes = [x.shape for x in design_variables]
            return monitor_values, (design_variable_shapes)

        def _simulate_rev(res, monitor_values_grad):
            """Runs adjoint simulation, returning VJP of design wrt monitor values."""
            design_variable_shapes = res
            self.adj_design_region_monitors = self._run_adjoint_simulation(
                monitor_values_grad
            )
            vjps = self._calculate_vjps(
                self.fwd_design_region_monitors,
                self.adj_design_region_monitors,
                design_variable_shapes,
            )
            return ([jnp.asarray(vjp) for vjp in vjps],)

        simulate.defvjp(_simulate_fwd, _simulate_rev)

        return simulate

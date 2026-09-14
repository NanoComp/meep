"""Post-install check run by cibuildwheel against the repaired wheel.

Deliberately small: it confirms that the extension modules load, that the
vendored native libraries resolve, that a time step actually advances the
fields, and that the MPB-backed eigenmode machinery is present, the part
most likely to go missing if libmpb was not picked up at configure time.
"""

import sys

import numpy as np

import meep as mp

print(f"meep {mp.__version__} on {sys.platform}, python {sys.version.split()[0]}")

# The MPB bindings are a separate extension module (_mpb.so).
from meep import mpb  # noqa: E402

print(f"mpb bindings: {mpb.ModeSolver.__name__} available")

resolution = 20
cell = mp.Vector3(8, 8, 0)
geometry = [
    mp.Block(
        mp.Vector3(mp.inf, 1, mp.inf),
        center=mp.Vector3(),
        material=mp.Medium(epsilon=12),
    )
]
sources = [
    mp.Source(
        mp.ContinuousSource(frequency=0.15), component=mp.Ez, center=mp.Vector3(-3, 0)
    )
]

sim = mp.Simulation(
    cell_size=cell,
    boundary_layers=[mp.PML(1.0)],
    geometry=geometry,
    sources=sources,
    resolution=resolution,
)
sim.run(until=25)

ez = sim.get_array(center=mp.Vector3(), size=mp.Vector3(8, 8), component=mp.Ez)
# Allow a one-cell tolerance rather than pinning the exact grid size: whether a
# slice includes its far edge is a rounding detail of the array metadata, not
# something a packaging smoke test should assert.
assert ez.ndim == 2, ez.shape
assert all(abs(n - 8 * resolution) <= 1 for n in ez.shape), ez.shape
assert np.isfinite(ez).all(), "non-finite fields"
assert np.abs(ez).max() > 1e-6, "fields never got excited"

eps = sim.get_array(center=mp.Vector3(), size=mp.Vector3(8, 8), component=mp.Dielectric)
assert eps.shape == ez.shape, (eps.shape, ez.shape)
assert np.isclose(eps.max(), 12.0), eps.max()

print(f"max|Ez| = {np.abs(ez).max():.4g}, max(eps) = {eps.max():.4g}")

# Importing meep.mpb only proves the bindings are present.  Decomposing a field
# into eigenmodes is what actually calls into libmpb, so this is the check that
# would catch a wheel built without MPB or with libmpb left unvendored.
fcen = 0.15
mode_sim = mp.Simulation(
    cell_size=mp.Vector3(8, 4, 0),
    boundary_layers=[mp.PML(1.0)],
    geometry=geometry,
    sources=[
        mp.EigenModeSource(
            mp.GaussianSource(fcen, fwidth=0.1),
            eig_band=1,
            size=mp.Vector3(0, 4),
            center=mp.Vector3(-2, 0),
        )
    ],
    resolution=resolution,
)
monitor = mode_sim.add_mode_monitor(
    fcen, 0, 1, mp.ModeRegion(center=mp.Vector3(2, 0), size=mp.Vector3(0, 4))
)
mode_sim.run(until_after_sources=30)

alpha = mode_sim.get_eigenmode_coefficients(monitor, [1]).alpha
forward = abs(alpha[0, 0, 0])
assert np.isfinite(forward), alpha
assert forward > 0, "eigenmode decomposition returned no forward power"

print(f"|alpha+|^2 = {forward ** 2:.4g}")

print("smoke test passed")

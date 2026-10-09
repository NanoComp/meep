"""Post-install check run by cibuildwheel against the repaired wheel.

Confirms that the extension modules load, that the vendored libraries resolve,
that time stepping advances the fields, and that the MPB-backed eigenmode
machinery works, the part most likely to go missing at configure time.
"""

import ctypes
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

import meep as mp

print(f"meep {mp.__version__} on {sys.platform}, python {sys.version.split()[0]}")

# The wheel carries only an MPI build, so this must be true even in one process.
assert mp.with_mpi(), "this wheel was not built with MPI"
assert mp.count_processors() == 1, mp.count_processors()

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
assert ez.ndim == 2, ez.shape
assert all(abs(n - 8 * resolution) <= 1 for n in ez.shape), ez.shape
assert np.isfinite(ez).all(), "non-finite fields"
assert np.abs(ez).max() > 1e-6, "fields never got excited"

eps = sim.get_array(center=mp.Vector3(), size=mp.Vector3(8, 8), component=mp.Dielectric)
assert eps.shape == ez.shape, (eps.shape, ez.shape)
assert np.isclose(eps.max(), 12.0), eps.max()

print(f"max|Ez| = {np.abs(ez).max():.4g}, max(eps) = {eps.max():.4g}")

# Importing meep.mpb only proves the bindings are present; decomposing a field
# is what calls into libmpb, so this catches an unvendored or missing libmpb.
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

# Serial HDF5 would still work, writing one rank at a time through h5file.cpp's
# exclusive-access path, so nothing but this check would notice it.
package = Path(mp.__file__).parent
# delocate keeps them inside the package, auditwheel in a sibling <dist>.libs.
# Other distributions vendor an HDF5 of their own into site-packages, h5py
# above all, so look only where this wheel's copy can be.
lib_dirs = [package / ".dylibs"] + [
    d for d in package.parent.glob("*.libs") if d.name.startswith("meep")
]
vendored = [lib for d in lib_dirs for lib in sorted(d.glob("libhdf5*"))]
assert vendored, f"no vendored libhdf5 found in {[str(d) for d in lib_dirs]}"
for lib in vendored:
    assert hasattr(
        ctypes.CDLL(str(lib)), "H5Pset_fapl_mpio"
    ), f"{lib.name} has no MPI-IO support"
    print(f"hdf5: {lib.name} with MPI-IO")

# Running in one process proves nothing about the launcher: the failure worth
# catching is ranks that each get their own MPI_COMM_WORLD and silently run the
# whole simulation N times, which only shows up under a real mpiexec.
mpiexec = shutil.which("mpiexec")
assert mpiexec is not None, "no mpiexec on PATH; the mpich dependency is missing"

# Reported on stderr, not stdout: Meep points every non-master rank's stdout at
# /dev/null, so a print() would only ever be seen from rank 0.
probe = (
    "import sys, meep as mp; "
    "sys.stderr.write(f'rank {mp.my_rank()} of {mp.count_processors()}\\n')"
)
out = subprocess.run(
    [mpiexec, "-np", "2", sys.executable, "-m", "mpi4py", "-c", probe],
    capture_output=True,
    text=True,
)
detail = f"stdout: {out.stdout!r}\nstderr: {out.stderr!r}"
assert out.returncode == 0, f"parallel run failed:\n{detail}"
# An ABI mismatch between the wheel and this mpiexec gives each rank its own
# MPI_COMM_WORLD, which shows up as two ranks that both call themselves 0 of 1.
assert (
    "rank 0 of 2" in out.stderr and "rank 1 of 2" in out.stderr
), f"expected 2 ranks in one communicator\n{detail}"
print("2 ranks in one communicator")

print("smoke test passed")

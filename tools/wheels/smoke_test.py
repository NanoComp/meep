"""Post-install check run by cibuildwheel against the repaired wheel.

Deliberately small: it confirms that the extension modules load, that the
vendored native libraries resolve, that a time step actually advances the
fields, and that the MPB-backed eigenmode machinery is present, the part
most likely to go missing if libmpb was not picked up at configure time.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

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

# If this wheel carries the MPI build, check it actually runs in parallel.
# Loading it in-process is not enough: the failure worth catching is ranks that
# each get their own MPI_COMM_WORLD and silently run the whole simulation N
# times, which only shows up under a real launcher.
parallel_so = Path(mp.__file__).parent / "_parallel" / "_meep.so"
if not parallel_so.exists():
    print("no parallel build in this wheel; skipping the MPI check")
else:
    mpiexec = shutil.which("mpiexec")
    if mpiexec is None:
        print("parallel build present but no mpiexec on PATH; skipping")
    else:
        probe = "import meep as mp; print(mp.MEEP_PARALLEL, mp.count_processors())"
        out = subprocess.run(
            [mpiexec, "-np", "2", sys.executable, "-m", "mpi4py", "-c", probe],
            capture_output=True,
            text=True,
            env={**os.environ, "MEEP_MPI": "1"},
        )
        # stderr carries the RuntimeWarning naming *why* the parallel build was
        # skipped, which is the only useful thing when this fails in CI.
        detail = f"stdout: {out.stdout!r}\nstderr: {out.stderr!r}"
        assert out.returncode == 0, f"parallel run failed:\n{detail}"
        assert (
            "True 2" in out.stdout
        ), f"expected the parallel build with 2 ranks\n{detail}"
        print("parallel build: 2 ranks in one communicator")

        # The opposite failure: MPI unavailable under a multi-rank launcher.
        # Falling back to the serial build there would hand every rank the whole
        # simulation, so Meep must refuse instead of warning. Simulated by
        # putting a module that raises ImportError ahead of the real mpi4py.
        with tempfile.TemporaryDirectory() as blocked:
            Path(blocked, "mpi4py.py").write_text(
                'raise ImportError("mpi4py blocked by the smoke test")\n'
            )
            out = subprocess.run(
                [mpiexec, "-np", "2", sys.executable, "-c", "import meep"],
                capture_output=True,
                text=True,
                env={**os.environ, "PYTHONPATH": blocked},
            )
            detail = f"stdout: {out.stdout!r}\nstderr: {out.stderr!r}"
            assert out.returncode != 0, (
                "a 2-rank job without mpi4py silently fell back to the serial "
                f"build instead of failing\n{detail}"
            )
            assert (
                "refusing" in out.stderr
            ), f"expected Meep to refuse the silent serial fallback\n{detail}"
        print("serial fallback under a multi-rank launcher: refused")

print("smoke test passed")

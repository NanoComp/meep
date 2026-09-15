"""Select the serial or the MPI build of Meep's extension modules at import.

The wheel carries two compiled copies of Meep: a serial one at meep/_meep.so
and an MPI one at meep/_parallel/_meep.so.  HAVE_MPI is a compile-time choice,
so one binary cannot be both; shipping both and choosing here is what lets a
single `pip install pymeep` serve either.

The parallel copy keeps the file name _meep.so on purpose.  A CPython extension
module's init function is named after the file (PyInit__meep), so renaming it to
_meep_mpi.so would make it unimportable; it goes in a subdirectory instead and
is loaded from an explicit path.

This module is copied into the package by setup.py; it is not part of the
autotools build.
"""

import importlib.util
import os
import sys
import warnings
from pathlib import Path

# Set by the launcher, not by MPI itself, so these say how many ranks were
# *asked* for even when MPI_Init has quietly produced singletons.
_LAUNCHER_SIZE_VARS = (
    "OMPI_COMM_WORLD_SIZE",  # Open MPI
    "PMI_SIZE",  # MPICH, Intel MPI
    "MPI_LOCALNRANKS",
    "SLURM_NTASKS",
)

_MODULES = (
    ("meep._meep", "_parallel/_meep.so"),
    ("meep.mpb._mpb", "_parallel/mpb/_mpb.so"),
)


def _requested() -> bool:
    """Whether the caller wants the parallel build."""
    forced = os.environ.get("MEEP_MPI")
    if forced is not None:
        return forced not in ("", "0", "no", "false")

    # Being launched under a multi-rank launcher is the reliable signal: the
    # launcher sets these itself, in every rank, which is how the ranks find
    # each other in the first place.
    if _launcher_size() > 1:
        return True

    # Weaker fallback. `python -m mpi4py` does not always have MPI imported by
    # the time user code runs (measured False under `-m mpi4py -c`), so this
    # catches only the case where something else initialised MPI first.
    mpi = sys.modules.get("mpi4py.MPI")
    return bool(mpi is not None and mpi.Is_initialized())


def _launcher_size() -> int:
    for var in _LAUNCHER_SIZE_VARS:
        value = os.environ.get(var)
        if value and value.isdigit():
            return int(value)
    return 0


def _check_launcher_agrees(comm_size: int) -> None:
    """Catch the silent N-independent-jobs failure.

    The wheel is built against the MPICH ABI. Launched by an incompatible
    mpirun, typically Open MPI's, MPI_Init gets no usable handshake and
    every rank becomes its own MPI_COMM_WORLD. Nothing errors: you simply get N
    copies of the whole simulation, each believing it is rank 0. That is far
    worse than a crash, so fail loudly instead.
    """
    launched = _launcher_size()
    if launched > 1 and comm_size == 1:
        raise RuntimeError(
            f"Meep's parallel build sees MPI_COMM_WORLD size 1, but the job was "
            f"launched with {launched} processes. The launcher's MPI does not "
            f"match the one this wheel is built against (MPICH ABI, "
            f"libmpi.so.12 / libmpi.12.dylib), so every rank would run the "
            f"entire simulation "
            f"independently. Launch with an MPICH-ABI mpiexec (the `mpich` "
            f"wheel provides one), or build Meep from source against your MPI: "
            f"MEEP_CONFIGURE_ARGS='--with-mpi' pip install --no-binary pymeep pymeep"
        )


def _degrade(reason: str) -> bool:
    """Fall back to the serial build, or refuse when that would fan out.

    Under a multi-rank launcher the serial build does not merely mean "slower":
    every rank runs the whole simulation believing it is alone, which is the
    silent failure _check_launcher_agrees exists to catch. A warning is too
    quiet for that, so the fallback becomes an error instead.
    """
    launched = _launcher_size()
    if launched > 1:
        raise RuntimeError(
            f"The parallel build was requested but {reason} The job was "
            f"launched with {launched} processes, so using the serial build "
            f"would run the entire simulation in every one of them; refusing "
            f"rather than doing that silently."
        )
    warnings.warn(
        f"MEEP_MPI was requested but {reason} Using the serial build.",
        RuntimeWarning,
        stacklevel=3,
    )
    return False


def install() -> bool:
    """Point meep._meep / meep.mpb._mpb at the parallel build.

    Returns True if the parallel build was loaded. Falls back to serial with a
    warning when MPI is wanted but unusable, and never raises for a missing
    dependency, only for the dangerous mismatch above.
    """
    if not _requested():
        return False

    here = Path(__file__).parent
    missing = [rel for _, rel in _MODULES if not (here / rel).exists()]
    if missing:
        return _degrade(
            f"this wheel has no parallel build ({', '.join(missing)} absent)."
        )

    # The parallel extensions link libmpi (libmpi.so.12 on Linux,
    # @rpath/libmpi.12.dylib on macOS), which is deliberately not vendored: it
    # has to be the same MPI the launcher uses. Importing mpi4py.MPI first pulls
    # libmpi into the process: mpi4py dlopens the one beside it in the
    # environment, after which the loader satisfies our reference from what is
    # already mapped. On macOS this is the only thing that resolves it, since
    # delocate strips the build tree's rpath during repair.
    try:
        from mpi4py import MPI
    except ImportError:
        return _degrade(
            "mpi4py is not installed, so libmpi cannot be located. Install "
            "pymeep[mpi]."
        )

    try:
        for name, relpath in _MODULES:
            spec = importlib.util.spec_from_file_location(name, here / relpath)
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
    except ImportError as exc:
        for name, _ in _MODULES:
            sys.modules.pop(name, None)
        return _degrade(
            f"the parallel build could not be loaded ({exc}). This usually "
            "means no MPICH-ABI libmpi is available; Open MPI is not "
            "ABI-compatible with it."
        )

    _check_launcher_agrees(MPI.COMM_WORLD.Get_size())
    return True

"""Build the `meep` Python package as a binary wheel.

The canonical build system is autotools.  Rather than duplicating it, this
shim drives ./configure && make in an out-of-tree build directory and then
packages the `python/meep` directory that python/Makefile.am's `meep:` target
already assembles (the .py sources plus the SWIG modules _meep.so / _mpb.so).

Environment variables:
  MEEP_VERSION            override the version written into wheel metadata
  MEEP_DEPS_PREFIX        prefix where libctl/harminv/mpb were installed
  MEEP_CONFIGURE_ARGS     extra arguments appended to ./configure
  MEEP_BUILD_JOBS         parallelism for make (default: os.cpu_count())
"""

import hashlib
import os
import re
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Distribution, setup
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.egg_info import egg_info as _egg_info

HERE = Path(__file__).parent.resolve()

# Where setup.py stashes the shared libraries that _meep.so links against, so
# that auditwheel/delocate can find them at a stable path during repair.
#
# cibuildwheel expands {project} in before-all and test-command but NOT in
# repair-wheel-command, so the repair step cannot name a path relative to the
# project. An absolute path both steps agree on avoids the problem entirely.
WHEEL_LIBS = Path(os.environ.get("MEEP_WHEEL_LIBS") or HERE / "build" / "wheel-libs")


def meep_version() -> str:
    """Read the version out of configure.ac and normalize it for PEP 440."""
    override = os.environ.get("MEEP_VERSION")
    if override:
        return override

    text = (HERE / "configure.ac").read_text()
    match = re.search(
        r"AC_INIT\(\[meep\],\s*\[m4_esyscmd\(\./version\.sh\s+([^)\s]+)", text
    )
    if not match:
        raise RuntimeError("could not parse the Meep version out of configure.ac")

    raw = match.group(1)
    # configure.ac spells pre-releases as 1.35.0-beta; PEP 440 wants 1.35.0b0.
    suffixes = {"alpha": "a0", "beta": "b0", "rc": "rc0"}
    for name, pep440 in suffixes.items():
        if raw.endswith(f"-{name}"):
            return raw[: -len(name) - 1] + pep440
    return raw


def require(program: str, hint: str) -> None:
    if shutil.which(program) is None:
        raise SystemExit(
            f"error: {program!r} is required to build Meep from source but was not "
            f"found on PATH.\n{hint}"
        )


# Environment that changes what ./configure concludes without appearing in its
# argument list.
_CONFIGURE_ENV = (
    "CC",
    "CXX",
    "CFLAGS",
    "CXXFLAGS",
    "CPPFLAGS",
    "LDFLAGS",
    "LIBS",
    "PKG_CONFIG_PATH",
    "MEEP_DEPS_PREFIX",
)


def configure_fingerprint(args) -> str:
    """Identify a configuration, so that a changed one is not silently reused.

    ./configure is skipped when the build tree already has a Makefile, which is
    what makes a repeated build fast. On its own that would also ignore a
    changed MEEP_CONFIGURE_ARGS or dependency prefix and quietly produce a wheel
    built to the previous configuration; comparing this against a stamp file
    turns that case into a reconfigure.
    """
    material = [str(a) for a in args]
    material += [f"{key}={os.environ.get(key, '')}" for key in _CONFIGURE_ENV]
    return hashlib.sha256("\0".join(material).encode()).hexdigest()


def patchelf(*args) -> str:
    out = subprocess.run(
        ["patchelf", *[str(a) for a in args]],
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def read_soname(lib: Path) -> str:
    """The library's SONAME, or "" when it cannot be read.

    OSError covers platforms with no patchelf at all; the caller falls back to
    the file name, which is right for Mach-O where there is no SONAME anyway.
    """
    try:
        return patchelf("--print-soname", lib)
    except (subprocess.CalledProcessError, OSError):
        return ""


def strip_binaries(paths) -> None:
    """Strip our own shared objects before anything rewrites them.

    Running strip *after* patchelf is a known way to end up with a library the
    loader rejects ("ELF load command address/offset not properly aligned"), and
    auditwheel's --strip does exactly that. Stripping first keeps the size win
    without the hazard, so the repair step no longer passes --strip.
    """
    if not sys.platform.startswith("linux") or shutil.which("strip") is None:
        return
    for path in paths:
        if path.is_file() and not path.is_symlink():
            subprocess.run(["strip", "--strip-unneeded", str(path)], check=False)


def run(cmd, cwd, env=None) -> None:
    printable = " ".join(str(c) for c in cmd)
    print(f"[meep-build] (cd {cwd} && {printable})", flush=True)
    subprocess.run([str(c) for c in cmd], cwd=str(cwd), env=env, check=True)


class build_ext(_build_ext):
    """Run the autotools build, then copy python/meep into the wheel.

    This hangs off build_ext rather than the more obvious build_py because
    `build` only runs build_py when the distribution has pure modules, and
    `packages` is empty here.  build_ext runs because BinaryDistribution
    reports has_ext_modules(); the base run() is skipped since setup() never
    declares an Extension for it to compile.
    """

    def run(self):
        package_dir = self.build_meep()

        target = Path(self.build_lib) / "meep"
        if target.exists():
            shutil.rmtree(target)
        # copy2 keeps the executable bit on the .so files.
        shutil.copytree(package_dir, target, copy_function=shutil.copy2)
        self.stage_shared_libraries(package_dir.parent.parent)

        # Everything we ship, stripped before repair rather than during it.
        strip_binaries(list(target.rglob("*.so")) + list(WHEEL_LIBS.glob("*.so*")))

    # internals -------------------------------------------------------------

    def build_meep(self) -> Path:
        # The build is Python-ABI specific (Python.h, libpython), so cibuildwheel
        # reusing one container for several interpreters must not reuse one tree.
        tag = f"{sysconfig.get_platform()}-{sys.implementation.cache_tag}"
        builddir = HERE / "build" / f"autotools-{tag}"
        builddir.mkdir(parents=True, exist_ok=True)

        if not (HERE / "configure").exists():
            require("autoreconf", "Install autoconf, automake and libtool.")
            # autogen.sh runs autoreconf three times "just in case"; once is
            # enough here because we never re-run it against a dirty tree.
            run(
                ["autoreconf", "--verbose", "--install", "--symlink", "--force"],
                cwd=HERE,
            )

        configure_args = self.configure_args(builddir)
        fingerprint = configure_fingerprint(configure_args)
        stamp = builddir / ".meep-configure-stamp"
        configured = (
            (builddir / "Makefile").exists()
            and stamp.is_file()
            and stamp.read_text() == fingerprint
        )
        if not configured:
            require(
                "swig", "Install SWIG 4.x (needed to generate the Python bindings)."
            )
            run([HERE / "configure", *configure_args], cwd=builddir)
            stamp.write_text(fingerprint)

        jobs = os.environ.get("MEEP_BUILD_JOBS") or str(os.cpu_count() or 1)
        run(["make", f"-j{jobs}"], cwd=builddir)

        package_dir = builddir / "python" / "meep"
        if not (package_dir / "__init__.py").exists():
            # configure only warns when it turns the Python interface off, so
            # the build "succeeds" and quietly yields a wheel with no extension
            # modules.  Point at the usual cause instead.
            raise RuntimeError(
                f"the autotools build did not produce {package_dir}.\n"
                f"Check {builddir / 'config.log'}: if configure logged "
                '"disabling Python wrappers", its numpy or Python headers '
                "were not found."
            )
        return package_dir

    def configure_args(self, builddir: Path) -> list:
        args = [
            "--enable-maintainer-mode",  # regenerate the SWIG wrappers
            "--enable-shared",
            "--disable-static",
            "--without-scheme",  # no Guile inside a wheel
            "--without-mpi",  # wheels cannot portably ship an MPI runtime
            f"--prefix={builddir / 'install'}",
            f"PYTHON={sys.executable}",
        ]

        # src/ compiles identically for every interpreter, since only the SWIG
        # wrappers see Python headers, so across four Pythons and two variants
        # per job, ccache turns six of the eight src/ builds into hits. Sharing
        # one build tree does not achieve this: re-running configure regenerates
        # sphere-quad.h via BUILT_SOURCES and make rebuilds src/ regardless.
        if shutil.which("ccache"):
            args.append("--enable-ccache")

        prefix = os.environ.get("MEEP_DEPS_PREFIX")
        if prefix:
            libctl = Path(prefix) / "share" / "libctl"
            if libctl.is_dir():
                args.append(f"--with-libctl={libctl}")

        args += self.split_extra_args()
        return args

    @staticmethod
    def split_extra_args() -> list:
        import shlex

        return shlex.split(os.environ.get("MEEP_CONFIGURE_ARGS", ""))

    @staticmethod
    def stage_shared_libraries(builddir: Path) -> None:
        """Copy libmeep/libpympb somewhere auditwheel and delocate can see.

        Staged files are named by their SONAME, not by their on-disk name: a
        consumer's DT_NEEDED says `libmeep.so.38` while libtool's .libs holds
        that only as a symlink to libmeep.so.38.0.0, and symlinks are skipped
        here.  Naming by SONAME is what lets the loader resolve against this
        directory instead of relying on the RPATH into the build tree.

        """
        WHEEL_LIBS.mkdir(parents=True, exist_ok=True)
        for subdir in ("src/.libs", "libpympb/.libs"):
            source = builddir / subdir
            if not source.is_dir():
                continue
            for lib in source.iterdir():
                if lib.is_symlink() or not re.search(
                    r"\.(so|dylib)(\.\d+)*$", lib.name
                ):
                    continue
                name = read_soname(lib) or lib.name
                shutil.copy2(lib, WHEEL_LIBS / name)


class egg_info(_egg_info):
    """Record `meep` as the import name.

    setuptools derives top_level.txt from `packages` and `ext_modules`, both of
    which are empty here, so it would otherwise write an empty file and tools
    that map an import name back to a distribution would not find pymeep.
    """

    def run(self):
        super().run()
        Path(self.egg_info, "top_level.txt").write_text("meep\n")


class BinaryDistribution(Distribution):
    """Force a platform-specific, ABI-tagged wheel.

    setup() is not given any ext_modules (autotools builds them), so setuptools
    would otherwise tag the wheel py3-none-any.
    """

    def has_ext_modules(self) -> bool:
        return True

    def is_pure(self) -> bool:
        return False


setup(
    version=meep_version(),
    packages=[],  # build_py copies the autotools output instead
    distclass=BinaryDistribution,
    cmdclass={"build_ext": build_ext, "egg_info": egg_info},
    zip_safe=False,
)

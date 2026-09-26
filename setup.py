"""Build the `meep` Python package as a binary wheel.

Rather than duplicating the autotools build, this shim runs ./configure && make
out of tree and packages the `python/meep` directory that python/Makefile.am's
`meep:` target assembles (the .py sources plus _meep.so / _mpb.so).

The extensions are built against MPI (MPICH ABI).  libmpi is never vendored:
it has to be the one the launcher uses, so the wheel depends on mpi4py to load
it before _meep.so does.

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

# Where the shared libraries _meep.so links against are stashed for
# auditwheel/delocate. Absolute, because cibuildwheel does not expand {project}
# in repair-wheel-command.
WHEEL_LIBS = Path(os.environ.get("MEEP_WHEEL_LIBS") or HERE / "build" / "wheel-libs")


# configure.ac spells pre-releases as 1.35.0-beta; PEP 440 wants 1.35.0b0.
_RELEASE_TAGS = {"alpha": "a0", "beta": "b0", "rc": "rc0"}


def meep_version() -> str:
    """Read the version out of configure.ac and normalize it for PEP 440.

    version.sh takes the release tag either glued to the version or as a second
    argument, and every form has to be recognised: reading `1.35.0 alpha` as
    plain `1.35.0` would publish a pre-release under the final release's
    version, which PyPI then refuses to let anyone correct.
    """
    override = os.environ.get("MEEP_VERSION")
    if override:
        return override

    text = (HERE / "configure.ac").read_text()
    match = re.search(
        r"AC_INIT\(\[meep\],\s*\[m4_esyscmd\(\./version\.sh([^)]*)\)", text
    )
    if not match:
        raise RuntimeError("could not parse the Meep version out of configure.ac")

    args = match.group(1).split()
    if not 1 <= len(args) <= 2:
        raise RuntimeError(f"unexpected version.sh arguments in configure.ac: {args}")

    version, _, tag = args[0].partition("-")
    if len(args) == 2:
        if tag:
            raise RuntimeError(
                f"configure.ac gives the release tag twice: {' '.join(args)}"
            )
        tag = args[1]

    if not re.fullmatch(r"\d+(\.\d+)*", version):
        raise RuntimeError(f"not a release number in configure.ac: {version!r}")
    if not tag:
        return version
    if tag not in _RELEASE_TAGS:
        raise RuntimeError(
            f"unknown release tag {tag!r} in configure.ac; "
            f"expected one of {', '.join(sorted(_RELEASE_TAGS))}"
        )
    return version + _RELEASE_TAGS[tag]


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


MACHO = sys.platform == "darwin"


def patchelf(*args) -> str:
    out = subprocess.run(
        ["patchelf", *[str(a) for a in args]],
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def otool(*args) -> str:
    out = subprocess.run(
        ["otool", *[str(a) for a in args]],
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout


def read_install_id(lib: Path) -> str:
    """The Mach-O install name (LC_ID_DYLIB), or "" for a bundle that has none.

    `otool -D` prints the file name first, so a one-line answer means no id.
    """
    lines = otool("-D", lib).splitlines()
    return lines[1].strip() if len(lines) > 1 else ""


def read_soname(lib: Path) -> str:
    """The name a consumer uses to ask for this library, or "" if unreadable.

    ELF records a bare SONAME; Mach-O records a path whose basename is what
    delocate names the copy in .dylibs. The caller falls back to the file name.
    """
    try:
        if MACHO:
            return Path(read_install_id(lib)).name
        return patchelf("--print-soname", lib)
    except (subprocess.CalledProcessError, OSError):
        return ""


def strip_binaries(paths) -> None:
    """Strip our own shared objects before anything rewrites them.

    Running strip *after* patchelf is a known way to end up with a library the
    loader rejects ("ELF load command address/offset not properly aligned"), and
    auditwheel's --strip does exactly that. Stripping first keeps the size win
    without the hazard, which is why the repair step does not pass --strip.
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
        strip_binaries(list(target.rglob("*.so")))
        self.stage_shared_libraries(package_dir.parent.parent)
        self.install_mpi_preload(target)

    def install_mpi_preload(self, target: Path) -> None:
        """Make `import meep` load libmpi before _meep.so asks for it.

        _meep.so links libmpi (libmpi.so.12, @rpath/libmpi.12.dylib), which the
        repair step deliberately leaves unvendored, so nothing resolves it on
        its own: on macOS delocate strips the build tree's rpath, and on Linux
        site-packages is not a library search path. Importing mpi4py first
        dlopens the launcher's libmpi, and the loader then satisfies our
        reference from what is already mapped.

        meep/__init__.py is SWIG output with a version line appended by
        python/Makefile.am; prepending here keeps this out of the autotools
        build entirely.
        """
        init = target / "__init__.py"
        prelude = "from mpi4py import MPI as _MPI  # noqa: F401  (loads libmpi)\n"
        init.write_text(prelude + init.read_text(encoding="utf-8"), encoding="utf-8")

    def build_meep(self) -> Path:
        # The build is Python-ABI specific (Python.h, libpython), so cibuildwheel
        # reusing one container for several interpreters must not reuse one tree.
        tag = f"{sysconfig.get_platform()}-{sys.implementation.cache_tag}"
        builddir = HERE / "build" / f"autotools-{tag}"
        builddir.mkdir(parents=True, exist_ok=True)

        if not (HERE / "configure").exists():
            require("autoreconf", "Install autoconf, automake and libtool.")
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
            env = dict(os.environ)
            # --with-mpi needs MPI compiler wrappers, and the C and C++ halves
            # have to come from the same implementation, so these are set
            # together or not at all. A caller-chosen pair is left alone: on a
            # cluster it names the site MPI's wrappers.
            chosen = [var for var in ("CC", "CXX") if env.get(var)]
            if not chosen:
                for wrapper in ("mpicc", "mpicxx"):
                    require(wrapper, "Install MPICH (pip install mpich).")
                env["CC"], env["CXX"] = "mpicc", "mpicxx"
            elif len(chosen) == 1:
                missing = "CXX" if chosen[0] == "CC" else "CC"
                raise SystemExit(
                    f"error: {chosen[0]}={env[chosen[0]]!r} is set but {missing} "
                    f"is not. Meep is built --with-mpi, so set both to the "
                    f"wrappers of one MPI, or neither to use the mpicc and "
                    f"mpicxx on PATH."
                )
            run([HERE / "configure", *configure_args], cwd=builddir, env=env)
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
            # libmpi is never vendored: it must be the one the launcher uses.
            "--with-mpi",
            f"--prefix={builddir / 'install'}",
            f"PYTHON={sys.executable}",
        ]

        # Only the SWIG wrappers see Python headers, so src/ compiles
        # identically for every interpreter and ccache serves the repeats.
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

        Staged under the name consumers ask for (DT_NEEDED says libmeep.so.38,
        which .libs holds only as a symlink to libmeep.so.38.0.0), so the repair
        step resolves here instead of through an RPATH into the build tree.
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
                dest = WHEEL_LIBS / (read_soname(lib) or lib.name)
                shutil.copy2(lib, dest)
                strip_binaries([dest])


class egg_info(_egg_info):
    """Record `meep` as the import name.

    setuptools derives top_level.txt from `packages` and `ext_modules`, both of
    which are empty here, so it would otherwise write an empty file and tools
    that map an import name back to a distribution would not find meep.
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
    packages=[],  # build_ext copies the autotools output instead
    distclass=BinaryDistribution,
    cmdclass={"build_ext": build_ext, "egg_info": egg_info},
    zip_safe=False,
)

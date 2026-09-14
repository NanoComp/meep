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


def want_mpi() -> bool:
    """Whether to build the second, MPI-enabled copy of the extensions."""
    return os.environ.get("MEEP_BUILD_MPI", "") not in ("", "0", "no", "false")


# ELF and Mach-O solve the same problem, "which library is this, and what does
# it load?", with different load commands and different tools, so the handful
# of helpers below are the only places that care which one is underfoot.
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


def install_name_tool(*args) -> None:
    """Rewrite Mach-O load commands, then restore the ad-hoc signature.

    Editing a Mach-O file invalidates its code signature, and arm64 refuses to
    map an image whose signature does not match its contents. delocate re-signs
    whatever it rewrites during repair, but these edits happen before repair and
    for the staged libraries they are the last edit before the loader sees them.
    """
    subprocess.run(
        ["install_name_tool", *[str(a) for a in args]],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        ["codesign", "--force", "--sign", "-", str(args[-1])],
        check=False,
        capture_output=True,
    )


def read_install_id(lib: Path) -> str:
    """The Mach-O install name (LC_ID_DYLIB), or "" for a bundle that has none.

    `otool -D` echoes the file name first and then the id, so a one-line answer
    means there was no id to print, which is the case for the .so bundles that
    CPython loads as extension modules.
    """
    lines = otool("-D", lib).splitlines()
    return lines[1].strip() if len(lines) > 1 else ""


def read_soname(lib: Path) -> str:
    """The name a consumer uses to ask for this library, or "" if unreadable.

    ELF records a bare SONAME. Mach-O records an install name that is normally a
    full path, but only its last component survives into the wheel, which is
    what delocate names the copy it drops in .dylibs, so the basename is the
    Mach-O answer to the same question.

    OSError covers a platform with neither tool; the caller falls back to the
    file name.
    """
    try:
        if MACHO:
            return Path(read_install_id(lib)).name
        return patchelf("--print-soname", lib)
    except (subprocess.CalledProcessError, OSError):
        return ""


def suffixed(soname: str, suffix: str) -> str:
    """Insert `suffix` into a library name, before the version and extension.

    libmeep.so.38 -> libmeep_mpi.so.38, and libmeep.38.dylib -> libmeep_mpi.38
    .dylib. Keyed on the name rather than the platform because macOS uses .so
    for the extension-module bundles and .dylib for the libraries they link.
    """
    if ".dylib" in soname:
        stem, dot, rest = soname.partition(".")
    else:
        stem, dot, rest = soname.partition(".so")
    return f"{stem}{suffix}{dot}{rest}"


def set_soname(lib: Path, name: str) -> None:
    """Rename a staged library from the loader's point of view."""
    if not MACHO:
        patchelf("--set-soname", name, lib)
        return
    # A Mach-O id is a path into a build tree that will not exist at runtime and
    # that delocate overwrites during repair anyway; only the basename carries
    # meaning here, so keep the directory the build gave it and swap the name.
    old = read_install_id(lib)
    install_name_tool("-id", str(Path(old).with_name(name)) if old else name, lib)


def read_needed(consumer: Path) -> list:
    """The libraries a binary asks the loader for.

    ELF gives bare SONAMEs; Mach-O gives the full path recorded in each
    LC_LOAD_DYLIB. `otool -L` leads with the file's own install name, which is
    not a dependency.
    """
    if not MACHO:
        return patchelf("--print-needed", consumer).split()

    install_id = read_install_id(consumer)
    needed = []
    for line in otool("-L", consumer).splitlines()[1:]:
        path = line.split("(", 1)[0].strip()
        if path and path != install_id:
            needed.append(path)
    return needed


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


def retarget_needed(consumers, mapping: dict) -> None:
    """Repoint each consumer at the renamed MPI libraries.

    Both variants build a libmeep under the same name, and both loaders resolve
    by name, so without this the two extension modules would bind to whichever
    copy was found first.
    """
    for consumer in consumers:
        for needed in read_needed(consumer):
            new_name = mapping.get(Path(needed).name if MACHO else needed)
            if not new_name:
                continue
            if MACHO:
                # Keep the directory. It points into a build tree that is gone
                # by the time anything loads this, but so does the serial
                # build's: delocate resolves either one by basename through the
                # DYLD_LIBRARY_PATH the repair command sets.
                install_name_tool(
                    "-change", needed, str(Path(needed).with_name(new_name)), consumer
                )
            else:
                patchelf("--replace-needed", needed, new_name, consumer)


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
        package_dir = self.build_meep("serial")

        target = Path(self.build_lib) / "meep"
        if target.exists():
            shutil.rmtree(target)
        # copy2 keeps the executable bit on the .so files.
        shutil.copytree(package_dir, target, copy_function=shutil.copy2)
        # Stripped here, as each binary lands and before anything rewrites it.
        # The staged libraries and the MPI copies are stripped the same way at
        # their own copy sites, so nothing is stripped twice or after patching.
        strip_binaries(list(target.rglob("*.so")))
        self.stage_shared_libraries(package_dir.parent.parent)

        if want_mpi():
            self.add_parallel_build(target)

        self.install_dispatch(target)

    # internals -------------------------------------------------------------

    def add_parallel_build(self, target: Path) -> None:
        """Build a second, MPI-enabled copy and place it under meep/_parallel.

        Both builds produce a libmeep under the same name, and neither repair
        tool can hold two of those in one wheel: auditwheel resolves DT_NEEDED
        by name, and delocate copies into .dylibs by basename and refuses a
        collision outright. Renaming here, before repair, keeps the whole
        thing inside cibuildwheel's ordinary one-build-one-repair flow.
        """
        if not (sys.platform.startswith("linux") or MACHO):
            raise SystemExit(
                f"error: MEEP_BUILD_MPI is supported on Linux and macOS, not "
                f"{sys.platform}. Unset it, or build against your own MPI with "
                "MEEP_CONFIGURE_ARGS='--with-mpi'."
            )

        package_dir = self.build_meep("mpi")
        builddir = package_dir.parent.parent

        parallel = target / "_parallel"
        parallel.mkdir(parents=True, exist_ok=True)
        # The file name stays _meep.so: the init symbol is PyInit__meep, so the
        # module cannot simply be renamed. See tools/wheels/_dispatch.py.
        shutil.copy2(package_dir / "_meep.so", parallel / "_meep.so")
        mpb_src = package_dir / "mpb" / "_mpb.so"
        if mpb_src.exists():
            (parallel / "mpb").mkdir(exist_ok=True)
            shutil.copy2(mpb_src, parallel / "mpb" / "_mpb.so")

        # Before retarget_needed rewrites them, never after.
        strip_binaries(list(parallel.rglob("*.so")))

        if MACHO:
            require("install_name_tool", "Install the Xcode command line tools.")
        else:
            require("patchelf", "Install patchelf (present in the manylinux images).")
        renamed = self.stage_shared_libraries(builddir, suffix="_mpi")
        # The staged MPI libs depend on each other too (libpympb needs libmeep),
        # so they are consumers of the rename as much as the extensions are.
        consumers = list(parallel.rglob("*.so"))
        consumers += [WHEEL_LIBS / new for new in renamed.values()]
        retarget_needed(consumers, renamed)

    def install_dispatch(self, target: Path) -> None:
        """Add the import-time serial/MPI selector to the package."""
        shutil.copy2(
            HERE / "tools" / "wheels" / "_dispatch.py", target / "_dispatch.py"
        )

        # meep/__init__.py is SWIG output with a version line appended by
        # python/Makefile.am. Prepending here rather than changing that rule
        # keeps the MPI wheel logic out of the autotools build entirely.
        init = target / "__init__.py"
        prelude = (
            "from . import _dispatch as _meep_dispatch\n"
            "MEEP_PARALLEL = _meep_dispatch.install()\n"
        )
        init.write_text(prelude + init.read_text(encoding="utf-8"), encoding="utf-8")

    def build_meep(self, variant: str) -> Path:
        # The build is Python-ABI specific (Python.h, libpython), so cibuildwheel
        # reusing one container for several interpreters must not reuse one tree.
        tag = f"{sysconfig.get_platform()}-{sys.implementation.cache_tag}-{variant}"
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

        configure_args = self.configure_args(builddir, variant)
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
            env = None
            if variant == "mpi":
                # mpicc/mpicxx come from whatever MPI is on PATH; the `mpich`
                # wheel installs them, and so does a system MPICH.
                require("mpicc", "Install MPICH (pip install mpich) for the MPI build.")
                env = {**os.environ, "CC": "mpicc", "CXX": "mpicxx"}
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

    def configure_args(self, builddir: Path, variant: str = "serial") -> list:
        args = [
            "--enable-maintainer-mode",  # regenerate the SWIG wrappers
            "--enable-shared",
            "--disable-static",
            "--without-scheme",  # no Guile inside a wheel
            # libmpi is never vendored: it must be the one the launcher uses.
            "--with-mpi" if variant == "mpi" else "--without-mpi",
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
    def stage_shared_libraries(builddir: Path, suffix: str = "") -> dict:
        """Copy libmeep/libpympb somewhere auditwheel and delocate can see.

        Staged files are named the way consumers ask for them, not the way
        libtool left them on disk: a DT_NEEDED says `libmeep.so.38` while
        .libs holds that only as a symlink to libmeep.so.38.0.0, and symlinks
        are skipped here.  That naming is what lets the repair step resolve
        against this directory instead of an RPATH into the build tree.

        `suffix` distinguishes the MPI copies, whose names are otherwise
        identical to the serial ones.  Returns {old name: new name}.
        """
        WHEEL_LIBS.mkdir(parents=True, exist_ok=True)
        renamed = {}
        for subdir in ("src/.libs", "libpympb/.libs"):
            source = builddir / subdir
            if not source.is_dir():
                continue
            for lib in source.iterdir():
                if lib.is_symlink() or not re.search(
                    r"\.(so|dylib)(\.\d+)*$", lib.name
                ):
                    continue
                soname = read_soname(lib) or lib.name
                name = soname
                if suffix:
                    name = suffixed(soname, suffix)
                    renamed[soname] = name
                dest = WHEEL_LIBS / name
                shutil.copy2(lib, dest)
                # Strip first: set_soname rewrites the file, and stripping a
                # rewritten one is the hazard strip_binaries describes.
                strip_binaries([dest])
                if suffix:
                    set_soname(dest, name)
        return renamed


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

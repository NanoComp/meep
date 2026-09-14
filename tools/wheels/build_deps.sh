#!/usr/bin/env bash
#
# Build Meep's native dependencies for a binary wheel.
#
# Run as cibuildwheel's `before-all` step: it installs the C libraries that are
# packaged for the platform, then builds the three NanoComp siblings that are
# not (libctl, harminv, MPB) into $MEEP_DEPS_PREFIX.
#
# Everything is built shared, because auditwheel/delocate vendor whatever ends
# up linked into _meep.so.
#
# Nothing here needs Guile.  Meep's Scheme interface is disabled in the wheel,
# which is what lets MPB drop Guile via --without-libctl; libctl needs Guile to
# *generate* utils/geom-ctl-io.c, but the release tarball ships that file
# already generated, so --without-guile works from a tarball, and only from a
# tarball, since a Git checkout fails with "No rule to make target geom-ctl-io.c".
#
# Environment variables:
#   MEEP_DEPS_PREFIX        install prefix (default /usr/local)
#   MEEP_SKIP_SYSTEM_DEPS   do not touch the system package manager
#   MEEP_BUILD_MPI          also install MPICH, for the parallel wheel variant
#   WORKDIR                 where tarballs are unpacked

set -euo pipefail

PREFIX="${MEEP_DEPS_PREFIX:-/usr/local}"
LIBCTL_VERSION="${LIBCTL_VERSION:-4.7.1}"   # >= 4.7.0 for mesh geometry (configure.ac)
HARMINV_VERSION="${HARMINV_VERSION:-1.4.3}"
MPB_VERSION="${MPB_VERSION:-1.12.0}"
MPICH_VERSION="${MPICH_VERSION:-5.0.1.post1}"
PATCHELF_VERSION="${PATCHELF_VERSION:-0.19.1.0}"

WORKDIR="${WORKDIR:-$(mktemp -d)}"
NPROC="$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) )"

export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PREFIX}/lib64/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="${PREFIX}/lib:${PREFIX}/lib64:${LD_LIBRARY_PATH:-}"
export PATH="${PREFIX}/bin:${PATH}"

log() { printf '\n=== %s ===\n' "$*"; }

# gfortran is not optional: LAPACK is Fortran, and harminv's configure runs
# AC_F77_WRAPPERS to work out its name mangling.  CI runners happen to ship it
# preinstalled, which is why Meep's own workflow never names it.
install_system_deps_linux() {
  if [ -n "${MEEP_SKIP_SYSTEM_DEPS:-}" ]; then
    log "skipping system libraries (MEEP_SKIP_SYSTEM_DEPS set)"
    return
  fi

  if command -v dnf >/dev/null; then
    log "installing system libraries (dnf)"
    dnf -y install epel-release
    # CRB/PowerTools carries the -devel packages for several of these on EL8/9.
    dnf -y config-manager --set-enabled powertools 2>/dev/null \
      || dnf -y config-manager --set-enabled crb 2>/dev/null \
      || true
    dnf -y install \
      autoconf automake libtool pkgconfig swig ccache \
      gcc-gfortran \
      fftw-devel \
      gsl-devel \
      openblas-devel \
      lapack-devel \
      hdf5-devel \
      libpng-devel \
      zlib-devel
  elif command -v apt-get >/dev/null; then
    # Not a manylinux target: this branch is for building and testing the wheel
    # machinery on a plain Debian/Ubuntu box (including WSL).
    log "installing system libraries (apt-get)"
    apt-get -y update
    DEBIAN_FRONTEND=noninteractive apt-get -y install \
      build-essential \
      autoconf automake libtool pkg-config swig ccache \
      gfortran \
      libfftw3-dev \
      libgsl-dev \
      liblapack-dev \
      libhdf5-dev \
      libpng-dev \
      zlib1g-dev \
      python3-dev
  else
    echo "no supported package manager found; set MEEP_SKIP_SYSTEM_DEPS=1 and install by hand" >&2
    exit 1
  fi
}

install_system_deps_macos() {
  if [ -n "${MEEP_SKIP_SYSTEM_DEPS:-}" ]; then
    log "skipping system libraries (MEEP_SKIP_SYSTEM_DEPS set)"
    return
  fi
  log "installing system libraries (brew)"
  # `gcc` is what provides gfortran on macOS.
  brew install autoconf automake libtool pkg-config swig ccache gcc fftw gsl hdf5 libpng
}

# fetch_and_build <repo> <version> [configure args...]
#
# Uses the release tarball rather than a Git checkout: it carries a pre-built
# configure (so no autoreconf) and the generated ctl-io sources that would
# otherwise need Guile.
fetch_and_build() {
  local repo="$1" version="$2"; shift 2
  local name="${repo}-${version}"
  local url="https://github.com/NanoComp/${repo}/releases/download/v${version}/${name}.tar.gz"

  log "building ${name}"
  mkdir -p "${WORKDIR}"
  if [ ! -d "${WORKDIR}/${name}" ]; then
    curl -sSL --fail -o "${WORKDIR}/${name}.tar.gz" "${url}"
    tar xzf "${WORKDIR}/${name}.tar.gz" -C "${WORKDIR}"
  fi

  pushd "${WORKDIR}/${name}" >/dev/null
  ./configure --prefix="${PREFIX}" --enable-shared --disable-static "$@"
  make -j"${NPROC}"
  make install
  popd >/dev/null
}

# Install MPICH from its PyPI wheel, which carries a complete toolchain:
# bin/mpicc, bin/mpicxx, bin/mpiexec, include/mpi.h and lib/libmpi.so.12.
#
# Deliberately the same MPICH end users get from `pip install mpich`, so the
# wheel is built against exactly the runtime it will be paired with. libmpi is
# never vendored into the wheel: it has to be the launcher's own, which is
# why the MPICH ABI matters: libmpi.so.12 is the ABI tag.
# Any interpreter will do, since only the C toolchain is wanted and not the
# Python package, but it must have pip. The manylinux images' /usr/bin/python3 does
# not, so the /opt/python interpreters come first.
find_pip_python() {
  local candidate
  for candidate in /opt/python/cp31*/bin/python python3 python; do
    candidate=$(command -v "${candidate}" 2>/dev/null || echo "${candidate}")
    if [ -x "${candidate}" ] && "${candidate}" -m pip --version >/dev/null 2>&1; then
      echo "${candidate}"
      return 0
    fi
  done
  echo "no python with pip found" >&2
  return 1
}

# manylinux_2_28 pins patchelf 0.17.2, which miscomputes segment alignment when
# it has to grow a binary: renaming a SONAME there yields a library the loader
# rejects with "ELF load command address/offset not properly aligned". Fixed
# upstream in 0.18. auditwheel shells out to whatever patchelf is on PATH, so
# installing a current one ahead of the system copy fixes both it and setup.py.
install_patchelf() {
  local py
  py=$(find_pip_python) || exit 1
  log "installing patchelf ${PATCHELF_VERSION} into ${PREFIX} (image ships $(patchelf --version 2>&1))"
  "${py}" -m pip install --quiet --prefix="${PREFIX}" "patchelf==${PATCHELF_VERSION}"
  hash -r
  patchelf --version
}

install_mpich() {
  local py
  py=$(find_pip_python) || exit 1

  log "installing MPICH from the PyPI wheel into ${PREFIX} (using ${py})"
  "${py}" -m pip install --quiet --prefix="${PREFIX}" "mpich==${MPICH_VERSION}"

  if [ ! -x "${PREFIX}/bin/mpicc" ]; then
    echo "the mpich wheel did not provide ${PREFIX}/bin/mpicc" >&2
    exit 1
  fi
  "${PREFIX}/bin/mpicc" -show || true
}

main() {
  case "$(uname -s)" in
    Linux)  install_system_deps_linux ;;
    Darwin) install_system_deps_macos ;;
    *) echo "unsupported platform: $(uname -s)" >&2; exit 1 ;;
  esac

  # Before anything is patched: both setup.py and auditwheel need a good one.
  if [ "$(uname -s)" = Linux ]; then
    install_patchelf
  fi

  # Meep's Python build needs only libctlgeom from libctl; the Scheme half is
  # gated on --with-scheme in Meep's configure.ac.
  fetch_and_build libctl "${LIBCTL_VERSION}" --without-guile
  fetch_and_build harminv "${HARMINV_VERSION}"

  # --without-libctl drops Guile and the mpb executable and leaves the libmpb
  # C library, which is all Meep links against (AC_CHECK_LIB(mpb, ...)).
  fetch_and_build mpb "${MPB_VERSION}" --without-libctl --with-hermitian-eps LIBS=-ldl

  # Only needed when setup.py is asked to build the second, MPI-enabled copy.
  if [ -n "${MEEP_BUILD_MPI:-}" ] && [ "${MEEP_BUILD_MPI}" != "0" ]; then
    install_mpich
  fi

  # `|| true` because there is no lib64 on every platform, and a failing ls
  # under `set -e` would take the whole cibuildwheel before-all step down.
  log "dependency prefix contents"
  ls -1 "${PREFIX}/lib" "${PREFIX}/lib64" 2>/dev/null | sort -u | head -40 || true
}

main "$@"

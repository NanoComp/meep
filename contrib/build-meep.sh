#!/bin/bash

# latest version of this script can be found at:
#   https://github.com/NanoComp/meep/blob/master/contrib/build-meep.sh

# detect if script is in a directory called src/
DESTDIR=$(pwd)
[ ${DESTDIR##*/} = src ] && DESTDIR=$(cd $(pwd)/..; pwd)
SRCDIR=${DESTDIR}/src

cat <<EOF

This sript will download or update sources, compile and install MEEP.
Please ensure the following final paths fit your needs:
    '${DESTDIR}/bin/meep'
    '${DESTDIR}/share/...'
    '${DESTDIR}/...'
    '${SRCDIR}/<sources>'

Press return to continue
EOF
read junk

set -ex

ubuntu=false
centos=false

distrib=$(lsb_release -r -s)
case "$distrib" in
    18.04) # ubuntu 18.04 bionic
        libpng=libpng-dev
        libpython=libpython3-dev
        ubuntu=true
        ;;
    16.04) # ubuntu 16.04 xenal
        libpng=libpng16-dev
        libpython=libpython3.5-dev
        ubuntu=true
        ;;
    7.*) # CentOS 7.x
        centos=true
        ;;
    *)
        echo "unsupported distribution '$(lsb_release -a)', edit and fix!"
        false
        ;;
esac

# sudo not needed for 'make install'
#SUDO=sudo
SUDO=

mkdir -p ${SRCDIR}
cd ${SRCDIR}

RPATH_FLAGS="-Wl,-rpath,${DESTDIR}/lib:/usr/lib/x86_64-linux-gnu/hdf5/openmpi"
LDFLAGS="-L${DESTDIR}/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi ${RPATH_FLAGS}"
CFLAGS="-I${DESTDIR}/include -I/usr/include/hdf5/openmpi"
CPPFLAGS=${CFLAGS}
PKG_CONFIG_PATH=${DESDTIR}/pkgconfig
export PKG_CONFIG_PATH
export PATH=${DESTDIR}/bin:${PATH}

gitclone ()
{
    repo=${1##*/}
    name=${repo%%.*}
    echo $repo $name
    if [ -d $name ]; then
        ( cd $name; git pull; )
    else
        git clone --depth=1 $1
    fi
}

autogensh ()
{
    sh autogen.sh PKG_CONFIG_PATH="${PKG_CONFIG_PATH}" RPATH_FLAGS="${RPATH_FLAGS}" LDFLAGS="${LDFLAGS}" CFLAGS="${CFLAGS}" CPPFLAGS="${CPPFLAGS}" \
        --disable-static --enable-shared --prefix="${DESTDIR}" \
        --with-libctl=${DESTDIR}/share/libctl \
        "$@"
}

if $ubuntu; then

    sudo apt-get update

    # If building on Ubuntu 18.04LTS, replace libpng16-dev with libpng-dev,
    # and libpython3.5-dev with libpython3-dev.
    sudo apt-get -y install     \
        build-essential         \
        gfortran                \
        libblas-dev             \
        liblapack-dev           \
        libgmp-dev              \
        swig                    \
        libgsl-dev              \
        autoconf                \
        pkg-config              \
        $libpng                 \
        git                     \
        guile-2.0-dev           \
        libfftw3-dev            \
        libhdf5-openmpi-dev     \
        hdf5-tools              \
        $libpython              \
        python3-numpy           \
        python3-scipy           \
        python3-pip             \
        ffmpeg                  \

    # The next line is only required on Ubuntu  16.04
    [ "$distrib" = xenial ] && sudo -H pip3 install --upgrade pip

    sudo -H pip3 install --no-cache-dir mpi4py
    export HDF5_MPI="ON"
    sudo -H pip3 install --no-binary=h5py h5py
    sudo -H pip3 install matplotlib>3.0.0
fi

if $centos; then

    sudo yum -y --enablerepo=extras install epel-release

    sudo yum -y install   \
        bison             \
        byacc             \
        cscope            \
        ctags             \
        cvs               \
        diffstat          \
        oxygen            \
        flex              \
        gcc               \
        gcc-c++           \
        gcc-gfortran      \
        gettext           \
        git               \
        indent            \
        intltool          \
        libtool           \
        patch             \
        patchutils        \
        rcs               \
        redhat-rpm-config \
        rpm-build         \
        subversion        \
        systemtap         \
        wget

    sudo yum -y install    \
        openblas-devel     \
        fftw3-devel        \
        libpng-devel       \
        gsl-devel          \
        gmp-devel          \
        pcre-devel         \
        libtool-ltdl-devel \
        libunistring-devel \
        libffi-devel       \
        gc-devel           \
        zlib-devel         \
        openssl-devel      \
        sqlite-devel       \
        bzip2-devel        \
        ffmpeg

    sudo yum -y install    \
        openmpi-devel      \
        hdf5-openmpi-devel \
        guile-devel

    export PATH=${PATH}:/usr/lib64/openmpi/bin
fi

mkdir -p $SRCDIR

cd $SRCDIR
gitclone https://github.com/NanoComp/harminv.git
cd harminv/
autogensh
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/NanoComp/libctl.git
cd libctl/
autogensh
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/NanoComp/h5utils.git
cd h5utils/
autogensh CC=mpicc
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/NanoComp/mpb.git
cd mpb/
autogensh CC=mpicc --with-hermitian-eps
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/HomerReid/libGDSII.git
cd libGDSII/
autogensh
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/NanoComp/meep.git
cd meep/
autogensh --with-mpi --with-openmp PYTHON=python3
make -j && $SUDO make install

# all done

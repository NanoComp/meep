#!/bin/bash

# Latest version of this script can be found at:
#   https://github.com/NanoComp/meep/blob/master/contrib/build-meep.sh

help ()
{
    cat << EOF

$1: Download MEEP sources and dependencies, compile, and install

Usage: $1 [options]
EOF
    sed -ne 's,[ \t]*\(-[^ \t]*\))[^#]*#[ \t]*\(.*\),    \1 \2,p' "$1"
    echo ""
    exit 1
}

[ -z "$1" ] && echo "(use -h for help)"

installdeps=true
with_guile=true
with_parallel=true

while [ ! -z "$1" ]; do
    case "$1" in
        -h)         # help
            help "$0"
            ;;
        -d)         # <installdir>  (default: current directory)
            DESTDIR="$2"
            shift
            ;;
        -s)         # use 'sudo' for 'make install'
            SUDO=sudo
            ;;
        -n)         # do not check for distribution dependencies
            installdeps=false
            ;;
        -g)        # build without guile
            with_guile=false
            ;;
        -se)        # build serial version
            with_parallel=false
            ;;
        *)
            echo "'$1' ?"
            help "$0"
            ;;
    esac
    shift
done


# detect wether DESTDIR is ending with src/
[ -z ${DESTDIR} ] && DESTDIR=$(pwd)
[ ${DESTDIR##*/} = src ] && DESTDIR=$(cd $(pwd)/..; pwd)
SRCDIR=${DESTDIR}/src

cat << EOF

This sript will download or update sources, compile and install MEEP.
Please ensure the following final paths fit your needs:
    '${DESTDIR}/bin/meep'
    '${DESTDIR}/lib/...'
    '${DESTDIR}/share/...'
    '${DESTDIR}/...'
    '${SRCDIR}/<sources>'

Press return to continue
EOF
#read junk

if ! lsb_release; then
    echo "Minimum requirements:"
    echo "    Ubuntu:"
    echo "        sudo apt-get -y install lsb-release sudo git"
    echo "    CentOS:"
    echo "        sudo yum -y install redhat-lsb-core sudo git"
    echo ""
    exit 1
fi

set -ex


mkdir -p ${SRCDIR}
cd ${SRCDIR}

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
        "$@"
}

configure()
{
    ./configure PKG_CONFIG_PATH="${PKG_CONFIG_PATH}" RPATH_FLAGS="${RPATH_FLAGS}" LDFLAGS="${LDFLAGS}" CFLAGS="${CFLAGS}" CPPFLAGS="${CPPFLAGS}" \
        --disable-static --enable-shared --prefix="${DESTDIR}" \
        "$@"
}

getlibctl()
{
    if $with_guile; then
        gitclone https://github.com/NanoComp/libctl.git
    else
        wget https://github.com/NanoComp/libctl/releases/download/v4.5.1/libctl-4.5.1.tar.gz
        tar -xzf libctl-4.5.1.tar.gz
        mv libctl-4.5.1 libctl
        rm libctl-4.5.1.tar.gz
        echo $PWD
    fi
}

if $installdeps; then

    ubuntu=false
    centos=false

    distrib=$(lsb_release -r -s)
    case "$distrib" in
        18.04) # ubuntu 18.04 bionic
            libpng=libpng-dev
            libpython=libpython3-dev
            ubuntu=true
            ;;
        16.04) # ubuntu 16.04 xenial
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

    if $ubuntu; then

        sudo apt-get update

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

        [ "$distrib" = 16.04 ] && sudo -H pip3 install --upgrade pip
        sudo -H pip3 install --no-cache-dir mpi4py
        export HDF5_MPI="ON"
        sudo -H pip3 install --no-binary=h5py h5py
        sudo -H pip3 install matplotlib>3.0.0

        RPATH_FLAGS="-Wl,-rpath,${DESTDIR}/lib:/usr/lib/x86_64-linux-gnu/hdf5/openmpi"
        LDFLAGS="-L${DESTDIR}/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi ${RPATH_FLAGS}"
        CFLAGS="-I${DESTDIR}/include -I/usr/include/hdf5/openmpi"

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
            guile-devel        \
            swig

    export PATH=${PATH}:/usr/lib64/openmpi/bin
    RPATH_FLAGS="-Wl,-rpath,${DESTDIR}/lib:/usr/lib64/openmpi/lib"
    LDFLAGS="-L${DESTDIR}/lib -L/usr/lib64/openmpi/lib ${RPATH_FLAGS}"
    CFLAGS="-I${DESTDIR}/include -I/usr/include/openmpi-x86_64/"
    fi
else
    RPATH_FLAGS="-Wl,-rpath,-L${DESTDIR}/lib"
    LDFLAGS="-L${DESTDIR}/lib"
    CFLAGS="-I${DESTDIR}/include"
fi

CPPFLAGS=${CFLAGS}
PKG_CONFIG_PATH=${DESTDIR}/pkgconfig
PATH=${DESTDIR}/bin:${PATH}


if $with_guile; then
    if_guile_flag=""
    guile_flag_mpb=""
    guile_flag_meep=""
else
    if_guile_flag="--without-guile"
    guile_flag_mpb="--without-libctl"
    guile_flag_meep="--without-scheme"
fi

if $with_parallel; then
    parallel_flags="--with-mpi --with-openmp"
else
    parallel_flags=""
fi

mkdir -p $SRCDIR


cd $SRCDIR
gitclone https://github.com/NanoComp/harminv.git
cd harminv/
autogensh ${CONFIGURE_harminv}
make -j && $SUDO make install

cd $SRCDIR
getlibctl
cd libctl
if $with_guile;then
   autogensh ${CONFIGURE_libctl}
else
   configure ${if_guile_flag} ${CONFIGURE_libctl}
fi
make -j && $SUDO make install


cd $SRCDIR
gitclone https://github.com/NanoComp/h5utils.git
cd h5utils/
autogensh ${if_guile_flag} ${CONFIGURE_h5utils}
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/NanoComp/mpb.git
cd mpb/
autogensh --with-hermitian-eps --with-libctl=${DESTDIR}/share/libctl ${guile_flag_mpb} ${CONFIGURE_mpb}
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/HomerReid/libGDSII.git
cd libGDSII/
autogensh ${CONFIGURE_libGDSII}
make -j && $SUDO make install

cd $SRCDIR
gitclone https://github.com/NanoComp/meep.git
cd meep/
autogensh ${parallel_flags} PYTHON=python3 --with-libctl=${DESTDIR}/share/libctl ${guile_flag_meep} ${CONFIGURE_meep}
make -j && $SUDO make install

# all done

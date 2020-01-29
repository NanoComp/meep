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

gitclone ()
{
    repo=${1##*/}
    name=${repo%%.*}
    echo $repo $name
    if [ -d $name ]; then
        ( cd $name; git pull; )
    else
        [ -z "$2" ] || branch="-b $2"
        git clone --depth=1 $1 $branch
    fi
}

autogensh ()
{
    LIB64="${DESTDIR}/lib"
    $centos && LIB64="${DESTDIR}/lib64"
    LLP="${LD_LIBRARY_PATH}:${LIB64}"
    sh autogen.sh PKG_CONFIG_PATH="${PKG_CONFIG_PATH}" RPATH_FLAGS="${RPATH_FLAGS}" LDFLAGS="${LDFLAGS}" CFLAGS="${CFLAGS}" CPPFLAGS="${CPPFLAGS}" LD_LIBRARY_PATH=${LLP} \
        --disable-static --enable-shared --prefix="${DESTDIR}" --libdir=${LIB64} \
        --with-libctl=${DESTDIR}/share/libctl \
        "$@"
}

compile=true
installdeps=true
unset DESTDIR

while [ ! -z "$1" ]; do
    case "$1" in
        -h)         # help
            help "$0"
            ;;
        -d)         # <installdir>  (mandatory)
            DESTDIR="$2"
            shift
            ;;
        -s)         # use 'sudo' for 'make install'
            SUDO=sudo
            ;;
        -S)         # source directory (default: <installdir>/src)
            SRCDIR="$2"
            shift
            ;;
        -n)         # skip checking for distribution dependencies
            installdeps=false
            ;;
        -c)         # skip build+install
            compile=false
            ;;
        *)
            echo "'$1' ?"
            help "$0"
            ;;
    esac
    shift
done


[ -z "${DESTDIR}" ] && { echo "-d option is missing" ; help "$0"; }
# detect wether DESTDIR is ending with src/
[ "${DESTDIR##*/}" = src ] && DESTDIR=$(cd $(pwd)/..; pwd)
[ -z "$SRCDIR" ] && SRCDIR=${DESTDIR}/src

if $ubuntu; then
    RPATH_FLAGS="-Wl,-rpath,${DESTDIR}/lib:/usr/lib/x86_64-linux-gnu/hdf5/openmpi"
    LDFLAGS="-L${DESTDIR}/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi ${RPATH_FLAGS}"
    CFLAGS="-I${DESTDIR}/include -I/usr/include/hdf5/openmpi"
fi

if $centos; then
    export CC=/usr/lib64/openmpi/bin/mpicc
    export CXX=/usr/lib64/openmpi/bin/mpicxx
    export PATH=${PATH}:/usr/lib64/openmpi/bin
    RPATH_FLAGS="-Wl,-rpath,${DESTDIR}/lib64:/usr/lib64/openmpi/lib"
    LDFLAGS="-L${DESTDIR}/lib64 -L/usr/lib64/openmpi/lib ${RPATH_FLAGS}"
    CFLAGS="-I${DESTDIR}/include -I/usr/include/openmpi-x86_64/"

    showenv()
    {
        echo export PYTHONPATH+=:${DESTDIR}/lib/python3.6/site-packages
        echo export PYTHONPATH+=:/usr/local/lib64/python3.6/site-packages
        echo export LD_LIBRARY_PATH+=:${DESTDIR}/lib64
        echo export LD_PRELOAD+=:/usr/lib64/openmpi/lib/libmpi.so
    }
    eval $(showenv)
fi


if [ ! -r /etc/os-release ]; then
    echo "Error: cannot read /etc/os-release"
    false
fi

set -e

ubuntu=false
centos=false

. /etc/os-release
distrib="${ID}${VERSION_ID}"
case "$distrib" in
    ubuntu18.04) # ubuntu 18.04 bionic
        libpng=libpng-dev
        libpython=libpython3-dev
        ubuntu=true
        ;;
    ubuntu16.04) # ubuntu 16.04 xenial
        libpng=libpng16-dev
        libpython=libpython3.5-dev
        ubuntu=true
        ;;
    centos7) # CentOS 7.x
        centos=true
        ;;
    *)
        echo "Error: unsupported distribution '${distrib}'"
        false
        ;;
esac

if $installdeps && $ubuntu; then

    sudo apt-get update

    sudo apt-get -y install     \
        git                     \
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
    sudo -H pip3 install matplotlib\>3.0.0
fi

if $installdeps && $centos; then

    sudo yum -y --enablerepo=extras install epel-release

    sudo yum -y install   \
        git               \
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
        python3            \
        python3-devel      \
        python36-numpy      \
        python36-scipy

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

    sudo -E pip3 install mpi4py
fi

CPPFLAGS=${CFLAGS}
PKG_CONFIG_PATH=${DESDTIR}/pkgconfig
export PKG_CONFIG_PATH
export PATH=${DESTDIR}/bin:${PATH}

if $compile; then

mkdir -p ${SRCDIR}

cd ${SRCDIR}
gitclone https://github.com/NanoComp/harminv.git
cd harminv/
autogensh
make -j && $SUDO make install

cd ${SRCDIR}
gitclone https://github.com/NanoComp/libctl.git
cd libctl/
autogensh
make -j && $SUDO make install

cd ${SRCDIR}
gitclone https://github.com/NanoComp/h5utils.git
cd h5utils/
autogensh CC=mpicc
make -j && $SUDO make install

cd ${SRCDIR}
gitclone https://github.com/NanoComp/mpb.git
cd mpb/
autogensh CC=mpicc --with-hermitian-eps
make -j && $SUDO make install

cd ${SRCDIR}
gitclone https://github.com/HomerReid/libGDSII.git
cd libGDSII/
autogensh
make -j && $SUDO make install

cd ${SRCDIR}
#gitclone https://github.com/NanoComp/meep.git 
gitclone https://github.com/d-a-v/meep.git fixCentos7Installer
cd meep/
autogensh --with-mpi --with-openmp PYTHON=python3
make -j && $SUDO make install

# all done

if $centos; then
     cd ${DESTDIR}/lib/python3.6/site-packages/meep/
     for i in ../../../../lib64/python3.6/site-packages/meep/*meep*; do
        ln -sf $i
     done
fi

fi # install

########
# test

test=/tmp/test-meep.py

cat << EOF > $test
import meep as mp
cell = mp.Vector3(16,8,0)
exit()
EOF


echo "------------ ENV (commands)"
showenv
echo "------------ ENV (result)"
echo export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
echo export PYTHONPATH=${PYTHONPATH}
echo export LD_PRELOAD=${LD_PRELOAD}
echo "------------ $test"
cat $test
echo "------------ EXEC python3 $test"
python3 $test

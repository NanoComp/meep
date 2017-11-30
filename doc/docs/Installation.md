---
# Installation
---

**Note**: Installing Meep from source can be challenging for novice users. As a simple workaround, the latest version of Meep preinstalled on Ubuntu can be accessed on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetus.com/launchsims.html).

The main effort in installing Meep lies in installing the various prerequisite packages. This requires some understanding of how to install software on Unix systems.

It is also possible to install Meep on Windows systems by using a free Unix-compatibility environment such as [Cygwin](http://www.cygwin.org/). For more information, see these [step-by-step instructions](http://novelresearch.weebly.com/installing-meep-in-windows-8-via-cygwin.html).

For those installing Meep on a parallel supercomputer, a note of caution: most supercomputers have multiple compilers installed, and different versions of libraries compiled with different compilers. Meep is written in C++, and it is almost impossible to mix C++ code compiled by different compilers &mdash; pick one set of compilers by one vendor and stick with it consistently.

[TOC]

Installation on Linux
-------------------------

For most [Linux distributions](https://en.wikipedia.org/wiki/Linux_distribution), there should be precompiled packages for most of Meep's prerequisites below, and we *highly* recommend installing those prerequisites using the available packages for your system whenever possible. Using precompiled packages means that you don't have to worry about how to install things manually. You are using packages which have already been tweaked to work well with your system, and usually your packages will be automatically upgraded when you upgrade the rest of your system.

The following precompiled packages are available: BLAS and LAPACK possibly as part of a package for [Atlas BLAS](https://en.wikipedia.org/wiki/Automatically_Tuned_Linear_Algebra_Software), Guile, MPI, and HDF5. One thing to be careful of is that many distributions split packages into two parts: one main package for the libraries and programs, and a **devel** package for [header files](https://en.wikipedia.org/wiki/Header_file) and other things needed to compile software using those libraries. You will need to install **both**. So, for example, you will probably need both a `guile` package (probably installed by default) and a `guile-dev` or `guile-devel` package (probably *not* installed by default), and similarly for HDF5 etcetera. You will probably also want to install a `libpng-dev` or `libpng-devel` package in order to compile the `h5topng` utility in [h5utils](https://github.com/stevengj/h5utils/blob/master/README.md).

The easiest installation is on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu_(operating_system)) which has precompiled packages for Meep:

```sh
apt-get install meep h5utils
```

Zipped tarballs of stable versions of the source are available on the [releases page](https://github.com/stevengj/meep/releases).

To build the latest version of Meep from source on Ubuntu, follow these [instructions](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05819.html).

Installation on MacOS 
-----------------------

Since [MacOS](https://en.wikipedia.org/wiki/MacOS) is, at its heart, a Unix system, one can, in principle compile and install Meep and all its prerequisites just as on any other Unix system. However, this process is much easier using the [Homebrew](https://en.wikipedia.org/wiki/Homebrew_(package_management_software)) package to install most of the prerequisites, since it will handle dependencies and other details for you. You will need [administrator privileges](http://support.apple.com/kb/PH3920) on your Mac.

The first steps are:

-   Install [Xcode](https://en.wikipedia.org/wiki/Xcode), the development/compiler package from Apple, free from the [Apple Xcode web page](https://developer.apple.com/xcode/).
-   Install Homebrew: download from the [Homebrew site](http://brew.sh/) and follow the instructions there.
-   Run the following commands in the terminal to compile and install the prerequisites. This may take a while to complete because it will install lots of other stuff first

```sh
brew doctor
brew install homebrew/science/hdf5 homebrew/science/openblas guile fftw h5utils
```

Now, install the Harminv, libctl, MPB, and Meep packages from source. Download [Harminv](https://github.com/stevengj/harminv/blob/master/README.md) and, in the `harminv` directory, do:

```sh
./configure && make && make install
```

Use the same commands for [libctl](https://libctl.readthedocs.io), then [MPB](https://mpb.readthedocs.io), then Meep.

You are done, and can now run Meep just by typing `meep`. You can run `make check` in the meep directory if you want to perform a self-test.

To build the latest version of Meep from source on MacOS Sierra, follow these [instructions](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05811.html).

Unix Installation Basics
------------------------

### Installation Paths

First, let's review some important information about installing software on Unix systems, especially in regards to installing software in non-standard locations. None of these issues are specific to Meep, but they've caused a lot of confusion among users.

Most of the software below, including Meep, installs under `/usr/local` by default. That is, libraries go in `/usr/local/lib`, programs in `/usr/local/bin`, etc. If you don't have `root` privileges on your machine, you may need to install somewhere else, e.g. under `$HOME/install` (the `install/` subdirectory of your home directory). Most of the programs below use a GNU-style `configure` script, which means that all you would do to install there would be:

```sh
 ./configure --prefix=$HOME/install
```

when configuring the program. The directories `$HOME/install/lib` etc. are created automatically as needed.

#### Paths for Configuring

There are two further complications. First, if you install in a non-standard location and `/usr/local` is considered non-standard by some proprietary compilers, you will need to tell the compilers where to find the libraries and header files that you installed. You do this by setting two environment variables:

```bash
 export LDFLAGS="-L/usr/local/lib"
 export CPPFLAGS="-I/usr/local/include"
```

Of course, substitute whatever installation directory you used. Do this **before** you run the `configure` scripts, etcetera. You may need to include multiple `-L` and `-I` flags separated by spaces if your machine has stuff installed in several non-standard locations. Bourne shell users (e.g. `bash` or `ksh`) should use the `export FOO=bar` syntax instead of `csh`'s `setenv FOO bar`, of course.

You might also need to update your `PATH` so that you can run the executables you installed although `/usr/local/bin/` is in the default `PATH` on many systems. e.g. if we installed in our home directory as described above, we would do:

```bash
 export PATH="$HOME/install/bin:$PATH"
```

#### Paths for Running (Shared Libraries)

Second, many of the packages installed below (e.g. Guile) are installed as shared libraries. You need to make sure that your runtime linker knows where to find these shared libraries. The bad news is that every operating system does this in a slightly different way. The good news is that, when you run `make install` for the packages involving shared libraries, the output includes the necessary instructions specific to your system, so pay close attention! It will say something like `add LIBDIR to the <foobar> environment variable`, where `LIBDIR` will be your library installation directory (e.g. `/usr/local/lib`) and `<foobar>` is some environment variable specific to your system (e.g. `LD_LIBRARY_PATH` on some systems, including Linux). For example, you might do:

```bash
 export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
```

Note that we just add to the library path variable, and don't replace it in case it contains stuff already. If you use Linux and have `root` privileges, you can instead simply run `/sbin/ldconfig`, first making sure that a line `/usr/local/lib` (or whatever) is in `/etc/ld.so.conf`.

If you don't want to type these commands every time you log in, you can put them in your `~/.cshrc` file or `~/.profile`, or `~/.bash_profile`, depending on your shell.

### Fun with Fortran

Meep, along with many of the libraries it calls, is written in C or C++, but it also calls libraries such as BLAS and LAPACK (see below) that are usually compiled from Fortran. This can cause some added difficulty because of the various linking schemes used by Fortran compilers. Our `configure` script attempts to detect the Fortran linking scheme automatically, but in order for this to work you must use the same Fortran compiler and options with Meep as were used to compile BLAS/LAPACK.

By default, Meep looks for a vendor Fortran compiler first (`f77`, `xlf`, etcetera) and then looks for GNU `g77`. In order to manually specify a Fortran compiler `foobar` you would configure it with `./configure F77=foobar ...`.

If, when you compiled BLAS/LAPACK, you used compiler options that alter the linking scheme (e.g. `g77`'s `-fcase-upper` or `-fno-underscoring`), you will need to pass the same flags to Meep via `./configure FFLAGS=...flags... ...`.

### Picking a Compiler

It is often important to be consistent about which compiler you employ. This is especially true for C++ software. To specify a particular C compiler `foo`, configure with `./configure CC=foo`; to specify a particular C++ compiler `foo++`, configure with `./configure CXX=foo++`; to specify a particular Fortran compiler `foo90`, configure with `./configure F77=foo90`.

### Linux and BSD Binary Packages

If you are installing on your personal Linux or BSD machine, then precompiled binary packages are likely to be available for many of these packages, and may even have been included with your system. On Debian systems, the packages are in `.deb` format and the built-in `apt-get` program can fetch them from a central repository. On Red Hat, SuSE, and most other Linux-based systems, binary packages are in RPM format.  OpenBSD has its "ports" system, and so on.

**Do not compile something from source if an official binary package is available.**  For one thing, you're just creating pain for yourself.  Worse, the binary package may already be installed, in which case installing a different version from source will just cause trouble.

One thing to watch out for is that libraries like LAPACK, Guile, HDF5, etcetera, will often come split up into two or more packages: e.g. a `guile` package and a `guile-devel` package. You need to install **both** of these to compile software using the library.

BLAS and LAPACK (recommended)
-----------------------------

BLAS and LAPACK libraries are required in order to install [Harminv](https://github.com/stevengj/harminv/blob/master/README.md). Harminv is not *required* for Meep, but is strongly recommended for use in resonant-mode computation.

Note also that Meep's usage of BLAS/LAPACK, via Harminv, is not generally performance-critical. So, it doesn't matter too much whether you install an especially optimized BLAS library. However, it makes a big difference if you also use [MPB](https://mpb.readthedocs.io).

### BLAS

The first thing you must have on your system is a BLAS implementation. "BLAS" stands for "Basic Linear Algebra Subroutines," and is a standard interface for operations like matrix multiplication. It is designed as a building-block for other linear-algebra applications, and is used both directly by LAPACK (see below). By using it, we can take advantage of many highly-optimized implementations of these operations that have been written to the BLAS interface. Note that you will need implementations of BLAS levels 1-3.

You can find more BLAS information, as well as a basic implementation, on its [homepage](http://www.netlib.org/blas/). Once you get things working with the basic BLAS implementation, it might be a good idea to try and find a more optimized BLAS code for your hardware. Vendor-optimized BLAS implementations are available as part of the Intel MKL, HP CXML, IBM ESSL, SGI sgimath, and other libraries. An excellent, high-performance, free-software BLAS implementation is  [OpenBLAS](http://www.openblas.net). Another is [ATLAS](http://math-atlas.sourceforge.net/).

Note that the generic BLAS does not come with a `Makefile`; compile it with something like: </nowiki>

```sh
  wget http://www.netlib.org/blas/blas.tgz
  gunzip blas.tgz
  tar xf blas.tar
  cd BLAS
  f77 -c -O3 *.f   # compile all of the .f files to produce .o files
  ar rv libblas.a *.o    #  combine the .o files into a library
  su -c "cp libblas.a /usr/local/lib"   # switch to root and install
```

Replace `-O3` with your favorite optimization options. On Linux, this could be `g77 -O3 -fomit-frame-pointer -funroll-loops -malign-double`. Note that MPB looks for the standard BLAS library with `-lblas`, so the library file should be called `libblas.a` and reside in a standard directory like `/usr/local/lib`. See also below for the `--with-blas=lib` option to MPB's `configure` script, to manually specify a library location.

### LAPACK

LAPACK, the Linear Algebra PACKage, is a standard collection of routines, built on BLAS, for more-complicated (dense) linear algebra operations like matrix inversion and diagonalization. You can download LAPACK from its [homepage](http://www.netlib.org/lapack).

Note that Meep looks for LAPACK by linking with `-llapack`. This means that the library must be called `liblapack.a` and be installed in a standard directory like `/usr/local/lib`. Alternatively, you can specify another directory via the `LDFLAGS` environment variable as described earlier. See also below for the `--with-lapack=''lib''` option to our `configure` script, to manually specify a library location.

We currently recommend installing OpenBLAS which includes LAPACK so you do not need to install it separately.

Harminv (recommended)
---------------------

To use Meep to extract resonant frequencies and decay rates, you must install [Harminv](https://github.com/stevengj/harminv/blob/master/README.md) which requires BLAS and LAPACK.

See the [Harminv installation](https://github.com/stevengj/harminv/blob/master/doc/installation.md) instructions.

Guile (recommended)
-----------------------

Guile is required in order to use the Scheme interface, and is strongly recommended. If you don't install it, you can only use the C++ interface.

Guile is an extension/scripting language implementation based on Scheme, and we use it to provide a rich, fully-programmable user interface with minimal effort. It's free, of course, and you can download it from the [Guile homepage](http://www.gnu.org/software/guile/). Guile is typically included with Linux systems.

- **Important:** Most Linux distributions come with Guile already installed. You can check by seeing whether you can run `guile --version` from the command line. In that case, do **not** install your own version of Guile from source &mdash; having two versions of Guile on the same system will cause problems. However, by default most distributions install only the Guile libraries and not the programming headers &mdash; to compile libctl and MPB, you should install the **guile-devel** or **guile-dev** package. Note: SWIG does [not support the latest version of Guile](https://github.com/swig/swig/issues/312) which requires that an older version (2.0.11) be used.  Meep contains a workaround for newer versions as well.

libctl (recommended)
--------------------

[libctl](https://libctl.readthedocs.io), which requires Guile, is required to use the Scheme interface, and is strongly recommended. If you don't install it, you can only use the C++ interface. libctl version **3.2.2 or later** is required.

Instead of using Guile directly, we separated much of the user interface code into a package called libctl, in the hope that this might be more generally useful. libctl automatically handles the communication between the program and Guile, converting complicated data structures and so on, to make it even easier to use Guile to control scientific applications. Download libctl from the [libctl page](https://libctl.readthedocs.io), unpack it, and run the usual `configure`, `make`, `make install` sequence. You'll also want to browse the [libctl manual](https://libctl.readthedocs.io), as this will give you a general overview of what the user interface will be like.

If you are not the system administrator of your machine, and/or want to install libctl somewhere else like your home directory, you can do so with the standard `--prefix=dir` option to `configure`. The default prefix is `/usr/local`. In this case, however, you'll need to specify the location of the libctl shared files for the MPB or Meep package, using the `--with-libctl=dir/share/libctl` option to our `configure` script.

MPI (parallel machines)
-----------------------

Optionally, Meep is able to run on a distributed-memory parallel machine, and to do this we use the standard message-passing interface (MPI). Most commercial supercomputers already have an MPI implementation installed. The recommended implementation is [Open MPI](http://www.open-mpi.org/). MPI is **not required** to compile the serial version of Meep.

In order for the MPI version of the Scheme interface to run successfully, we have a slightly nonstandard requirement: each process must be able to read from the disk. This way, Guile can boot for each process and they can all read your control file in parallel. Most commercial supercomputers satisfy this requirement. On the other hand, the C++ interface to Meep does not have this requirement.

If you use Meep with MPI, you should compile HDF5 with MPI support as well (see below).

As described below, when you configure Meep with MPI support (`--with-mpi`), it installs itself as **meep-mpi**.

HDF5 (recommended)
------------------

Meep outputs its fields and other volumetric data in the HDF5 format, so you must install the HDF5 libraries if you want to visualize the fields.

HDF is a widely-used, free, portable library and file format for multi-dimensional scientific data, developed in the National Center for Supercomputing Applications (NCSA) at the University of Illinois. You can get HDF and learn about it on the [HDF homepage](https://www.hdfgroup.org).

There are two incompatible versions of HDF, HDF4 and HDF5 (no, not HDF1 and HDF2). We require the newer version, HDF5, which is supported by a number scientific of visualization tools, including [h5utils](https://github.com/stevengj/h5utils/blob/master/README.md) utilities.

HDF5 includes parallel I/O support under MPI, which can be enabled by configuring it with `--enable-parallel`. You may also have to set the `CC` environment variable to `mpicc`. Unfortunately, the parallel HDF5 library then does not work with serial code, so you have may have to choose one or the other.

We have some hacks in Meep and MPB so that they can do parallel I/O even with the serial HDF5 library. These hacks work okay when you are using a small number of processors, but on large supercomputers we strongly recommend using the parallel HDF5.

**Note:** If you have a version of HDF5 compiled with MPI parallel I/O support, then you need to use the MPI compilers to link to it, even when you are compiling the serial versions of Meep or MPB.  Just use `./configure CC=mpicc CXX=mpic++` or whatever your MPI compilers are when configuring.

Meep
----

Once you've installed all of the prerequisites, you can install Meep via:

```sh
./configure
make
sudo make install
```

Assuming you've set your `LDFLAGS` etcetera, the configure script should find all of the libraries you've installed and, with luck, compile successfully. The `sudo` in the last command uses administrator privileges to install the binaries in standard system directories. Alternatively, you can just use `make install` if you have used `--prefix` to change the installation directory to something like your home directory. This is described below. To make sure Meep is working, you can run its test suite via:

```sh
make check
```

The configure script accepts several flags to modify its behavior.

**`--prefix=dir`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Install into `dir/bin`, etcetera, as described above.


**`--with-mpi`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Attempt to compile a [parallel version of Meep](Parallel_Meep.md) using MPI; the resulting program will be installed as `meep-mpi.` Requires MPI to be installed, as described above. Does **not** compile the serial Meep. If you want that, you will have to `make distclean` and install the serial version separately. Note that the configure script attempts to automatically detect how to compile MPI programs, but this may fail if you have an unusual version of MPI or if you have several versions of MPI installed and you want to select a particular one. You can control the version of MPI selected by setting the `MPICXX` variable to the name of the compiler to use and the `MPILIBS` variable to any additional libraries that must be linked (e.g., `./configure MPICXX=foompiCC MPILIBS=-lfoo ...`).

**`--with-libctl=dir`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If libctl was installed in a nonstandard location (i.e. neither `/usr` nor `/usr/local`), you need to specify the location of the libctl directory, *`dir`*. This is either `prefix/share/libctl`, where `prefix` is the installation prefix of libctl, or the original libctl source code directory. To configure *without* the libctl/Guile interface, use `--without-libctl`.

**`--with-blas=lib`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The `configure` script automatically attempts to detect accelerated BLAS libraries, like DXML (DEC/Alpha), SCSL and SGIMATH (SGI/MIPS), ESSL (IBM/PowerPC), ATLAS, and PHiPACK. You can, however, force a specific library name to try via `--with-blas=lib`.

**`--with-lapack=lib`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Cause the `configure` script to look for a LAPACK library called *`lib`*. The default is to use `-llapack`.

**`--enable-debug`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Compile for debugging, adding extra runtime checks and so on.

**`--enable-shared`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Install the Meep libraries as shared libraries (i.e. dynamically linked) rather than as static libraries. This is off by default because shared libraries require the user to configure their runtime linker paths correctly (see "Paths for Running" above).

**`--without-hdf5`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Install Meep without support for the HDF5 libraries (this means you won't be able to output fields and so on).

**`--enable-portable-binary`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
By default, Meep's `configure` script picks compiler flags to optimize Meep as much as possible for the machine you are compiling on. If you wish to run the *same compiled executable* on other machines, however, you need to tell it not to pick compiler flags that use features specific to your current processor. In this case you should pass `--enable-portable-binary` to `configure`. (This option is mainly useful for people creating binary packages for Debian, Fedora, etcetera.)

**`--with-gcc-arch=arch`, `--without-gcc-arch`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
By default, Meep's configure script tries to guess the gcc `-march` flag for the system you are compiling on using `-mtune` instead when `--enable-portable-binary` is specified. If it guesses wrong, or if you want to specify a different architecture, you can pass it here. If you want to omit `-march`/`-mtune` flags entirely, pass `--without-gcc-arch`.

Python Meep
-----------

An alpha version of Python Meep is available, currently only for the serial (i.e., non-MPI) version. This can be installed on Ubuntu in two ways: [building from source](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05850.html) or as a [pre-compiled package using Conda](https://gist.github.com/ChristopherHogan/c22bb7cc7248f7595c4b09d0e3b16735).

Meep for Developers
-------------------

If you want to modify the Meep source code, you will want to have a number of additional packages, most importantly:

-   The [Git](https://git-scm.com/) version-control system.

Once you have Git, you can grab the latest development version of Meep with:

```sh
 git clone https://github.com/stevengj/meep.git
```

This gives you a fresh, up-to-date Meep repository in a directory `meep`. See [git-scm.com](https://git-scm.com/) for more information on using Git; perhaps the most useful command is `git pull`, which you can execute periodically to get any new updates to the development version.

Git will give you an absolutely minimal set of sources; to create a usable Meep directory, you should run:

```sh
sh autogen.sh
make
```

in the `meep` directory. And subsequently, if you are editing the sources you should include `--enable-maintainer-mode` whenever you reconfigure. To do this, however, you will need a number of additional packages beyond those listed above:

-   GNU [autoconf](https://www.gnu.org/software/autoconf/autoconf.html), [automake](https://www.gnu.org/software/automake/), and [libtool](https://www.gnu.org/software/libtool/libtool.html) &mdash; these are used to create the Makefiles and configure scripts, and to build shared libraries.
-   [SWIG](http://www.swig.org/) &mdash; the Scheme/libctl interface to Meep is largely generated by a program called *SWIG* (Simple Wrapper and Interface Generator).

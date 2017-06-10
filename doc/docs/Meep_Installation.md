---
# Meep Installation
---

Much like [MPB](http://ab-initio.mit.edu/wiki/index.php/MPB_Installation), the main effort in installing Meep lies in installing the various prerequisite packages. This requires some understanding of how to install software on Unix systems.

It is also possible to install Meep on [Windows](https://en.wikipedia.org/wiki/Microsoft_Windows) systems by using a free Unix-compatibility environment such as [Cygwin](http://www.cygwin.org/). For more information, see these [step-by-step instructions](http://novelresearch.weebly.com/installing-meep-in-windows-8-via-cygwin.html).

For those installing Meep on a parallel supercomputer, a note of caution: most supercomputers have more than one different compiler installed, and different versions of libraries compiled with different compilers. Meep is written in C++, and it is almost impossible to mix C++ code compiled by different compilers—pick one set of compilers by one vendor and stick with it consistently!

**Note**: The latest, pre-installed versions of Meep and MPB running on Ubuntu can also be accessed on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetuscloud.com/launchsims.html).

Installation on GNU/Linux
-------------------------

For most [GNU/Linux distributions](https://en.wikipedia.org/wiki/Linux_distribution), there should be precompiled packages for most of Meep's prerequisites below, and we highly recommend installing those prerequisites using the available packages for your system whenever possible. Using precompiled packages means that you don't have to worry about how to install things manually, you are using packages which have already been tweaked to work well with your system, and (usually) your packages will be automatically upgraded when you upgrade the rest of your system.

BLAS and LAPACK (possibly as part of a package for [Atlas BLAS](https://en.wikipedia.org/wiki/Automatically_Tuned_Linear_Algebra_Software)), Guile, MPI, and HDF5 are all almost certain to be available as precompiled packages. One thing to be careful of is that many distributions split packages into two parts: one main package for the libraries and programs, and a "devel" package for [header files](https://en.wikipedia.org/wiki/Header_file) and other things needed to compile software using those libraries. You will need to install *both*. So, for example, you will probably need both a `guile` package (probably installed by default) and a `guile-dev` or `guile-devel` package (probably *not* installed by default), and similarly for HDF5 etcetera. You will probably also want to install a `libpng-dev` or `libpng-devel` package in order to compile the `h5topng` utility in [H5utils](http://ab-initio.mit.edu/wiki/index.php/H5utils).

The easiest installation is on [Debian GNU/Linux](https://en.wikipedia.org/wiki/Debian_GNU/Linux) (and [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu_(operating_system))), which has precompiled packages for Meep. On Debian, just do:

```
apt-get install meep h5utils
```


to install everything. (See also [Meep Download](Meep_Download.md).)

Installation on MacOS X
-----------------------

Since [MacOS X](https://en.wikipedia.org/wiki/Mac_OS_X) is, at its heart, a Unix system, one can, in principle compile and install Meep and all its prerequisites just as on any other Unix system. However, life will be much easier if you use the [Homebrew](https://en.wikipedia.org/wiki/Homebrew_(package_management_software)) package to install most of the prerequisites, since it will handle dependencies and other details for you. (You will need [administrator privileges](http://support.apple.com/kb/PH3920) on your Mac.)

The first steps are:

-   Install [Xcode](https://en.wikipedia.org/wiki/Xcode), the development/compiler package from Apple, free from the [Apple Xcode web page](https://developer.apple.com/xcode/).
-   Install Homebrew: download from the [Homebrew web site](http://brew.sh/) and follow the instructions there.
-   Run the following commands in the terminal to compile and install Meep prerequisites (this may take a while to complete because it will install lots of other stuff first)

```
brew doctor
brew install homebrew/science/hdf5 homebrew/science/openblas guile fftw h5utils
```


Now, install the Harminv, Libctl, MPB, and Meep packages from source. Download [Harminv](http://ab-initio.mit.edu/wiki/index.php/harminv) and, in the `harminv` directory, do:

```
./configure && make && make install
```


Use the same commands for [Libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl), then [MPB](http://ab-initio.mit.edu/wiki/index.php/MPB), then Meep.

**Note:** Meep currently does not compile on OSX 10.9, due to conflicts with the latest C++ standard. We hope to fix this with a new release shortly.

You are done, and can now run Meep just by typing `meep`. You can run `make` `check` in the meep directory if you want to perform a self-test.

BLAS and LAPACK (recommended)
-----------------------------

BLAS and LAPACK libraries are required in order to install [Harminv](http://ab-initio.mit.edu/wiki/index.php/harminv). Harminv is not *required* for Meep, but is strongly recommended for use in resonant-mode computation.

Note also that Meep's usage of BLAS/LAPACK, via Harminv, is not generally performance-critical. So, it doesn't matter too much whether you install an especially optimized BLAS library. (It makes a big difference if you also use [MPB](http://ab-initio.mit.edu/wiki/index.php/MPB), though.)

Harminv (recommended)
---------------------

To use Meep to extract resonant frequencies and decay rates, you must install [Harminv](http://ab-initio.mit.edu/wiki/index.php/harminv) (which requires BLAS and LAPACK).

See the [Harminv installation](http://ab-initio.mit.edu/wiki/index.php/Harminv_installation) instructions.

GNU Guile (recommended)
-----------------------

Guile is required in order to use the ordinary libctl-based front-end to Meep, and is strongly recommended. (If you don't install it, you can only use the C++ interface to Meep.)

-   **Important:** Meep currently requires Guile version 1.6 or later. You can see which version of Guile you have by running `guile` `--version`.

libctl (recommended)
--------------------

[libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl), which requires Guile, is required to use the ordinary libctl-based front-end to Meep, and is strongly recommended. (If you don't install it, you can only use the C++ interface to Meep.) Meep requires version **3.2 or later** of libctl.

MPI (parallel machines)
-----------------------

Optionally, Meep is able to run on a distributed-memory parallel machine, and to do this we use the standard MPI message-passing interface. You can learn about MPI from the [MPI Home Page](http://www-unix.mcs.anl.gov/mpi/). Most commercial supercomputers already have an MPI implementation installed; two free MPI implementations that you can install yourself are [MPICH](http://www-unix.mcs.anl.gov/mpi/mpich/) and [LAM](http://www.lam-mpi.org/); a promising new implementation is [Open MPI](http://www.open-mpi.org/). MPI is *not required* to compile the ordinary, uniprocessor version of our software.

In order for the MPI version of our [libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl)/Scheme front-end to run successfully, we have a slightly nonstandard requirement: each process must be able to read from the disk. (This way, Guile can boot for each process and they can all read your control file in parallel.) Many (most?) commercial supercomputers, Linux [Beowulf](http://www.beowulf.org) clusters, etcetera, satisfy this requirement. On the other hand, the C++ interface to Meep does not have this requirement.

If you use Meep with MPI, you should compile HDF5 with MPI support as well (see below).

As described below, when you configure Meep with MPI support (`--with-mpi`), it installs itself as `meep-mpi`.

HDF5 (recommended)
------------------

Meep outputs its fields and other volumetric data in the HDF5 format, so you must install the HDF5 libraries if you want to visualize the fields.

Meep
----

Once you've installed all of the prerequisites, you can install Meep via:

```
./configure
make
su -c "make install"
```


Assuming you've set your `LDFLAGS` etcetera, the configure script should find all of the libraries you've installed and (with luck) compile successfully. The `su` in the last command switches to root for installation; you can just use `make` `install` if you have used `--prefix` (below) to change the installation directory to something like your home directory. To make sure Meep is working, you can run its test suite via:

```
make check
```


The configure script accepts several flags to modify its behavior:

`--prefix=`*`dir`*
Install into *`dir`*`/bin`, etcetera, as described above.

```
--with-mpi
```

Attempt to compile a [parallel version of Meep](Parallel_Meep.md) using MPI; the resulting program will be installed as `meep-mpi`. Requires MPI to be installed, as described above. Does *not* compile the serial Meep; if you want that, you will have to `make` `distclean` and install the uniprocessor Meep separately. Note that the configure script attempts to automatically detect how to compile MPI programs, but this may fail if you have an unusual version of MPI or if you have several versions of MPI installed and you want to select a particular one. You can control the version of MPI selected by setting the `MPICXX` variable to the name of the compiler to use and the `MPILIBS` variable to any additional libraries that must be linked (e.g. `./configure` `MPICXX=foompiCC` `MPILIBS="-lfoo"` ...).

`--with-libctl=`*`dir`*
If libctl was installed in a nonstandard location (i.e. neither `/usr` nor `/usr/local`), you need to specify the location of the libctl directory, *`dir`*. This is either *`prefix`*`/share/libctl`, where *`prefix`* is the installation prefix of libctl, or the original libctl source code directory. To configure *without* the libctl/Guile interface, use `--without-libctl`.

`--with-blas=`*`lib`*
The `configure` script automatically attempts to detect accelerated BLAS libraries, like DXML (DEC/Alpha), SCSL and SGIMATH (SGI/MIPS), ESSL (IBM/PowerPC), ATLAS, and PHiPACK. You can, however, force a specific library name to try via `--with-blas=`*`lib`*.

`--with-lapack=`*`lib`*
Cause the `configure` script to look for a LAPACK library called *`lib`* (the default is to use `-llapack`).

```
--enable-debug
```

Compile for debugging, adding extra runtime checks and so on.

```
--enable-shared
```

Install the meep libraries as shared libraries (i.e. dynamically linked) rather than as static libraries. This is off by default because shared libraries require the user to configure their runtime linker paths correctly (see "Paths for Running" above).

```
--without-hdf5
```

Install Meep without support for the HDF5 libraries (this means you won't be able to output fields and so on).

```
--enable-portable-binary
```

By default, Meep's `configure` script picks compiler flags to optimize Meep as much as possible for the machine you are compiling on. If you wish to run the *same compiled executable* on other machines, however, you need to tell it not to pick compiler flags that use features specific to your current processor. In this case you should pass `--enable-portable-binary` to `configure`. (This option is mainly useful for people creating binary packages for Debian, Fedora, etcetera.)

`--with-gcc-arch=`*`arch`*, --without-gcc-arch
By default, Meep's configure script tries to guess the gcc `-march` flag for the system you are compiling on (using `-mtune` instead when `--enable-portable-binary` is specified). If it guesses wrong, or if you want to specify a different architecture, you can pass it here. If you want to omit `-march`/`-mtune` flags entirely, pass `--without-gcc-arch`.

Meep for developers
-------------------

If you want to modify the Meep source code, you will want to have a number of additional packages, most importantly:

-   The [Git](http://git-scm.com/) version-control system.

Once you have Git, you can grab the latest development version of Meep with:

`git clone `[`git://github.com/stevengj/meep`](git://github.com/stevengj/meep)

This gives you a fresh, up-to-date Meep repository in a directory `meep`. See [git-scm.com](http://git-scm.com/) for more information on using Git; perhaps the most useful command is `git` `pull`, which you can execute periodically to get any new updates to the development version.

Git will give you an absolutely minimal set of sources; to create a usable Meep directory, you should run:

```
sh autogen.sh
make
```


in the `meep` directory. (And subsequently, if you are editing the sources you should include `--enable-maintainer-mode` whenever you reconfigure.) To do this, however, you will need a number of additional packages beyond those listed above:

-   GNU [autoconf](http://www.gnu.org/software/autoconf/), [automake](http://sources.redhat.com/automake/), and [libtool](http://www.gnu.org/software/libtool/libtool.html) — these are used to create the Makefiles and configure scripts, and to build shared libraries.
-   [SWIG](http://www.swig.org/) — the Scheme/libctl interface to Meep is largely generated by a program called *SWIG* (Simple Wrapper and Interface Generator). We currently require SWIG version 1.3.25 or later. Moreover, if you are using 1.3.27 or earlier, you must patch the file `Source/Modules/guile.cxx` with [this bug fix](http://cvs.sourceforge.net/viewcvs.py/swig/SWIG/Source/Modules/guile.cxx?r1=1.33&r2=1.34).

[Category:Meep](Meep.md)

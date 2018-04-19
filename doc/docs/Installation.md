---
# Installation
---

**Note**: Installing Meep from source can be challenging for novice users. As a simple workaround, the latest version of Meep preinstalled on Ubuntu can be accessed on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetus.com/launchsims.html). For easy access to the Python interface, we provide a binary installation in the form of Conda packages. Details can be found [below](#conda-packages).

[TOC]

Conda Packages
---------------

The recommended way to install PyMeep is using the [Conda](https://conda.io/docs/) package manager. Binary packages for serial and parallel PyMeep on Linux and macOS are currently available (64 bit architectures only), and are updated with each MEEP release. The easiest way to get started is to install [Miniconda](https://conda.io/miniconda.html), which comes with everything necessary to create Python environments with Conda. For example, to install Miniconda with Python 3 on Linux:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p <desired_prefix>
export PATH=<desired_prefix>/bin:$PATH
```

Next, we create a Conda environment for PyMeep to isolate it from other Python libraries that may be installed.

```bash
conda create -n mp -c chogan -c defaults -c conda-forge pymeep
```

This creates an environment called "mp" (you can name this anything you like) with PyMeep and all its dependencies. This will default to the version of Python in your Miniconda installation (Python 3 for us since we installed Miniconda3), but if you want to work with Python 2, just add `python=2` to the end of the command. We hope to move everything to conda-forge to simplify the channel selection, but currently we need to pull dependencies from three different channels (the -c arguments), and the order they are specified in is very important.

Next, we need to activate the environment before we can start using it.

```bash
source activate mp
```

Now, `python -c 'import meep'` should work, and you can try running some of the examples in the `meep/python/examples` directory.

Installing parallel PyMeep follows the same pattern, but the package is called `pymeep-parallel`.

```bash
conda create -n pmp -c chogan -c defaults -c conda-forge pymeep-parallel
source activate pmp
```

The environment includes `mpi4py`, so you can run an MPI job with 4 processes like this:

```bash
mpirun -np 4 python <script_name>.py
```

If you run into issues, make sure your `PYTHONPATH` environment variable is unset.

*Note:* For pymeep-parallel on macOS, a [bug](https://github.com/open-mpi/ompi/issues/2956) in openmpi requires that the environment variable `TMPDIR` be set to a short path like `/tmp`. Without this workaround, you may see errors similar to this:

```bash
[laptop:68818] [[53415,0],0] ORTE_ERROR_LOG: Bad
parameter in file ../../orte/orted/pmix/pmix_server.c at line 264

[laptop:68818] [[53415,0],0] ORTE_ERROR_LOG: Bad
parameter in file ../../../../../orte/mca/ess/hnp/ess_hnp_module.c at line
666
```

Installation on Linux
-------------------------

For most [Linux distributions](https://en.wikipedia.org/wiki/Linux_distribution), there should be precompiled packages for most of Meep's prerequisites below, and we *highly* recommend installing those prerequisites using the available packages for your system whenever possible. Using precompiled packages means that you don't have to worry about how to install things manually. You are using packages which have already been tweaked to work well with your system, and usually your packages will be automatically upgraded when you upgrade the rest of your system.

The following precompiled packages are available: BLAS and LAPACK possibly as part of a package for [Atlas BLAS](https://en.wikipedia.org/wiki/Automatically_Tuned_Linear_Algebra_Software), Guile, MPI, and HDF5. One thing to be careful of is that many distributions split packages into two parts: one main package for the libraries and programs, and a **devel** package for [header files](https://en.wikipedia.org/wiki/Header_file) and other things needed to compile software using those libraries. You will need to install **both**. So, for example, you will probably need both a `guile` package (probably installed by default) and a `guile-dev` or `guile-devel` package (probably *not* installed by default), and similarly for HDF5 etcetera. You will probably also want to install a `libpng-dev` or `libpng-devel` package in order to compile the `h5topng` utility in [h5utils](https://github.com/stevengj/h5utils/blob/master/README.md).

The easiest installation is on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu_(operating_system)) which has precompiled packages for Meep:

```sh
apt-get install meep h5utils
```

Gzipped tarballs of stable versions of the source are available on the [releases page](https://github.com/stevengj/meep/releases). See [Build From Source](Build_From_Source) for more information.

Installation on macOS 
-----------------------

Since [macOS](https://en.wikipedia.org/wiki/macOS) is, at its heart, a Unix system, one can, in principle compile and install Meep and all its prerequisites just as on any other Unix system. However, this process is much easier using the [Homebrew](https://en.wikipedia.org/wiki/Homebrew_(package_management_software)) package to install most of the prerequisites, since it will handle dependencies and other details for you. You will need [administrator privileges](http://support.apple.com/kb/PH3920) on your Mac.

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

Use the same commands for [libctl](https://libctl.readthedocs.io), then [MPB](https://mpb.readthedocs.io), then Meep. See [Build From Source](Buil_From_Source.md) for more detailed information.

You are done, and can now run Meep just by typing `meep`. You can run `make check` in the meep directory if you want to perform a self-test.

To build the latest version of Meep from source on macOS Sierra, follow these [instructions](https://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg05811.html).


---
# Download
---

[TOC]

GitHub Source Repository
------------------------

The [source repository](https://github.com/NanoComp/meep) is hosted on GitHub along with gzipped tarballs of [official (stable) releases](https://github.com/NanoComp/meep/releases).

Refer to [NEWS](https://github.com/NanoComp/meep/blob/master/NEWS.md) for a list of the latest changes, and be sure to read [Installation](Installation.md) for how to compile and install it.

To receive notifications when new versions are released, subscribe to the [meep-announce](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce) mailing list.

Precompiled Packages for Ubuntu
-------------------------------

Precompiled packages of [Meep version 1.17.1](https://github.com/NanoComp/meep/releases/tag/v1.17.1) (January 2021) with Python interface will be available for [Ubuntu 21.04 ("Hirsute Hippo")](https://packages.ubuntu.com/hirsute/python3-meep) in April 2021. We recommend Ubuntu as Meep and all of its dependencies will be able to be installed using just one line:

```sh
sudo apt-get install python3-meep h5utils
```

You will also be able to install the [parallel version of Meep](https://packages.ubuntu.com/hirsute/python3-meep-openmpi) which is based on [OpenMPI](https://www.open-mpi.org/) using:

```sh
sudo apt-get install python3-meep-openmpi
```

These upcoming Meep packages for Ubuntu 21.04 are derived from the Debian 11 ("Bullseye") packages ([serial](https://packages.debian.org/bullseye/python3-meep) and [parallel](https://packages.debian.org/bullseye/python3-meep-openmpi)). Debian 11 is currently the testing distribution.
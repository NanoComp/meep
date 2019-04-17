---
# Download
---

[TOC]

GitHub Source Repository
------------------------

The source repository is available via [GitHub](https://github.com/NanoComp/meep).

The current stable release can be obtained from:

-   <https://github.com/NanoComp/meep/releases>

Refer to [NEWS](https://github.com/NanoComp/meep/blob/master/NEWS.md) for a list of the latest changes, and be sure to read the [Installation](Installation.md) section for how to compile and install it.

To receive notifications when new versions are released, subscribe to the **meep-announce** mailing list:

-   [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce)

Precompiled Packages for Ubuntu
-------------------------------

Precompiled packages of Meep version 1.3 (September 2017) are available for [Ubuntu](https://packages.ubuntu.com/search?keywords=meep). We recommend Ubuntu as Meep and all of its dependencies can be installed using just one line:

```sh
sudo apt-get install meep h5utils
```

You can also install the parallel version of Meep which is based on [OpenMPI](https://www.open-mpi.org/) using:

```sh
sudo apt-get install meep-openmpi
```

The Meep package for Ubuntu is in the process of being updated and will likely appear in Ubuntu 19.10 as derived from the [Debian package](https://packages.debian.org/search?searchon=names&keywords=meep).

Amazon Web Services (AWS)
-------------------------

The latest stable version of Meep preinstalled on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) 16.04 can be accessed for free on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as an [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS) provided by [Simpetus](http://www.simpetus.com/launchsims.html).

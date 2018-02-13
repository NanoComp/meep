---
# Download
---

The latest development source is available on [GitHub](https://github.com/stevengj/meep).

The current stable release is **version 1.4.3** which can be downloaded from:

-   <https://github.com/stevengj/meep/releases>

Refer to the [NEWS file](https://github.com/stevengj/meep/blob/master/NEWS.md) to see what's new in this version, and be sure to read the [Installation](Installation.md) section for how to compile and install it.

Please subscribe to the **meep-announce** mailing list to receive notifications when new versions are released:

-   [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce)

Meep on Amazon Web Services (AWS)
---------------------------------

The most recent, stable version of Meep preinstalled on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) can be accessed for free on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as an [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS) provided by [Simpetus](http://www.simpetus.com/launchsims.html).

Precompiled Meep packages for Ubuntu
------------------------------------

Precompiled packages of Meep are available for [Ubuntu](https://packages.ubuntu.com/search?keywords=meep). We recommend Ubuntu as Meep and all of its dependencies can be installed using just one line:

```sh
sudo apt-get install meep h5utils
```

You can also install the [parallel version of Meep](http://packages.debian.org/testing/science/meep-mpi-default) which is based on [Open MPI](https://www.open-mpi.org/) using:

```sh
sudo apt-get install meep-mpi-default
```
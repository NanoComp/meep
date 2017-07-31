---
# Download
---

The latest development sources are available on [GitHub](https://github.com/stevengj/meep).

The latest release of Meep is **version 1.3** which can be downloaded from:

-   <https://github.com/stevengj/meep/releases>

Meep is free software under the [GNU GPL](License_and_Copyright.md).

Refer to the [Release Notes](Release_Notes.md) to see what's new in this version, and be sure to read the [Installation](Installation.md) manual for how to compile and install it.

Please subscribe to the **meep-announce** mailing list to receive notifications when new versions are released:

-   [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce)

Meep on Amazon Web Services (AWS)
---------------------------------

The latest version of Meep preinstalled on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) can be accessed for free on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as an [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS) provided by [Simpetus](http://www.simpetuscloud.com/launchsims.html).

Precompiled Meep packages for Ubuntu
------------------------------------

Precompiled packages of Meep are available for [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) as **meep**. The Ubuntu 16.10 package is available in the [science](https://packages.ubuntu.com/yakkety/meep) repository. We highly recommend using Ubuntu as the Meep software and all of its dependencies can be installed using just one line:

```
sudo apt-get install meep h5utils
```

You can also install the [parallel version of Meep](http://packages.debian.org/testing/science/meep-mpi-default) which is based on [Open-MPI](https://www.open-mpi.org/) using:

```
sudo apt-get install meep-mpi-default
```

Python User Interface
----------------

A native Python user interface for Meep is currently under development by [Simpetus](http://www.simpetuscloud.com).

A non-native version was developed by the [Photonics Research Group](http://photonics.intec.ugent.be/) at Ghent University, in particular thanks to Emmanuel Lambert, Martin Fiers, Shavkat Nizamov, Martijn Tassaert, Peter Bienstman, and Wim Bogaerts. This [project](http://f.dominec.eu/meep/) is now maintained by Filip Dominec.
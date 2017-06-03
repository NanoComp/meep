---
title: Meep Download
permalink: /Meep_Download/
---

The latest release of Meep is **version 1.3** (which requires [libctl](http://ab-initio.mit.edu/wiki/index.php/Libctl) version 3.2 or later), which may be downloaded from:

-   <http://ab-initio.mit.edu/meep/meep-1.3.tar.gz>

You can also download the latest development sources from [Meep on Github](https://github.com/stevengj/meep).

Older releases may be found at <http://ab-initio.mit.edu/meep/old>

Meep is [free software under the GNU GPL](Meep_License_and_Copyright.md) and comes with NO WARRANTY of any kind (see the [license](Meep_License_and_Copyright.md)).

Refer to the [Meep release notes](Meep_release_notes.md) to see what's new in this version, and be sure to read the [Meep Installation](Meep_Installation.md) manual for how to compile and install it.

Please subscribe to the `meep-announce` mailing list to receive a message when Meep is updated:

-   [meep-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/meep-announce)

Python interface
----------------

A Python interface for Meep is currently under development. A separate version was developed by the [Photonics Research Group](http://photonics.intec.ugent.be/) at Ghent University, in particular thanks to Emmanuel Lambert, Martin Fiers, Shavkat Nizamov, Martijn Tassaert, Peter Bienstman, and Wim Bogaerts. It can be downloaded from:

-   [Python-Meep](https://launchpad.net/python-meep)

Additional information including [installation instructions and tutorials](http://f.dominec.eu/meep/) is provided by F. Dominec.

Precompiled Meep packages for Debian and Ubuntu
-----------------------------------------------

A convenient precompiled package of Meep is also available in [Debian GNU/Linux](https://en.wikipedia.org/wiki/Debian_GNU/Linux) and [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu), specifically as the `meep` package (currently only in Debian `testing` or `unstable`):

-   <http://packages.debian.org/testing/science/meep>

We highly recommend using Debian or Ubuntu, as in Debian or Ubuntu the Meep software and all of its dependencies can be installed simply by typing one line:

```
apt-get install meep h5utils
```


You can also install the [parallel version of Meep](http://packages.debian.org/testing/science/meep-mpi) using:

```
apt-get install meep-mpi
```


Currently, the Debian `meep-mpi` package uses [MPICH](https://en.wikipedia.org/wiki/MPICH), so if you want to use [Open MPI](https://en.wikipedia.org/wiki/Open_MPI) then you need to build your own.

Meep on Amazon Web Services (AWS)
---------------------------------

The latest, pre-installed versions of Meep and MPB running on Ubuntu can also be accessed on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetuscloud.com/launchsims.html).

[![Latest Docs](https://readthedocs.org/projects/pip/badge/?version=latest)](http://meep.readthedocs.io/en/latest/Meep/)
[![Build Status](https://travis-ci.org/stevengj/meep.svg?branch=master)](https://travis-ci.org/stevengj/meep)

Meep (or MEEP) is a free finite-difference time-domain (FDTD)
simulation software package developed at MIT to model electromagnetic
systems.  You can download Meep at the
home page:

* http://ab-initio.mit.edu/meep/

and the latest documentation is available on [readthedocs](http://meep.readthedocs.io/en/latest/Meep/).

To compile directly from the git repository, you need to run
```
sh autogen.sh
make
```
in the cloned directory in order to generate the necessary files.  (You will need GNU autotools and SWIG installed, among other things.)

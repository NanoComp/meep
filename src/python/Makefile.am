# Copyright (C) 2005-2017 Massachusetts Institute of Technology  
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2, or (at your option)
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

# TODO: Remove hard coded paths
# TODO: Use Makefile.am
# TODO: Integrated with existing build system
# TODO: Py2 and Py3

MEEP_ROOT=/home/chris/meep
BUILD_ROOT=./build
PY_ROOT=/home/chris/miniconda2/envs/mp

all:
	$(SWIG) -python -c++ -Wextra -I/home/chris/dev/simpetus/meep/src meep.i
	$(CXX) -c -fPIC -I/home/chris/dev/simpetus/meep/src -I$(PY_ROOT)/include/python2.7 -I$(PY_ROOT)/lib/python2.7/site-packages/numpy/core/include -o _meep.o meep_wrap.cxx
	$(CXX) -shared -Wl,--no-as-needed -L$(MEEP_ROOT)/lib -L$(PY_ROOT)/lib -lpython2.7 -lmeep -o $(BUILD_ROOT)/meep/_meep.so _meep.o
	cp meep.py $(BUILD_ROOT)/meep/meep.py

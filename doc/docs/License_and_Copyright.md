---
# License and Copyright
---

Meep is copyright © 2005–2017, Massachusetts Institute of Technology.

Meep is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or at your option any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the GNU web site:

[http://www.gnu.org/copyleft/gpl.html](http://www.gnu.org/copyleft/gpl.html)

As a clarification, we should note that Scheme control (ctl) files, written by the user which do not contain code distributed with Meep and loaded at runtime by the Meep software, are *not* derived works of Meep and do *not* fall thereby under the restrictions of the GNU General Public License. On the other hand, C++ programs linked with the Meep libraries *are* derived works, and you must obey the terms of the GPL if you wish to distribute programs based in part on Meep. You are not affected for programs you do not distribute.

Referencing
-----------

We kindly request that you use the following reference in any publication for which you use Meep:

- A.F. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J.D. Joannopoulos, and S.G. Johnson, [MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method](http://dx.doi.org/doi:10.1016/j.cpc.2009.11.008), Computer Physics Communications, vol. 181, pp. 687-702 (2010).

If you want a one-sentence description of the algorithm for inclusion in a publication, we recommend something like:

- Simulations were performed with the finite-difference time-domain (FDTD) method [ref FDTD], using a freely available software package [ref Meep].

As a general reference on the FDTD method you might use, for example:

- A. Taflove and S.C. Hagness, Computational Electrodynamics: The Finite-Difference Time-Domain Method, Artech: Norwood, MA, (2005).

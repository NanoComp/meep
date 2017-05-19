/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/
#include <stdio.h>
#include <math.h>

#include "meep.hpp"
#include "mpSession.hpp"

using namespace meep;
using namespace meepSession;

int main(int argc, char *argv[])
{
  mpSession S;

  double s=11.0;           // size of computational cell, not including PML
  double dpml=1.0;         // thickness of PML layers
  double sxy = s+2.0*dpml; // cell size, including PML

  //(set! geometry-lattice (make lattice (size sxy sxy no-size)))
  S.set_geometry_lattice(sxy, sxy);

  //(set! pml-layers (list (make pml (thickness dpml))))
  //S.set_pml_layers( pml(dpml) );

  int resolution=10;
  double fcen = 0.8;
  double df   = 0.02;
  vec kdir(1.0,1.0);
  vec k = kdir*2.0*pi*fcen/abs(kdir);

  //S.add_source( 

  //S.add_step_func( 
  double T=400.0;
  S.run_until(T);

}

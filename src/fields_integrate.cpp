/* Copyright (C) 2004 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "meep.h"
#include "meep_internals.h"

namespace meep {

/* The following two functions convert a vec to the nearest ivec
   in the dielectric (odd-coordinate) grid, either rounding down (floor)
   or up (ceil).  In the special case where a component of the vec is
   *exactly* on a component of the ivec, we add the corresponding
   component of equal_shift (which should be either -2, 0, or +2).
   (equal_shift is there to prevent us from counting edge points twice.) */
   
static ivec vec2diel_floor(const vec &v, double a, const ivec &equal_shift) {
  ivec iv(v.dim);
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    iv.set_direction(d, 1+2*int(floor(v.in_direction(d)*a-.5)));
    if (iv.in_direction(d) == v.in_direction(d))
      iv.set_direction(d, iv.in_direction(d) + equal_shift.in_direction(d));
  }
  return iv;
}
static ivec vec2diel_ceil(const vec &v, double a, const ivec &equal_shift) {
  ivec iv(v.dim);
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    iv.set_direction(d, 1+2*int(ceil(v.in_direction(d)*a-.5)));
    if (iv.in_direction(d) == v.in_direction(d))
      iv.set_direction(d, iv.in_direction(d) + equal_shift.in_direction(d));
  }
  return iv;
}

/* Generic function for computing integrals of fields, and
   integral-like things, over a volume WHERE.  The job of this
   function is to call INTEGRAND() for each chunk that intersects
   WHERE, passing it the chunk, the range of integer coordinates to
   integrate, the integration weights for the boundary points, and the
   bloch phase shift, translational shift, and symmetry operation to
   transform the chunk to the actual integration location.  (N.B.
   we apply the symmetry first to the chunk, *then* the shift.)

   We also pass the integrand dV0 and dV1, such that the integration
   "volume" dV is dV0 + dV1 * iloopR, where iloopR is the loop variable
   (starting from 0 at the starting integer coord and incrementing by
   1) corresponding to the direction R.  Note that, in the
   LOOP_OVER_IVECS macro, iloopR corresponds to the loop variable
   loop_i2 in Dcyl (cylindrical coordinates).  In other coordinates,
   dV1 is 0.  Note also that by "volume" dV we mean the integration
   unit corresponding to the dimensionality of WHERE (e.g. an area
   if WHERE is 2d, etc.)

   In particular, the integration coordinates are on the odd-indexed
   "epsilon grid", which has the virtue that it is disjoint for each
   chunk and each chunk has enough information to interpolate all of
   its field components onto this grid without communication.  Another
   virtue of this grid is that it is invariant under all of our symmetry
   transformations, so we can uniquely decide which transformed chunk
   gets to integrate which grid point.

   The integration weights are chosen to correspond to integrating the
   linear interpolation of the function values from these grid points.

   For a simple example of an integrand routine, see the
   tests/integrate.cpp file. */
void fields::integrate(field_integrand integrand, void *integrand_data,
		       const geometric_volume &where)
{
  // Argument checks:
  if (where.dim != v.dim)
    abort("Invalid dimensions %d for WHERE in fields::integrate", where.dim);
  LOOP_OVER_DIRECTIONS(v.dim, d) 
    if (where.in_direction_max(d) - where.in_direction_min(d) >
        user_volume.boundary_location(High, d)
	- user_volume.boundary_location(Low, d)) 
      abort("Cannot handle integration width larger than cell width in %s direction!\n", direction_name(d));


  /* Find the corners (is and ie) of the smallest bounding box for
     where, on the grid of odd-coordinate ivecs (i.e. the
     "epsilon grid"). */
  ivec is(vec2diel_floor(where.get_min_corner(), v.a, zero_ivec(v.dim)));
  ivec ie(vec2diel_ceil(where.get_max_corner(), v.a, zero_ivec(v.dim)));
  vec s0(v.dim), e0(v.dim), s1(v.dim), e1(v.dim);
  
  /* Integration weights at boundaries (c.f. long comment at bottom). */
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    double w0, w1;
    w0 = 1. - where.in_direction_min(d)*v.a + 0.5*is.in_direction(d);
    w1 = 1. + where.in_direction_max(d)*v.a - 0.5*ie.in_direction(d);
    if (ie.in_direction(d) >= is.in_direction(d) + 3*2) {
      s0.set_direction(d, w0*w0 / 2);
      s1.set_direction(d, 1 - (1-w0)*(1-w0) / 2);
      e0.set_direction(d, w1*w1 / 2);
      e1.set_direction(d, 1 - (1-w1)*(1-w1) / 2);
    }
    else if (ie.in_direction(d) == is.in_direction(d) + 2*2) {
      s0.set_direction(d, w0*w0 / 2);
      s1.set_direction(d, 1 - (1-w0)*(1-w0) / 2 - (1-w1)*(1-w1) / 2);
      e0.set_direction(d, w1*w1 / 2);
      e1.set_direction(d, s1.in_direction(d));
    }
    else if (where.in_direction_min(d) == where.in_direction_max(d)) {
      s0.set_direction(d, w0);
      s1.set_direction(d, w1);
      e0.set_direction(d, w1);
      e1.set_direction(d, w0);
    }
    else if (ie.in_direction(d) == is.in_direction(d) + 1*2) {
      s0.set_direction(d, w0*w0 / 2 - (1-w1)*(1-w1) / 2);
      e0.set_direction(d, w1*w1 / 2 - (1-w0)*(1-w0) / 2);
      s1.set_direction(d, e0.in_direction(d));
      e1.set_direction(d, s0.in_direction(d));
    }
    else
      abort("bug: impossible(?) integration boundaries");
  }


  // loop over symmetry transformations of the chunks:
  for (int sn = 0; sn < S.multiplicity(); ++sn) {
    geometric_volume vS = S.transform(v.surroundings(), sn);
    vec L(v.dim);
    ivec iL(v.dim);

    // n.b. we can't just S.transform(lattice_vector,sn), 'cause of origin
    LOOP_OVER_DIRECTIONS(v.dim, d) {
      direction dS = S.transform(d, -sn).d;
      L.set_direction(d, lattice_vector(dS).in_direction(dS));
      iL.set_direction(d, ilattice_vector(dS).in_direction(dS));
    }

    // figure out range of lattice shifts for which vS intersects where:
    ivec min_ishift(v.dim), max_ishift(v.dim);
    LOOP_OVER_DIRECTIONS(v.dim, d) {
      if (boundaries[High][S.transform(d, -sn).d] == Periodic) {
	min_ishift.set_direction(d, 
	 int(floor((where.in_direction_min(d) - vS.in_direction_max(d))
		   / L.in_direction(d))));
	max_ishift.set_direction(d,
	 int(ceil((where.in_direction_max(d) - vS.in_direction_min(d))
		  / L.in_direction(d))));
      }
      else {
	min_ishift.set_direction(d, 0);
	max_ishift.set_direction(d, 0);
      }
    }
    
    // loop over lattice shifts
    ivec ishift(min_ishift);
    int ishiftn = 0;
    do {
      complex<double> ph = 1.0;
      vec shift(v.dim);
      ivec shifti(v.dim);
      LOOP_OVER_DIRECTIONS(v.dim, d) {
	shift.set_direction(d, L.in_direction(d) * ishift.in_direction(d));
	shifti.set_direction(d, iL.in_direction(d) * ishift.in_direction(d));
	ph *= pow(eikna[d], ishift.in_direction(d));
      }

      for (int i = 0; i < num_chunks; ++i) {
	if (!chunks[i]->is_mine()) continue;
	// Chunk integration boundaries:
	geometric_volume gvS = S.transform(chunks[i]->gv, sn);
	ivec iscS(max(is-shifti, vec2diel_ceil(gvS.get_min_corner(),
					       v.a, one_ivec(v.dim) * 2)));
	ivec iecS(min(ie-shifti, vec2diel_floor(gvS.get_max_corner(),
						v.a, zero_ivec(v.dim))));
	if (iscS <= iecS) {
	  // Determine weights at chunk integration boundaries:
	  ivec isc(S.transform(iscS, -sn)), iec(S.transform(iecS, -sn));
	  vec s0c(v.dim), s1c(v.dim), e0c(v.dim), e1c(v.dim);
	  iscS += shifti;
	  iecS += shifti;
	  LOOP_OVER_DIRECTIONS(v.dim, d) {
	    direction dS = S.transform(d, sn).d;
	    if (iscS.in_direction(dS) == is.in_direction(dS)) {
	      s0c.set_direction(d, s0.in_direction(dS));
	      s1c.set_direction(d, s1.in_direction(dS));
	    }
	    else if (iscS.in_direction(dS) == is.in_direction(dS) + 2) {
	      s0c.set_direction(d, s1.in_direction(dS));
	      s1c.set_direction(d, 1.0);
	    }
	    else {
	      s0c.set_direction(d, 1.0);
	      s1c.set_direction(d, 1.0);
	    }
	    if (iecS.in_direction(dS) == ie.in_direction(dS)) {
	      e0c.set_direction(d, e0.in_direction(dS));
	      e1c.set_direction(d, e1.in_direction(dS));
	    }
	    else if (iecS.in_direction(dS) == ie.in_direction(dS) - 2) {
	      e0c.set_direction(d, e1.in_direction(dS));
	      e1c.set_direction(d, 1.0);
	    }
	    else {
	      e0c.set_direction(d, 1.0);
	      e1c.set_direction(d, 1.0);
	    }
	    if (iecS.in_direction(dS) == iscS.in_direction(dS)) {
	      double w = min(s0c.in_direction(d), e0c.in_direction(d));
	      s0c.set_direction(d, w); e0c.set_direction(d, w);
	      s1c.set_direction(d, w); e1c.set_direction(d, w);
	    }
	    else if (iecS.in_direction(dS) == iscS.in_direction(dS) + 1*2) {
	      double w = min(s0c.in_direction(d), e1c.in_direction(d));
	      s0c.set_direction(d, w); e1c.set_direction(d, w);
	      w = min(s1c.in_direction(d), e0c.in_direction(d));
	      s1c.set_direction(d, w); e0c.set_direction(d, w);
	    }
	    else if (iecS.in_direction(dS) == iscS.in_direction(dS) + 2*2) {
	      double w = min(s1c.in_direction(d), e1c.in_direction(d));
	      s1c.set_direction(d, w); e1c.set_direction(d, w);
	    }

	    // swap endpoints/weights if in wrong order due to S.transform
	    if (isc.in_direction(d) > iec.in_direction(d)) {
	      int iswap = isc.in_direction(d);
	      isc.set_direction(d, iec.in_direction(d));
	      iec.set_direction(d, iswap);
	      double swap = s0c.in_direction(d);
	      s0c.set_direction(d, e0c.in_direction(d));
	      e0c.set_direction(d, swap);
	      swap = s1c.in_direction(d);
	      s1c.set_direction(d, e1c.in_direction(d));
	      e1c.set_direction(d, swap);
	    }
	  }
	  

	  // Determine integration "volumes" dV0 and dV1;
	  double dV0 = 1.0, dV1 = 0.0;
	  LOOP_OVER_DIRECTIONS(v.dim, d)
	    if (where.in_direction_max(d) > where.in_direction_min(d))
	      dV0 *= v.inva;
	  if (v.dim == Dcyl) {
	    dV0 *= 2*pi * (S.transform(chunks[i]->v[isc], sn) + shift)
	      .in_direction(R);
	    dV1 = 2*pi * S.transform(chunks[i]->v[unit_ivec(v.dim,R)*2], sn)
	      .in_direction(R);
	  }
	 
	  integrand(chunks[i], 
		    isc, iec,
		    s0c, s1c, e0c, e1c,
		    dV0, dV1,
		    shift, ph,
		    S, sn,
		    integrand_data);
	}
      }
      
      
      ++ishiftn;
      LOOP_OVER_DIRECTIONS(v.dim, d) {
	if (ishift.in_direction(d) + 1 <= max_ishift.in_direction(d)) {
	  ishift.set_direction(d, ishift.in_direction(d) + 1);
	  break;
	}
	ishift.set_direction(d, min_ishift.in_direction(d));
      }
    } while (ishift != min_ishift);
  }
}

/*

                     Integration Weights

We want the integral from a to b, assuming linear interpolation of fn
(function values on grid points n).  Most interior points have weight
1, but the points just inside and outside the boundaries have
different weights.  Call the weights for the points just *outside* the
starting and ending boundaries s0 and e0, respectively, and weights
for the points just *inside* the boundaries s1 and e1.  Then we have
to handle the following cases:

1) a and b separated by at least 2 grid points, e.g.:

x  |    x       x       x  |    x
0  a    1       2       3  b    4

first segment: f(x) = f0 (1 - x) + f1 x
      -- \int_a^1 f(x) dx = f0 (1 - a)^2/2 + f1 (1 - a^2) / 2

last segment: f(x) = f3 (4 - x) + f4 (x - 3)
     -- \int_3^b f(x) dx = f3 [1 - (4-b)^2] / 2 + f4 (b - 3)^2 / 2

integral = f0 (1 - a)^2/2         <---- f0 s0
         + f1 (1 - a^2/2)         <---- f1 s1
         + f2
         + f3 (1 - (4-b)^2 / 2)   <---- f3 e1
         + f4 (b - 3)^2 / 2       <---- f4 e0

In terms of starting and ending weights:
	 w0 = 1 - a
	 w1 = b - 3

	 s0 = w0^2 / 2
	 s1 = 1 - (1 - w0)^2 / 2
	 e0 = w1^2 / 2
	 e1 = 1 - (1 - w1)^2 / 2

2) one grid point between a and b.

x  |    x  |    x
0  a    1  b    2

integral = f0 (1 - a)^2 / 2
         + f1 [(1 - a^2) + 1 - (2 - b)^2] / 2
         + f3 (b - 1)^2 / 2

s0 = w0^2 / 2
e0 = w1^2 / 2
s1' = e1' = 1 - (1 - w0)^2 / 2 - (1 - w1)^2 / 2

3) no grid points between a and b.

x  |      |    x
0  a      b    1

integral = f0 [ (1-a)^2 - (1-b)^2 ] / 2 + f1 [ b^2 - a^2 ] / 2
         = f0 [ w0^2 - (1-w1)^2 ]  / 2 + f1 [ w1^2 - (1-w0)^2 ] / 2

s0 = e1 = w0^2/2 - (1-w1)^2/2
e0' = s1' = w1^2/2 - (1-w0)^2/2

4) as (3), but a = b: interpolation, not integration:

 -- want: f0 * w0 + f1 * w1

s0' = w0
e0' = w1 = 1 - w0

*/

} // namespace meep

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

/****************************************************************************

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
s1 = e1 = 1 - (1 - w0)^2 / 2 - (1 - w1)^2 / 2

3) no grid points between a and b.

x  |      |    x
0  a      b    1

integral = f0 [ (1-a)^2 - (1-b)^2 ] / 2 + f1 [ b^2 - a^2 ] / 2
         = f0 [ w0^2 - (1-w1)^2 ]  / 2 + f1 [ w1^2 - (1-w0)^2 ] / 2

s0 = e1 = w0^2/2 - (1-w1)^2/2
e0 = s1 = w1^2/2 - (1-w0)^2/2

4) as (3), but a = b: interpolation, not integration:

 -- want: f0 * w0 + f1 * w1

s0 = w0
e0 = w1 = 1 - w0

                        --------------

          Integration Weights in Cylindrical Coordinates
                   FIXME: implement this below?

Ideally, we should have different weights for the R direction of
cylindrical coordinates, i.e. for integrating f(r) r dr, because again
we want to perfectly integrate any linear f(r).  Thus, the integration
weights will depend upon r.  Note, however, that we also have an r in
the dV, so we will have to divide the weights by this factor.

1) a and b separated by at least 2 grid points, e.g.:

x  |    x       x       x  |    x
i  a   i+1     i+2     i+3  b  i+4

(where r = i * inva).

linear interpolation in [i,i+1): f(x) = f_i (i+1 - x) + f_{i+1} (x-i)

   want: \int_a^b f(x) x dx

in terms of starting and ending weights:
	 w0 = (i+1) - a
	 w1 = b - (i+3)

integral = f_i [-w0^3 / 3 + (i+1) w0^2 / 2]                       <- s0 i
      + f_{j=i+1} [w0^3 / 3 - (j+1) w0^2 / 2 + j w0 + j/2 + 1/6]  <- s1 (i+1)
      + f_{j=i+2} j                                               <- 1 (i+2)
      + f_{j=i+3} [-w1^3 / 3 - (j-1) w1^2 / 2 + j w1 + j/2 - 1/6] <- e1 (i+3)
      + f_{j=i+4} [w1^3 / 3 + (j-1) w1^2 / 2]                     <- e0 (i+4)

      (thanks to Maple for doing the annoying algebra)
  (yes, I have tested that it correctly integrates linear f(r))

Note that the coefficients need to be divided by i, i+1, etcetera to
get s0, s1, etcetera; this gives an interior-point weight of 1 as
before.  For i->infinity, this should converge to the weights from
before.  Avoiding division by zero is more tricky, because the weight
at j=0 is not necessarily zero, due to the interpolation.  It might be
better to pre-include the dV in the weight for edge elements, with
appropriate logic in the IVEC_LOOP_WEIGHT macro.  Tricky.

The above is also not correct for integrals that cross x=0, because
it should really be the integral of f(x) |x|.  Even interior points
probably need special handling in that case.  For sanity, we would
just divide the integration region into positive and negative r and
integrate them separately somehow.  Grrr.

2) one grid point between a and b.

x  |    x   |    x
i  a   i+1  b   i+2

integral = f_i [-w0^3 / 3 + (i+1) w0^2 / 2]                 <- s0 i
      + f_{j=i+1} [w0^3 / 3 - (j+1) w0^2 / 2 + j w0 +
                  -w1^3 / 3 - (j-1) w1^2 / 2 + j w1]        <- {s1,e1} (i+1)
      + f_{j=i+2} [w1^3 / 3 + (j-1) w1^2 / 2]               <- e0 (i+2)

3) no grid points between a and b.

x  |      |    x
i  a      b   i+1

integral = f_i [-w0^3/3 + (i+1) w0^2/2
              + -w1^3/3 - (i-1) w1^2/2 + i w1 - i/2 - 1/6]  <- s0 i
   + f_{j=i+1} [ w0^3/3 - (j+1) w0^2/2 + j w0
              +  w1^3/3 + (j-1) w1^2/2 - j/2 + 1/6]         <- e0 (i+1)

4) as (3), but a = b: interpolation, not integration: same as above

****************************************************************************/

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

static inline int iabs(int i) { return (i < 0 ? -i : i); }

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

   In particular, the integration coordinates are calculated on the
   Yee grid for component cgrid.  cgrid == Dielectric is a good choice
   if you want to integrate a combination of multiple field components,
   because all of the field components can be interpolated onto this
   grid without communication between chunks.

   The integration weights are chosen to correspond to integrating the
   linear interpolation of the function values from these grid points.

   For a simple example of an integrand routine, see the
   tests/integrate.cpp file. */
void fields::integrate(field_integrand integrand, void *integrand_data,
		       const geometric_volume &where, 
		       component cgrid,
		       bool use_symmetry)
{
  if (coordinate_mismatch(v.dim, component_direction(cgrid)))
    abort("Invalid fields::integrate grid type %s for dimensions %s\n",
	  component_name(cgrid), dimension_name(v.dim));
  if (where.dim != v.dim)
    abort("Invalid dimensions %d for WHERE in fields::integrate", where.dim);
  
  /*
    We handle integration on an arbitrary component grid by shifting
    to the dielectric grid and then shifting back.  The integration
    coordinates are internally calculated on the odd-indexed
    "dielectric grid", which has the virtue that it is disjoint for
    each chunk and each chunk has enough information to interpolate all
    of its field components onto this grid without communication.
    Another virtue of this grid is that it is invariant under all of
    our symmetry transformations, so we can uniquely decide which
    transformed chunk gets to integrate which grid point. 
  */
  vec yee_c(v.yee_shift(Dielectric) - v.yee_shift(cgrid));
  ivec iyee_c(v.iyee_shift(Dielectric) - v.iyee_shift(cgrid));
  geometric_volume wherec(where + yee_c);


  /* Find the corners (is and ie) of the smallest bounding box for
     wherec, on the grid of odd-coordinate ivecs (i.e. the
     "epsilon grid"). */
  ivec is(vec2diel_floor(wherec.get_min_corner(), v.a, zero_ivec(v.dim)));
  ivec ie(vec2diel_ceil(wherec.get_max_corner(), v.a, zero_ivec(v.dim)));
  
  /* Integration weights at boundaries (c.f. long comment at top). */
  vec s0(v.dim), e0(v.dim), s1(v.dim), e1(v.dim);
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    double w0, w1;
    w0 = 1. - wherec.in_direction_min(d)*v.a + 0.5*is.in_direction(d);
    w1 = 1. + wherec.in_direction_max(d)*v.a - 0.5*ie.in_direction(d);
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
    else if (wherec.in_direction_min(d) == wherec.in_direction_max(d)) {
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
  for (int sn = 0; sn < (use_symmetry ? S.multiplicity() : 1); ++sn) {
    component cS = S.transform(cgrid, -sn);
    ivec iyee_cS(S.transform_unshifted(iyee_c, -sn));

    geometric_volume vS = S.transform(v.surroundings(), sn);
    vec L(v.dim);
    ivec iL(v.dim);

    // n.b. we can't just S.transform(lattice_vector,sn), 'cause of signs
    LOOP_OVER_DIRECTIONS(v.dim, d) {
      direction dS = S.transform(d, -sn).d;
      L.set_direction(d, fabs(lattice_vector(dS).in_direction(dS)));
      iL.set_direction(d, iabs(ilattice_vector(dS).in_direction(dS)));
    }

    // figure out range of lattice shifts for which vS intersects wherec:
    ivec min_ishift(v.dim), max_ishift(v.dim);
    LOOP_OVER_DIRECTIONS(v.dim, d) {
      if (boundaries[High][S.transform(d, -sn).d] == Periodic) {
	min_ishift.set_direction(d, 
	 int(floor((wherec.in_direction_min(d) - vS.in_direction_max(d))
		   / L.in_direction(d))));
	max_ishift.set_direction(d,
	 int(ceil((wherec.in_direction_max(d) - vS.in_direction_min(d))
		  / L.in_direction(d))));
      }
      else {
	min_ishift.set_direction(d, 0);
	max_ishift.set_direction(d, 0);
      }
    }
    
    // loop over lattice shifts
    ivec ishift(min_ishift);
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
	geometric_volume gvS(v.dim);

	if (use_symmetry)
	  gvS = S.transform(chunks[i]->gv, sn);
	else {
	  /* If we're not using symmetry, it's because (as in src_vol)
	     we don't care about correctly counting the points in the
	     volume.  Rather, we just want to make sure to get *all*
	     of the chunk points that intersect where.  Hence, add a little
	     padding to make sure we don't miss any points due to rounding. */
	  vec pad(one_ivec(v.dim) * v.a * 1e-3);
	  gvS = geometric_volume(chunks[i]->v.loc(Dielectric,0) - pad,
				 chunks[i]->v.loc(Dielectric,
						  chunks[i]->v.ntot()-1) +pad);
	}

	ivec iscS(max(is-shifti, vec2diel_ceil(gvS.get_min_corner(),
					       v.a, one_ivec(v.dim) * 2)));
	ivec iecS(min(ie-shifti, vec2diel_floor(gvS.get_max_corner(),
						v.a, zero_ivec(v.dim))));
	if (iscS <= iecS) {
	  // Determine weights at chunk integration boundaries:
	  ivec isc(S.transform(iscS, -sn)), iec(S.transform(iecS, -sn));
	  vec s0c(v.dim,1.0), s1c(v.dim,1.0), e0c(v.dim,1.0), e1c(v.dim,1.0);
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
	    }
	    if (iecS.in_direction(dS) == ie.in_direction(dS)) {
	      e0c.set_direction(d, e0.in_direction(dS));
	      e1c.set_direction(d, e1.in_direction(dS));
	    }
	    else if (iecS.in_direction(dS) == ie.in_direction(dS) - 2) {
	      e0c.set_direction(d, e1.in_direction(dS));
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
	    if (wherec.in_direction_max(d) > wherec.in_direction_min(d))
	      dV0 *= v.inva;
	  if (v.dim == Dcyl) {
	    dV1 = dV0 * 2*pi * v.inva;
	    dV0 *= 2*pi * fabs((S.transform(chunks[i]->v[isc], sn) + shift
				- yee_c).in_direction(R));
	  }
	 
	  integrand(chunks[i], cS,
		    isc - iyee_cS, iec - iyee_cS,
		    s0c, s1c, e0c, e1c,
		    dV0, dV1,
		    shift, ph,
		    S, sn,
		    integrand_data);
	}
      }
      
      
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

} // namespace meep

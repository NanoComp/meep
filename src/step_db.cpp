/* Copyright (C) 2005-2008 Massachusetts Institute of Technology
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
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"

#define RESTRICT

namespace meep {

void fields::step_db(field_type ft) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_db(ft);
}

void fields_chunk::step_db(field_type ft) {
  if (ft != B_stuff && ft != D_stuff)
    abort("bug - step_db should only be called for B or D");

  bool have_pml = false;
  FOR_FT_COMPONENTS(ft, cc)
    if (s->sigsize[cycle_direction(v.dim,component_direction(cc),1)] > 1)
	have_pml = true;
  if (have_pml) FOR_FT_COMPONENTS(ft, cc) DOCMP
    if (f[cc][cmp]) {
      if (!f_prev[cc][cmp])
        f_prev[cc][cmp] = new double[v.ntot()];
      // update_e_from_d requires previous D-P, not D, if D-P is present
      memcpy(f_prev[cc][cmp],
	     f_minus_p[cc][cmp] ? f_minus_p[cc][cmp] : f[cc][cmp],
	     v.ntot()*sizeof(double));
    }

  DOCMP FOR_FT_COMPONENTS(ft, cc)
    if (f[cc][cmp]) {
      const component c_p=plus_component[cc], c_m=minus_component[cc];
      const direction d_deriv_p = plus_deriv_direction[cc];
      const direction d_deriv_m = minus_deriv_direction[cc];
      const direction d_c = component_direction(cc);
      const bool have_p = have_plus_deriv[cc];
      const bool have_m = have_minus_deriv[cc];
      const direction dsig0 = cycle_direction(v.dim,d_c,1);
      const bool have_pml = s->sigsize[dsig0] > 1;
      const direction dsig = have_pml ? dsig0 : NO_DIRECTION;
      int stride_p = have_p?v.stride(d_deriv_p):0;
      int stride_m = have_m?v.stride(d_deriv_m):0;
      double *f_p = have_p?f[c_p][cmp]:NULL;
      double *f_m = have_m?f[c_m][cmp]:NULL;
      double *the_f = f[cc][cmp];
      
      if (ft == D_stuff) { // strides are opposite sign for H curl
	stride_p = -stride_p;
	stride_m = -stride_m;
      }

      if (v.dim == Dcyl) switch (d_c) {
      case R:
	f_p = NULL; // im/r Fz term will be handled separately
	break;
      case P:
	break; // curl works normally for phi component
      case Z: {
	f_m = NULL; // im/r Fr term will be handled separately
	
	/* Here we do a somewhat cool hack: the update of the z
	   component gives a 1/r d(r Fp)/dr term, rather than
	   just the derivative dg/dr expected in step_curl.
	   Rather than duplicating all of step_curl to handle
	   this bloody derivative, however, we define a new
	   array f_rderiv_int which is the integral of 1/r d(r Fp)/dr,
	   so that we can pass it to the unmodified step_curl
	   and get the correct derivative.  (More precisely,
	   the derivative and integral are replaced by differences
	   and sums, but you get the idea). */
	if (!f_rderiv_int) f_rderiv_int = new double[v.ntot()];
	double ir0 = (v.origin_r() + rshift) * v.a 
	  + 0.5 * v.iyee_shift(c_p).in_direction(R);
	for (int iz = 0; iz <= v.nz(); ++iz) f_rderiv_int[iz] = 0;
	int sr = v.nz() + 1;
	for (int ir = 1; ir <= v.nr(); ++ir) {
	  double rinv = 1.0 / ((ir+ir0)-0.5);
	  for (int iz = 0; iz <= v.nz(); ++iz) {
	    int idx = ir*sr + iz;
	    f_rderiv_int[idx] = f_rderiv_int[idx - sr] +
	      rinv * (f_p[idx] * (ir+ir0) - f_p[idx - sr] * ((ir-1)+ir0));
	  }
	}
	f_p = f_rderiv_int;
	break;
      }
      default: abort("bug - non-cylindrical field component in Dcyl");
      }
      
      step_curl(the_f, cc, f_p, f_m, stride_p, stride_m, v, Courant, 
		dsig, s->sig[dsig], s->siginv[dsig],
		dt, s->conductivity[cc][d_c], s->condinv[cc][d_c]);
    }

  // in cylindrical coordinates, we now have to add the i*m/r terms... */
  if (v.dim == Dcyl && m != 0) DOCMP FOR_FT_COMPONENTS(ft, cc) {
    const direction d_c = component_direction(cc);
    if (f[cc][cmp] && (d_c == R || d_c == Z)) {
      const component c_g = d_c==R ? plus_component[cc] : minus_component[cc];
      const double *g = f[c_g][1-cmp];
      double *the_f = f[cc][cmp];
      const double *cndinv = s->condinv[cc][d_c];
      const direction dsig = cycle_direction(v.dim,d_c,1);
      const double the_m = 
	m * (1-2*cmp) * (1-2*(ft==B_stuff)) * (1-2*(d_c==R)) * Courant;
      const double ir0 = (v.origin_r() + rshift) * v.a 
	+ 0.5 * v.iyee_shift(cc).in_direction(R);
      int sr = v.nz() + 1;
      if (cndinv) { // conductivity, possibly including PML
	for (int ir = ir0 <= 0.5; ir <= v.nr(); ++ir) {
	  double rinv = the_m / (ir+ir0);
	  for (int iz = 0; iz <= v.nz(); ++iz) {
	    int idx = ir*sr + iz;
	    the_f[idx] += rinv * g[idx] * cndinv[idx];
	  }
	}
      }
      else if (s->sigsize[dsig] > 1) { // PML
	const double *siginv = s->siginv[dsig];
	int dk = v.iyee_shift(cc).in_direction(dsig);
	for (int ir = ir0 <= 0.5; ir <= v.nr(); ++ir) {
	  double rinv = the_m / (ir+ir0);
	  for (int iz = 0; iz <= v.nz(); ++iz) {
	    int idx = ir*sr + iz;
	    the_f[idx] += rinv * g[idx] * siginv[dk + 2*(dsig==Z ? iz : ir)];
	  }
	}
      }
      else { // no PML, no conductivity
	for (int ir = ir0 <= 0.5; ir <= v.nr(); ++ir) {
	  double rinv = the_m / (ir+ir0);
	  for (int iz = 0; iz <= v.nz(); ++iz) {
	    int idx = ir*sr + iz;
	    the_f[idx] += rinv * g[idx];
	  }
	}
      }
    }
  }

  // deal with annoying r=0 boundary conditions for m=0 and m=1
  if (v.dim == Dcyl && v.origin_r() == 0.0) DOCMP {
    if (m == 0 && ft == D_stuff && f[Dz][cmp]) {
      // d(Dz)/dt = (1/r) * d(r*Hp)/dr
      double *the_f = f[Dz][cmp];
      const double *g = f[Hp][cmp];
      const double *cndinv = s->condinv[Dz][Z];
      const direction dsig = cycle_direction(v.dim,Z,1);
      if (cndinv) // conductivity, possibly including PML
	for (int iz = 0; iz < v.nz(); ++iz) 
	  the_f[iz] += g[iz] * (Courant * 4) * cndinv[iz];
      else if (s->sigsize[dsig] > 1) { // PML
	const double *siginv = s->siginv[dsig];
	int dk = v.iyee_shift(Dz).in_direction(dsig);
	for (int iz = 0; iz < v.nz(); ++iz) 
	  the_f[iz] += g[iz] * (Courant * 4) * siginv[dk + 2*(dsig==Z)*iz];
      }
      else // no PML, no conductivity
	for (int iz = 0; iz < v.nz(); ++iz) 
	  the_f[iz] += g[iz] * (Courant * 4);
      // Note: old code was missing factor of 4??

      for (int iz = 0; iz <= v.nz(); ++iz) f[Dp][cmp][iz] = 0.0;
    }
    else if (m == 0 && ft == B_stuff && f[Br][cmp])
      for (int iz = 0; iz <= v.nz(); ++iz) f[Br][cmp][iz] = 0.0;
    else if (fabs(m) == 1) {
      // D_stuff: d(Dp)/dt = d(Hr)/dz - d(Hz)/dr
      // B_stuff: d(Br)/dt = d(Ep)/dz - i*m*Ez/r
      component cc = ft == D_stuff ? Dp : Br;
      direction d_c = component_direction(cc);
      double *the_f = f[cc][cmp];
      if (!the_f) continue;
      const double *f_p = f[ft == D_stuff ? Hr : Ep][cmp];
      const double *f_m = ft == D_stuff ? f[Hz][cmp]
	: (f[Ez][1-cmp] + (v.nz()+1));
      const double *cndinv = s->condinv[cc][d_c];
      const direction dsig = cycle_direction(v.dim,d_c,1);
      int sd = ft == D_stuff ? +1 : -1;
      double f_m_mult = ft == D_stuff ? 2 : (1-2*cmp);
      if (cndinv) // conductivity, possibly including PML
	for (int iz = (ft == D_stuff); iz < v.nz() + (ft == D_stuff); ++iz)
	  the_f[iz] += (sd*Courant) * (f_p[iz]-f_p[iz-sd] - f_m_mult*f_m[iz])
	    * cndinv[iz];
      else if (s->sigsize[dsig] > 1) { // PML
	const double *siginv = s->siginv[dsig];
	int dk = v.iyee_shift(cc).in_direction(dsig);
	for (int iz = (ft == D_stuff); iz < v.nz() + (ft == D_stuff); ++iz)
	  the_f[iz] += (sd*Courant) * (f_p[iz]-f_p[iz-sd] - f_m_mult*f_m[iz])
	    * siginv[dk + 2*(dsig==Z)*iz];
      }
      else // no PML, no conductivity
	for (int iz = (ft == D_stuff); iz < v.nz() + (ft == D_stuff); ++iz)
	  the_f[iz] += (sd*Courant) * (f_p[iz]-f_p[iz-sd] - f_m_mult*f_m[iz]);

      if (ft == D_stuff)
	for (int iz = 0; iz <= v.nz(); ++iz) f[Dz][cmp][iz] = 0.0;
    }
    else if (m != 0) { // m != {0,+1,-1}
      /* I seem to recall David telling me that this was for numerical
	 stability of some sort - the larger m is, the farther from
	 the origin we need to be before we can use nonzero fields
	 ... note that this is a fixed number of pixels for a given m,
	 so it should still converge.  Still, this is weird... */
      double rmax = fabs(m) - int(v.origin_r()*v.a+0.5);
      if (ft == D_stuff)
	for (int r = 0; r <= v.nr() && r < rmax; r++) {
          const int ir = r*(v.nz()+1);
          for (int z=0;z<=v.nz();z++) f[Dp][cmp][ir+z] = f[Dz][cmp][ir+z] = 0;
        }
      else
	for (int r = 0; r <= v.nr() && r < rmax; r++) {
          const int ir = r*(v.nz()+1);
          for (int z=0;z<=v.nz();z++) f[Br][cmp][ir+z] = 0;
        }
    }
  }
}

} // namespace meep

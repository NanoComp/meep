/* Copyright (C) 2006 Massachusetts Institute of Technology
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

#include "meep.hpp"
#include "meep_internals.hpp"

#define RESTRICT

namespace meep {

static inline double it(int cmp, double *(f[2]), int ind) {
  return (f[1-cmp]) ? (1-2*cmp)*f[1-cmp][ind] : 0;
}

inline int rstart_0(const volume &v, double m) {
  return (int) max(0.0, m - (int)(v.origin_r()*v.a+0.5) - 1.0);
}
inline int rstart_1(const volume &v, double m) {
  return (int) max(1.0, m - (int)(v.origin_r()*v.a+0.5));
}

void fields::step_h() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_h();
}

void fields_chunk::step_h() {
  const volume v = this->v;
  if (v.dim != Dcyl) {
    const int n1 = num_each_direction[0];
    const int n2 = num_each_direction[1];
    const int n3 = num_each_direction[2];
    const int s1 = stride_each_direction[0];
    const int s2 = stride_each_direction[1];
    const int s3 = stride_each_direction[2];
    DOCMP FOR_MAGNETIC_COMPONENTS(cc)
      if (f[cc][cmp]) {
        const int yee_idx = v.yee_index(cc);
        const component c_p=plus_component[cc], c_m=minus_component[cc];
        const direction d_deriv_p = plus_deriv_direction[cc];
        const direction d_deriv_m = minus_deriv_direction[cc];
        const bool have_p = have_plus_deriv[cc];
        const bool have_m = have_minus_deriv[cc];
        const bool have_p_pml = have_p && s->C[d_deriv_p][cc];
        const bool have_m_pml = have_m && s->C[d_deriv_m][cc];
        const int stride_p = (have_p)?v.stride(d_deriv_p):0;
        const int stride_m = (have_m)?v.stride(d_deriv_m):0;
        // The following lines "promise" the compiler that the values of
        // these arrays won't change during this loop.
        RESTRICT const double *C_m = (have_m_pml)?s->C[d_deriv_m][cc] + yee_idx:NULL;
        RESTRICT const double *C_p = (have_p_pml)?s->C[d_deriv_p][cc] + yee_idx:NULL;
        RESTRICT const double *decay_m = (!have_m_pml)?NULL:
          s->Cdecay[d_deriv_m][cc][component_direction(cc)] + yee_idx;
        RESTRICT const double *decay_p = (!have_p_pml)?NULL:
          s->Cdecay[d_deriv_p][cc][component_direction(cc)] + yee_idx;
        RESTRICT const double *f_p = (have_p)?f[c_p][cmp] + v.yee_index(c_p):NULL;
        RESTRICT const double *f_m = (have_m)?f[c_m][cmp] + v.yee_index(c_m):NULL;
        RESTRICT double *the_f = f[cc][cmp] + yee_idx;
        RESTRICT double *the_f_p_pml = f_p_pml[cc][cmp] + yee_idx;
        RESTRICT double *the_f_m_pml = f_m_pml[cc][cmp] + yee_idx;
#include "step_h.hpp"
      }
  } else if (v.dim == Dcyl) {
    int ir0 = int((v.origin_r() + rshift) * v.a + 0.5);
    DOCMP {
      // Propogate Hr
      if (s->C[Z][Hr])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/((int)(v.origin_r()*v.a+0.5) + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Czhr = s->C[Z][Hr][z+ir];
              const double ooop_Czhr = s->Cdecay[Z][Hr][R][z+ir];
              double dhrp = Courant*(-it(cmp,f[Ez],z+ir)*mor);
              double hrz = f[Hr][cmp][z+ir] - f_p_pml[Hr][cmp][z+ir];
              f_p_pml[Hr][cmp][z+ir] += dhrp;
              f[Hr][cmp][z+ir] += dhrp +
                ooop_Czhr*(Courant*(f[Ep][cmp][z+ir+1]-f[Ep][cmp][z+ir]) - Czhr*hrz);
            }
          }
      else
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/((int)(v.origin_r()*v.a + 0.5) + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Hr][cmp][z+ir] += Courant*
                ((f[Ep][cmp][z+ir+1]-f[Ep][cmp][z+ir]) - it(cmp,f[Ez],z+ir)*mor);
          }
      // Propogate Hp
      if (s->C[Z][Hp] || s->C[R][Hp])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Czhp = (s->C[Z][Hp])?s->C[Z][Hp][z+ir]:0;
              const double Crhp = (s->C[R][Hp])?s->C[R][Hp][z+ir]:0;
              const double ooop_Czhp = (s->Cdecay[Z][Hp][P]) ?
                s->Cdecay[Z][Hp][P][z+ir]:1.0;
              const double dhpz = ooop_Czhp*(-Courant*(f[Er][cmp][z+ir+1]-f[Er][cmp][z+ir])
                                             - Czhp*f_p_pml[Hp][cmp][z+ir]);
              const double hpr = f[Hp][cmp][z+ir]-f_p_pml[Hp][cmp][z+ir];
              f_p_pml[Hp][cmp][z+ir] += dhpz;
              f[Hp][cmp][z+ir] += dhpz +
                ooop_Czhp*(Courant*(f[Ez][cmp][z+irp1]-f[Ez][cmp][z+ir]) - Crhp*hpr);
            }
          }
      else 
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Hp][cmp][z+ir] += Courant*
                ((f[Ez][cmp][z+irp1]-f[Ez][cmp][z+ir])
                 - (f[Er][cmp][z+ir+1]-f[Er][cmp][z+ir]));
        }
      // Propogate Hz
      if (s->C[R][Hz])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
          double oorph = 1.0/(ir0 + r+0.5);
          double morph = m/((int)(v.origin_r()*v.a+0.5) + r+0.5);
          const int ir = r*(v.nz()+1);
          const int irp1 = (r+1)*(v.nz()+1);
          for (int z=1;z<=v.nz();z++) {
            const double Crhz = s->C[R][Hz][z+ir];
            const double ooop_Crhz = s->Cdecay[R][Hz][Z][z+ir];
            const double dhzr =
              ooop_Crhz*(-Courant*(f[Ep][cmp][z+irp1]*(ir0 + r+1.)-
                             f[Ep][cmp][z+ir]*(ir0 + r))*oorph
                         - Crhz*f_p_pml[Hz][cmp][z+ir]);
            f_p_pml[Hz][cmp][z+ir] += dhzr;
            f[Hz][cmp][z+ir] += dhzr + Courant*(it(cmp,f[Er],z+ir)*morph);
          }
        }
      else
        for (int r=rstart_0(v,m);r<v.nr();r++) {
          double oorph = 1.0/(ir0 + r+0.5);
          double morph = m/((int)(v.origin_r()*v.a+0.5) + r+0.5);
          const int ir = r*(v.nz()+1);
          const int irp1 = (r+1)*(v.nz()+1);
          for (int z=1;z<=v.nz();z++)
            f[Hz][cmp][z+ir] += Courant*
              (it(cmp,f[Er],z+ir)*morph
               - (f[Ep][cmp][z+irp1]*(ir0 + r+1.)-
                  f[Ep][cmp][z+ir]*(ir0 + r))*oorph);
        }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0) {
        // Nothing needed for H.
      } else if (m == 1 && v.origin_r() == 0.0) {
        if (s->C[Z][Hr])
          for (int z=0;z<v.nz();z++) {
            const double Czhr = s->C[Z][Hr][z];
            const double ooop_Czhr = s->Cdecay[Z][Hr][R][z];
            const double dhrp = Courant*(-it(cmp,f[Ez],z+(v.nz()+1))/* /1.0 */);
            const double hrz = f[Hr][cmp][z] - f_p_pml[Hr][cmp][z];
            f_p_pml[Hr][cmp][z] += dhrp;
            f[Hr][cmp][z] += dhrp +
              ooop_Czhr*(Courant*(f[Ep][cmp][z+1]-f[Ep][cmp][z]) - Czhr*hrz);
          }
        else
          for (int z=0;z<v.nz();z++)
            f[Hr][cmp][z] += Courant*
              ((f[Ep][cmp][z+1]-f[Ep][cmp][z]) - it(cmp,f[Ez],z+(v.nz()+1))/* /1.0 */);
      } else {
        for (int r=0;r<=v.nr() && (int)(v.origin_r()*v.a+0.5) + r < m;r++) {
          const int ir = r*(v.nz()+1);
          for (int z=0;z<=v.nz();z++) f[Hr][cmp][z+ir] = 0;
          if (f_p_pml[Hr][cmp])
            for (int z=0;z<=v.nz();z++) f_p_pml[Hr][cmp][z+ir] = 0;
        }
      }
    }
  } else {
    abort("Can't step H in these dimensions.");
  }
}

} // namespace meep

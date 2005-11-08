/* Copyright (C) 2005 Massachusetts Institute of Technology
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
#include "ran.hpp"

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

void fields::step_d() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_d();
}

void fields_chunk::step_d() {
  const volume v = this->v;
  if (v.dim != Dcyl) {
    const int n1 = num_each_direction[0];
    const int n2 = num_each_direction[1];
    const int n3 = num_each_direction[2];
    const int s1 = stride_each_direction[0];
    const int s2 = stride_each_direction[1];
    const int s3 = stride_each_direction[2];
    DOCMP FOR_D_COMPONENTS(cc)
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
#include "step_d.hpp"
      }
  } else if (v.dim == Dcyl) {
    int ir0 = int((v.origin_r() + rshift) * v.a + 0.5);
    DOCMP {
      // Propogate Dp
      if (f[Dp][cmp])
        if (s->C[Z][Dp] || s->C[R][Dp])
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czep = (s->C[Z][Dp])?s->C[Z][Dp][z+ir]:0;
              const double Crep = (s->C[R][Dp])?s->C[R][Dp][z+ir]:0;
              const double ooop_Czep = (s->Cdecay[Z][Dp][P]) ?
                s->Cdecay[Z][Dp][P][z+ir] : 1.0;
              const double ooop_Crep = (s->Cdecay[R][Dp][P]) ?
                s->Cdecay[R][Dp][P][z+ir] : 1.0;
              const double depz = ooop_Czep*(Courant*(f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                                             - Czep*f_p_pml[Dp][cmp][z+ir]);
              const double epr = f[Dp][cmp][z+ir] - f_p_pml[Dp][cmp][z+ir];
              f_p_pml[Dp][cmp][z+ir] += depz;
              f[Dp][cmp][z+ir] += depz +
                ooop_Crep*(Courant*(-(f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1])) - Crep*epr);
            }
          }
        else
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Dp][cmp][z+ir] += Courant*((f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                                     - (f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1]));
          }
      // Propogate Dz
      if (f[Dz][cmp])
        if (s->C[R][Dz])
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/(ir0 + r);
	    double mor = m/((int)(v.origin_r()*v.a + 0.5) + r);
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Crez = s->C[R][Dz][z+ir];
              const double ooop_Crez = (s->Cdecay[R][Dz][Z]) ?
                s->Cdecay[R][Dz][Z][z+ir] : 1.0;
              const double dezr = ooop_Crez*
                (Courant*(f[Hp][cmp][z+ir]*(ir0 + r+0.5)-
                    f[Hp][cmp][z+irm1]*(ir0 + r-0.5))*oor
                 - Crez*f_p_pml[Dz][cmp][z+ir]);
              f_p_pml[Dz][cmp][z+ir] += dezr;
              f[Dz][cmp][z+ir] += dezr + Courant*(-it(cmp,f[Hr],z+ir)*mor);
            }
          }
        else
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/(ir0 + r);
	    double mor = m/((int)(v.origin_r()*v.a + 0.5) + r);
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Dz][cmp][z+ir] += Courant*
                ((f[Hp][cmp][z+ir]*(ir0 + r+0.5)-
                  f[Hp][cmp][z+irm1]*(ir0 + r-0.5))*oor
                 - it(cmp,f[Hr],z+ir)*mor);
          }
      // Propogate Dr
      if (f[Dr][cmp])
        if (s->C[Z][Dr])
          for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin_r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czer = s->C[Z][Dr][z+ir];
              const double ooop_Czer = (s->Cdecay[Z][Dr][R]) ?
                s->Cdecay[Z][Dr][R][z+ir] : 1.0;
              double derp = Courant*(it(cmp,f[Hz],z+ir)*morph);
              double erz = f[Dr][cmp][z+ir] - f_p_pml[Dr][cmp][z+ir];
              f_p_pml[Dr][cmp][z+ir] += derp;
              f[Dr][cmp][z+ir] += derp + ooop_Czer*
                (-Courant*(f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]) - Czer*erz);
            }
          }
        else
          for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin_r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Dr][cmp][z+ir] += Courant*
                (it(cmp,f[Hz],z+ir)*morph - (f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]));
          }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0 && v.origin_r() == 0.0 && f[Dz][cmp]) {
        for (int z=0;z<=v.nz();z++)
          f[Dz][cmp][z] += Courant*(f[Hp][cmp][z] + it(cmp,f[Hr],z)*m);
      } else if (m == 1 && v.origin_r() == 0.0 && f[Dp][cmp]) {
        if (s->C[Z][Dp])
          for (int z=1;z<=v.nz();z++) {
            const double Czep = s->C[Z][Dp][z];
            const double ooop_Czep = (s->Cdecay[Z][Dp][P]) ?
              s->Cdecay[Z][Dp][P][z] : 1.0;
            const double depz = ooop_Czep*(Courant*(f[Hr][cmp][z]-f[Hr][cmp][z-1])
                                                  - Czep*f_p_pml[Dp][cmp][z]);
            f_p_pml[Dp][cmp][z] += depz;
            f[Dp][cmp][z] += depz + Courant*(-f[Hz][cmp][z]*2.0);
          }
        else
          for (int z=1;z<=v.nz();z++)
            f[Dp][cmp][z] += Courant*((f[Hr][cmp][z]-f[Hr][cmp][z-1]) - f[Hz][cmp][z]*2.0);
      } else {
        for (int r=0;r<=v.nr() && (int)(v.origin_r()*v.a+0.5) + r < m;r++) {
          const int ir = r*(v.nz()+1);
          for (int z=0;z<=v.nz();z++) f[Dp][cmp][z+ir] = 0;
          if (f_p_pml[Dp][cmp])
            for (int z=0;z<=v.nz();z++) f_p_pml[Dp][cmp][z+ir] = 0;
          for (int z=0;z<=v.nz();z++) f[Dz][cmp][z+ir] = 0;
          if (f_p_pml[Dz][cmp])
            for (int z=0;z<=v.nz();z++) f_p_pml[Dz][cmp][z+ir] = 0;
        }
      }
    }
  } else {
    abort("Unsupported dimension.\n");
  }
}

} // namespace meep

/* Copyright (C) 2003 Massachusetts Institute of Technology
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

#include "meep.h"
#include "meep_internals.h"

#define RESTRICT

void fields::step() {
  am_now_working_on(Stepping);
  phase_material();

  //for (int i=0;i<num_chunks;i++)
  //  master_printf("Field is now %lg\n", chunks[i]->peek_field(Ex,vec2d(1.55,0.6)));
  step_h();
  step_h_source();
  step_boundaries(H_stuff);
  // because step_boundaries overruns the timing stack...
  am_now_working_on(Stepping);

  step_d();
  step_e_source();
  step_boundaries(D_stuff);

  prepare_step_polarization_energy();
  half_step_polarization_energy();
  update_e_from_d();
  //step_boundaries(E_stuff);  FIXME: Need to enable this when I encode anisotropy, etc.

  // because step_boundaries overruns the timing stack...
  am_now_working_on(Stepping);

  half_step_polarization_energy();

  update_polarization_saturation();
  step_polarization_itself();

  update_fluxes();
  t += 1;
  finished_working();
}

double fields_chunk::peek_field(component c, const vec &where) {
  double w[8];
  ivec ilocs[8];
  v.interpolate(c,where, ilocs, w);
  if (v.contains(ilocs[0]) && f[c][0]) {
    double hello = 0.0;
    if (is_mine()) hello = f[c][0][v.index(c,ilocs[0])];
    broadcast(n_proc(), &hello, 1);
    return hello;
  }
  //abort("Got no such %s field at %lg %lg!\n",
  //      component_name(c), v[ilocs[0]].x(), v[ilocs[0]].y());
  return 0.0;
}

void fields::phase_material() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->phase_material(phasein_time);
  phasein_time--;
}

void fields_chunk::phase_material(int phasein_time) {
  if (new_ma && phasein_time > 0) {
    // First multiply the electric field by epsilon...
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(c)
        if (f[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f[c][cmp][i] /= ma->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_p_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_p_pml[c][cmp][i] /= ma->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_m_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_m_pml[c][cmp][i] /= ma->inveps[c][component_direction(c)][i];
    }
    // Then change epsilon...
    ma->mix_with(new_ma, 1.0/phasein_time);
    // Finally divide the electric field by the new epsilon...
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(c)
        if (f[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f[c][cmp][i] *= ma->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_p_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_p_pml[c][cmp][i] *= ma->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_m_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_m_pml[c][cmp][i] *= ma->inveps[c][component_direction(c)][i];
    }
  }
}

inline double it(int cmp, double *(f[2]), int ind) { return (1-2*cmp)*f[1-cmp][ind]; }

inline int rstart_0(const volume &v, int m) {
  return (int) max(0.0, m - (int)(v.origin.r()*v.a+0.5) - 1.0);
}
inline int rstart_1(const volume &v, int m) {
  return (int) max(1.0, (double)m - (int)(v.origin.r()*v.a+0.5));
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
        const bool have_p_pml = have_p && ma->C[d_deriv_p][cc];
        const bool have_m_pml = have_m && ma->C[d_deriv_m][cc];
        const int stride_p = (have_p)?v.stride(d_deriv_p):0;
        const int stride_m = (have_m)?v.stride(d_deriv_m):0;
        // The following lines "promise" the compiler that the values of
        // these arrays won't change during this loop.
        RESTRICT const double *C_m = (have_m_pml)?ma->C[d_deriv_m][cc] + yee_idx:NULL;
        RESTRICT const double *C_p = (have_p_pml)?ma->C[d_deriv_p][cc] + yee_idx:NULL;
        RESTRICT const double *decay_m = (!have_m_pml)?NULL:
          ma->Cdecay[d_deriv_m][cc][component_direction(cc)] + yee_idx;
        RESTRICT const double *decay_p = (!have_p_pml)?NULL:
          ma->Cdecay[d_deriv_p][cc][component_direction(cc)] + yee_idx;
        RESTRICT const double *f_p = (have_p)?f[c_p][cmp] + v.yee_index(c_p):NULL;
        RESTRICT const double *f_m = (have_m)?f[c_m][cmp] + v.yee_index(c_m):NULL;
        RESTRICT double *the_f = f[cc][cmp] + yee_idx;
        RESTRICT double *the_f_p_pml = f_p_pml[cc][cmp] + yee_idx;
        RESTRICT double *the_f_m_pml = f_m_pml[cc][cmp] + yee_idx;
#include "step_h.h"
      }
  } else if (v.dim == Dcyl) {
    DOCMP {
      // Propogate Hr
      if (ma->C[Z][Hr])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/((int)(v.origin.r()*v.a+0.5) + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Czhr = ma->C[Z][Hr][z+ir];
              const double ooop_Czhr = ma->Cdecay[Z][Hr][R][z+ir];
              double dhrp = c*(-it(cmp,f[Ez],z+ir)*mor);
              double hrz = f[Hr][cmp][z+ir] - f_p_pml[Hr][cmp][z+ir];
              f_p_pml[Hr][cmp][z+ir] += dhrp;
              f[Hr][cmp][z+ir] += dhrp +
                ooop_Czhr*(c*(f[Ep][cmp][z+ir+1]-f[Ep][cmp][z+ir]) - Czhr*hrz);
            }
          }
      else
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Hr][cmp][z+ir] += c*
                ((f[Ep][cmp][z+ir+1]-f[Ep][cmp][z+ir]) - it(cmp,f[Ez],z+ir)*mor);
          }
      // Propogate Hp
      if (ma->C[Z][Hp] || ma->C[R][Hp])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Czhp = (ma->C[Z][Hp])?ma->C[Z][Hp][z+ir]:0;
              const double Crhp = (ma->C[R][Hp])?ma->C[R][Hp][z+ir]:0;
              const double ooop_Czhp = (ma->Cdecay[Z][Hp][P]) ?
                ma->Cdecay[Z][Hp][P][z+ir]:1.0;
              const double ooop_Crhp = (ma->Cdecay[R][Hp][P]) ?
                ma->Cdecay[R][Hp][P][z+ir]:1.0;
              const double dhpz = ooop_Czhp*(-c*(f[Er][cmp][z+ir+1]-f[Er][cmp][z+ir])
                                             - Czhp*f_p_pml[Hp][cmp][z+ir]);
              const double hpr = f[Hp][cmp][z+ir]-f_p_pml[Hp][cmp][z+ir];
              f_p_pml[Hp][cmp][z+ir] += dhpz;
              f[Hp][cmp][z+ir] += dhpz +
                ooop_Czhp*(c*(f[Ez][cmp][z+irp1]-f[Ez][cmp][z+ir]) - Crhp*hpr);
            }
          }
      else 
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Hp][cmp][z+ir] += c*
                ((f[Ez][cmp][z+irp1]-f[Ez][cmp][z+ir])
                 - (f[Er][cmp][z+ir+1]-f[Er][cmp][z+ir]));
        }
      // Propogate Hz
      if (ma->C[R][Hz])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
          double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
          double morph = m*oorph;
          const int ir = r*(v.nz()+1);
          const int irp1 = (r+1)*(v.nz()+1);
          for (int z=1;z<=v.nz();z++) {
            const double Crhz = ma->C[R][Hz][z+ir];
            const double ooop_Crhz = ma->Cdecay[R][Hz][Z][z+ir];
            const double dhzr =
              ooop_Crhz*(-c*(f[Ep][cmp][z+irp1]*((int)(v.origin.r()*v.a+0.5) + r+1.)-
                             f[Ep][cmp][z+ir]*((int)(v.origin.r()*v.a+0.5) + r))*oorph
                         - Crhz*f_p_pml[Hz][cmp][z+ir]);
            const double hzp = f[Hz][cmp][z+ir] - f_p_pml[Hz][cmp][z+ir];
            f_p_pml[Hz][cmp][z+ir] += dhzr;
            f[Hz][cmp][z+ir] += dhzr + c*(it(cmp,f[Er],z+ir)*morph);
          }
        }
      else
        for (int r=rstart_0(v,m);r<v.nr();r++) {
          double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
          double morph = m*oorph;
          const int ir = r*(v.nz()+1);
          const int irp1 = (r+1)*(v.nz()+1);
          for (int z=1;z<=v.nz();z++)
            f[Hz][cmp][z+ir] += c*
              (it(cmp,f[Er],z+ir)*morph
               - (f[Ep][cmp][z+irp1]*((int)(v.origin.r()*v.a+0.5) + r+1.)-
                  f[Ep][cmp][z+ir]*((int)(v.origin.r()*v.a+0.5) + r))*oorph);
        }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0) {
        // Nothing needed for H.
      } else if (m == 1 && v.origin.r() == 0.0) {
        if (ma->C[Z][Hr])
          for (int z=0;z<v.nz();z++) {
            const double Czhr = ma->C[Z][Hr][z];
            const double ooop_Czhr = ma->Cdecay[Z][Hr][R][z];
            const double dhrp = c*(-it(cmp,f[Ez],z+(v.nz()+1))/* /1.0 */);
            const double hrz = f[Hr][cmp][z] - f_p_pml[Hr][cmp][z];
            f_p_pml[Hr][cmp][z] += dhrp;
            f[Hr][cmp][z] += dhrp +
              ooop_Czhr*(c*(f[Ep][cmp][z+1]-f[Ep][cmp][z]) - Czhr*hrz);
          }
        else
          for (int z=0;z<v.nz();z++)
            f[Hr][cmp][z] += c*
              ((f[Ep][cmp][z+1]-f[Ep][cmp][z]) - it(cmp,f[Ez],z+(v.nz()+1))/* /1.0 */);
      } else {
        for (int r=0;r<=v.nr() && (int)(v.origin.r()*v.a+0.5) + r < m;r++) {
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
        const bool have_p_pml = have_p && ma->C[d_deriv_p][cc];
        const bool have_m_pml = have_m && ma->C[d_deriv_m][cc];
        const int stride_p = (have_p)?v.stride(d_deriv_p):0;
        const int stride_m = (have_m)?v.stride(d_deriv_m):0;
        // The following lines "promise" the compiler that the values of
        // these arrays won't change during this loop.
        RESTRICT const double *C_m = (have_m_pml)?ma->C[d_deriv_m][cc] + yee_idx:NULL;
        RESTRICT const double *C_p = (have_p_pml)?ma->C[d_deriv_p][cc] + yee_idx:NULL;
        RESTRICT const double *decay_m = (!have_m_pml)?NULL:
          ma->Cdecay[d_deriv_m][cc][component_direction(cc)] + yee_idx;
        RESTRICT const double *decay_p = (!have_p_pml)?NULL:
          ma->Cdecay[d_deriv_p][cc][component_direction(cc)] + yee_idx;
        RESTRICT const double *f_p = (have_p)?f[c_p][cmp] + v.yee_index(c_p):NULL;
        RESTRICT const double *f_m = (have_m)?f[c_m][cmp] + v.yee_index(c_m):NULL;
        RESTRICT double *the_f = f[cc][cmp] + yee_idx;
        RESTRICT double *the_f_p_pml = f_p_pml[cc][cmp] + yee_idx;
        RESTRICT double *the_f_m_pml = f_m_pml[cc][cmp] + yee_idx;
#include "step_d.h"
      }
  } else if (v.dim == Dcyl) {
    DOCMP {
      // Propogate Dp
      if (f[Dp][cmp])
        if (ma->C[Z][Dp] || ma->C[R][Dp])
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
            const double mor = m*oor;
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czep = (ma->C[Z][Dp])?ma->C[Z][Dp][z+ir]:0;
              const double Crep = (ma->C[R][Dp])?ma->C[R][Dp][z+ir]:0;
              const double ooop_Czep = (ma->Cdecay[Z][Dp][P]) ?
                ma->Cdecay[Z][Dp][P][z+ir] : 1.0;
              const double ooop_Crep = (ma->Cdecay[R][Dp][P]) ?
                ma->Cdecay[R][Dp][P][z+ir] : 1.0;
              const double depz = ooop_Czep*(c*(f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                                             - Czep*f_p_pml[Dp][cmp][z+ir]);
              const double epr = f[Dp][cmp][z+ir] - f_p_pml[Dp][cmp][z+ir];
              f_p_pml[Dp][cmp][z+ir] += depz;
              f[Dp][cmp][z+ir] += depz +
                ooop_Crep*(c*(-(f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1])) - Crep*epr);
            }
          }
        else
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Dp][cmp][z+ir] += c*((f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                                     - (f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1]));
          }
      // Propogate Dz
      if (f[Dz][cmp])
        if (ma->C[R][Dz])
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Crez = ma->C[R][Dz][z+ir];
              const double ooop_Crez = (ma->Cdecay[R][Dz][Z]) ?
                ma->Cdecay[R][Dz][Z][z+ir] : 1.0;
              const double dezr = ooop_Crez*
                (c*(f[Hp][cmp][z+ir]*((int)(v.origin.r()*v.a+0.5) + r+0.5)-
                    f[Hp][cmp][z+irm1]*((int)(v.origin.r()*v.a+0.5) + r-0.5))*oor
                 - Crez*f_p_pml[Dz][cmp][z+ir]);
              const double ezp = f[Dz][cmp][z+ir]-f_p_pml[Dz][cmp][z+ir];
              f_p_pml[Dz][cmp][z+ir] += dezr;
              f[Dz][cmp][z+ir] += dezr + c*(-it(cmp,f[Hr],z+ir)*mor);
            }
          }
        else
          for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Dz][cmp][z+ir] += c*
                ((f[Hp][cmp][z+ir]*((int)(v.origin.r()*v.a+0.5) + r+0.5)-
                  f[Hp][cmp][z+irm1]*((int)(v.origin.r()*v.a+0.5) + r-0.5))*oor
                 - it(cmp,f[Hr],z+ir)*mor);
          }
      // Propogate Dr
      if (f[Dr][cmp])
        if (ma->C[Z][Dr])
          for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czer = ma->C[Z][Dr][z+ir];
              const double ooop_Czer = (ma->Cdecay[Z][Dr][R]) ?
                ma->Cdecay[Z][Dr][R][z+ir] : 1.0;
              double derp = c*(it(cmp,f[Hz],z+ir)*morph);
              double erz = f[Dr][cmp][z+ir] - f_p_pml[Dr][cmp][z+ir];
              f_p_pml[Dr][cmp][z+ir] += derp;
              f[Dr][cmp][z+ir] += derp + ooop_Czer*
                (-c*(f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]) - Czer*erz);
            }
          }
        else
          for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Dr][cmp][z+ir] += c*
                (it(cmp,f[Hz],z+ir)*morph - (f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]));
          }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0 && v.origin.r() == 0.0 && f[Dz][cmp]) {
        for (int z=0;z<=v.nz();z++)
          f[Dz][cmp][z] += c*(f[Hp][cmp][z] + it(cmp,f[Hr],z)*m);
      } else if (m == 1 && v.origin.r() == 0.0 && f[Dp][cmp]) {
        if (ma->C[Z][Dp])
          for (int z=1;z<=v.nz();z++) {
            const double Czep = ma->C[Z][Dp][z];
            const double ooop_Czep = (ma->Cdecay[Z][Dp][P]) ?
              ma->Cdecay[Z][Dp][P][z] : 1.0;
            const double depz = ooop_Czep*(c*(f[Hr][cmp][z]-f[Hr][cmp][z-1])
                                                  - Czep*f_p_pml[Dp][cmp][z]);
            const double epr = f[Dp][cmp][z] - f_p_pml[Dp][cmp][z];
            f_p_pml[Dp][cmp][z] += depz;
            f[Dp][cmp][z] += depz + c*(-f[Hz][cmp][z]*2.0);
          }
        else
          for (int z=1;z<=v.nz();z++)
            f[Dp][cmp][z] += c*((f[Hr][cmp][z]-f[Hr][cmp][z-1]) - f[Hz][cmp][z]*2.0);
      } else {
        for (int r=0;r<=v.nr() && (int)(v.origin.r()*v.a+0.5) + r < m;r++) {
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

void fields::update_e_from_d() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->update_e_from_d();
}

#define FOR_POLARIZATIONS(init, p) for (polarization *p = init; p; p = p->next)

void fields_chunk::update_e_from_d() {
  const int ntot = ma->v.ntot();
#include "step_e.h"
}

#include "config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

void fields::step_boundaries(field_type ft) {
  am_now_working_on(MpiTime);
  // Do the metals first!
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) chunks[i]->zero_metal(ft);
  // First copy outgoing data to buffers...
  int *wh = new int[num_chunks];
  for (int i=0;i<num_chunks;i++) wh[i] = 0;
  for (int i=0;i<num_chunks;i++)
    for (int j=0;j<num_chunks;j++)
      if (chunks[j]->is_mine()) {
        const int pair = j+i*num_chunks;
        for (int n=0;n<comm_sizes[ft][pair];n++)
          comm_blocks[ft][pair][n] =
            *(chunks[j]->connections[ft][Outgoing][wh[j]++]);
      }
  // Communicate the data around!
#if 0 // This is the blocking version, which should always be safe!
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      const int pair = j+i*num_chunks;
      DOCMP {
        send(chunks[j]->n_proc(), chunks[i]->n_proc(),
             comm_blocks[ft][pair], comm_sizes[ft][pair]);
      }
    }
#endif
#ifdef HAVE_MPI
  const int maxreq = num_chunks*num_chunks;
  MPI_Request *reqs = new MPI_Request[maxreq];
  MPI_Status *stats = new MPI_Status[maxreq];
  int reqnum = 0;
  int *tagto = new int[count_processors()];
  for (int i=0;i<count_processors();i++) tagto[i] = 0;
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      if (chunks[j]->n_proc() != chunks[i]->n_proc()) {
        const int pair = j+i*num_chunks;
        if (comm_sizes[ft][pair] > 0) {
          if (chunks[j]->is_mine())
            MPI_Isend(comm_blocks[ft][pair], comm_sizes[ft][pair],
                      MPI_DOUBLE, chunks[i]->n_proc(),
                      tagto[chunks[i]->n_proc()]++,
                      MPI_COMM_WORLD, &reqs[reqnum++]);
          if (chunks[i]->is_mine())
            MPI_Irecv(comm_blocks[ft][pair], comm_sizes[ft][pair],
                      MPI_DOUBLE, chunks[j]->n_proc(),
                      tagto[chunks[j]->n_proc()]++,
                      MPI_COMM_WORLD, &reqs[reqnum++]);
        }
      }
    }
  delete[] tagto;
  if (reqnum > maxreq) abort("Too many requests!!!\n");
  if (reqnum > 0) MPI_Waitall(reqnum, reqs, stats);
  delete[] reqs;
  delete[] stats;
#endif
  
  // Finally, copy incoming data to the fields themselves!
  for (int i=0;i<num_chunks;i++) {
    int wh = 0;
    if (chunks[i]->is_mine())
      for (int j=0;j<num_chunks;j++) {
        const int pair = j+i*num_chunks;
        int n;
        for (n=0;n<comm_num_complex[ft][pair];n+=2) {
          const double phr = real(chunks[i]->connection_phases[ft][wh/2]);
          const double phi = imag(chunks[i]->connection_phases[ft][wh/2]);
          *(chunks[i]->connections[ft][Incoming][wh]) =
            phr*comm_blocks[ft][pair][n] - phi*comm_blocks[ft][pair][n+1];
          *(chunks[i]->connections[ft][Incoming][wh+1]) =
            phr*comm_blocks[ft][pair][n+1] + phi*comm_blocks[ft][pair][n];
          wh += 2;
        }
        for (;n<comm_num_complex[ft][pair]+comm_num_negate[ft][pair];n++) {
          *(chunks[i]->connections[ft][Incoming][wh]) = -comm_blocks[ft][pair][n];
          wh++;
        }
        for (;n<comm_sizes[ft][pair];n++) {
          *(chunks[i]->connections[ft][Incoming][wh]) = comm_blocks[ft][pair][n];
          wh++;
        }
      }
  }
  delete[] wh;
  finished_working();
}

void fields::step_h_source() {
  const double tim = time();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_h_source(chunks[i]->h_sources, tim);
}

void fields_chunk::step_h_source(const src *s, double time) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(time);
  if (A == 0.0) {
    step_h_source(s->next, time);
    return;
  }
  FOR_COMPONENTS(c)
    if (f[c][0] && is_magnetic(c)) {
      f[c][0][s->i] += real(A*s->A[c]);
      if (!is_real) f[c][1][s->i] += imag(A*s->A[c]);
    }
  step_h_source(s->next, time);
}

void fields::step_e_source() {
  const double tim = time()+0.5*inva*c;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_e_source(chunks[i]->e_sources, tim);
}

void fields_chunk::step_e_source(const src *s, double time) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(time);
  if (A == 0.0) {
    step_e_source(s->next, time);
    return;
  }
  FOR_E_AND_D(c,dc)
    if (f[c][0]) {
      f[dc][0][s->i] += real(A*s->A[c]);
      if (!is_real) f[dc][1][s->i] += imag(A*s->A[c]);
    }
  step_e_source(s->next, time);
}

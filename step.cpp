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

#include "dactyl.h"
#include "dactyl_internals.h"

void fields::step_old() {
  am_now_working_on(Stepping);
  phase_material();

  //for (int i=0;i<num_chunks;i++)
  //  master_printf("Field is now %lg\n", chunks[i]->peek_field(Ex,vec2d(1.55,0.6)));
  step_h_old();
  step_h_source();
  step_boundaries(H_stuff);
  // because step_boundaries overruns the timing stack...
  am_now_working_on(Stepping);

  prepare_step_polarization_energy();
  half_step_polarization_energy();
  step_e();
  step_e_source();
  step_e_polarization();
  step_boundaries(E_stuff);
  // because step_boundaries overruns the timing stack...
  am_now_working_on(Stepping);

  half_step_polarization_energy();

  update_polarization_saturation();
  step_polarization_itself();

  update_fluxes();
  t += 1;
  finished_working();
}

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

  prepare_step_polarization_energy();
  half_step_polarization_energy();
  step_e();
  step_e_source();
  step_e_polarization();
  step_boundaries(E_stuff);
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
  int inds[8];
  double w[8];
  v.interpolate(c,where,inds,w);
  const int ind = inds[0];
  if (ind >= 0 && ind < v.ntot()) {
    double hello = 0.0;
    if (is_mine()) {
      hello = f[c][0][ind];
    }
    //master_printf("Found the field here at %lg %lg!\n",
    //              v.loc(c,ind).x(), v.loc(c,ind).y());
    broadcast(n_proc(), &hello, 1);
    return hello;
  }
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
        if (f[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_pml[c][cmp][i] /= ma->inveps[c][component_direction(c)][i];
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
        if (f[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_pml[c][cmp][i] *= ma->inveps[c][component_direction(c)][i];
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

void fields::step_h_old() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_h_old();
}

void fields_chunk::step_h_old() {
  DOCMP {
    // Propogate Hz
    if ((ma->C[X][Hz] || ma->C[Y][Hz])&& f[Hz][cmp])
      for (int x=0;x<v.nx();x++) {
        const int ix = x*(v.ny()+1);
        const int ixp1 = (x+1)*(v.ny()+1);
        for (int y=0;y<v.ny();y++) {
          const double Cxhz = (ma->C[X][Hz])?ma->C[X][Hz][y+ix]:0;
          const double Cyhz = (ma->C[Y][Hz])?ma->C[Y][Hz][y+ix]:0;
          const double ooop_Cxhz = (ma->Cdecay[X][Hz][Z]) ?
            ma->Cdecay[X][Hz][Z][y+ix]:1.0;
          const double ooop_Cyhz = (ma->Cdecay[Y][Hz][Z]) ?
            ma->Cdecay[Y][Hz][Z][y+ix] : 1.0;
          const double dhzx = ooop_Cxhz*
            (c*(f[Ey][cmp][y+ixp1] - f[Ey][cmp][y+ix]) - Cxhz*f_pml[Hz][cmp][y+ix]);
          const double hzy = f[Hz][cmp][y+ix] - f_pml[Hz][cmp][y+ix];
          f_pml[Hz][cmp][y+ix] += dhzx;
          f[Hz][cmp][y+ix] += dhzx +
            ooop_Cyhz*(c*(-(f[Ex][cmp][y+ix+1] - f[Ex][cmp][y+ix])) - Cyhz*hzy);
        }
      }
    else if (f[Hz][cmp])
      for (int x=0;x<v.nx();x++) {
        const int ix = x*(v.ny()+1);
        const int ixp1 = (x+1)*(v.ny()+1);
        for (int y=0;y<v.ny();y++)
          f[Hz][cmp][y+ix] += c*((f[Ey][cmp][y+ixp1] - f[Ey][cmp][y+ix]) -
                                 (f[Ex][cmp][y+ix+1] - f[Ex][cmp][y+ix]));
      }
    // Propogate Hx
    if (ma->C[Y][Hx] && f[Hx][cmp])
      for (int x=1;x<=v.nx();x++) {
        const int ix = x*(v.ny()+1);
        for (int y=0;y<v.ny();y++) {
          const double Cyhx = ma->C[Y][Hx][y+ix];
          const double ooop_Cyhx = ma->Cdecay[Y][Hx][X][y+ix];
          f[Hx][cmp][y+ix] += ooop_Cyhx*
            (c*(f[Ez][cmp][y+ix+1] - f[Ez][cmp][y+ix]) - Cyhx*f[Hx][cmp][y+ix]);
        }
      }
    else if (f[Hx][cmp])
      for (int x=1;x<=v.nx();x++) {
        const int ix = x*(v.ny()+1);
        for (int y=0;y<v.ny();y++)
          f[Hx][cmp][y+ix] += c*(f[Ez][cmp][y+ix+1] - f[Ez][cmp][y+ix]);
      }
    // Propogate Hy
    if (ma->C[X][Hy] && f[Hy][cmp])
      for (int x=0;x<v.nx();x++) {
        const int ix = x*(v.ny()+1);
        const int ixp1 = (x+1)*(v.ny()+1);
        for (int y=1;y<=v.ny();y++) {
          const double Cxhy = ma->C[X][Hy][y+ix];
          const double ooop_Cxhy = ma->Cdecay[X][Hy][Y][y+ix];
          f[Hy][cmp][y+ix] += ooop_Cxhy*
            (c*(-(f[Ez][cmp][y+ixp1] - f[Ez][cmp][y+ix])) - Cxhy*f[Hy][cmp][y+ix]);
        }
      }
    else if (f[Hy][cmp])
      for (int x=0;x<v.nx();x++) {
        const int ix = x*(v.ny()+1);
        const int ixp1 = (x+1)*(v.ny()+1);
        for (int y=1;y<=v.ny();y++)
          f[Hy][cmp][y+ix] += c*(-(f[Ez][cmp][y+ixp1] - f[Ez][cmp][y+ix]));
      }
  }
}

static inline bool cross_negative(direction a, direction b) {
  return ((3+b-a)%3) == 2;
}

static inline direction cross(direction a, direction b) {
  if (a < R && b < R) return (direction)((3+2*a-b)%3);
  return (direction) (2 + (3+2*(a-2)-(b-2))%3);
}

void fields_chunk::step_h() {
  const volume v = this->v;
  if (v.dim == D1) {
    DOCMP {
      if (ma->C[Z][Hy])
        for (int z=0;z<v.nz();z++) {
          const double C = ma->C[Z][Hy][z];
          const double ooop_C = ma->Cdecay[Z][Hy][Y][z];
          f[Hy][cmp][z] += ooop_C*(-c*(f[Ex][cmp][z+1] - f[Ex][cmp][z])
                                   - C*f[Hy][cmp][z]);
        }
      else
        for (int z=0;z<v.nz();z++)
          f[Hy][cmp][z] += -c*(f[Ex][cmp][z+1] - f[Ex][cmp][z]);
    }
  } else if (v.dim == D2) {
    DOCMP FOR_MAGNETIC_COMPONENTS(cc)
      if (f[cc][cmp]) {
        const int yee_idx = v.yee_index(cc);
        const component c_p=plus_component[cc], c_m=minus_component[cc];
        const direction d_deriv_p = plus_deriv_direction[cc];
        const direction d_deriv_m = minus_deriv_direction[cc];
        const bool have_p = have_plus_deriv[cc];
        const bool have_m = have_minus_deriv[cc];
        const bool have_p_pml = have_p && have_pml_in_direction[d_deriv_p];
        const bool have_m_pml = have_m && have_pml_in_direction[d_deriv_m];
        const int stride_p = (have_p)?v.stride(d_deriv_p):0;
        const int stride_m = (have_m)?v.stride(d_deriv_m):0;
        // The following lines "promise" the compiler that the values of
        // these arrays won't change during this loop.
        const double *C_m = (have_m_pml)?ma->C[d_deriv_m][cc] + yee_idx:NULL;
        const double *C_p = (have_p_pml)?ma->C[d_deriv_p][cc] + yee_idx:NULL;
        const double *decay_m = (!have_m_pml)?NULL:
          ma->Cdecay[d_deriv_m][cc][component_direction(cc)] + yee_idx;
        const double *decay_p = (!have_p_pml)?NULL:
          ma->Cdecay[d_deriv_p][cc][component_direction(cc)] + yee_idx;
        const double *f_p = (have_p)?f[c_p][cmp] + v.yee_index(c_p):NULL;
        const double *f_m = (have_m)?f[c_m][cmp] + v.yee_index(c_m):NULL;
        double *the_f = f[cc][cmp] + yee_idx;
        double *the_f_pml = f_pml[cc][cmp] + yee_idx;
        const int n1 = num_each_direction[0];
        const int n2 = num_each_direction[1];
        const int n3 = num_each_direction[2];
        const int s1 = stride_each_direction[0];
        const int s2 = stride_each_direction[1];
        const int s3 = stride_each_direction[2];
#include "step_h.h"
        /*for (int i1 = 0; i1 < n1; i1++)
        for (int i2 = 0; i2 < n2; i2++)
        for (int i3 = 0, ind = i1*s1+i2*s2; i3 < n3; i3++, ind += s3) {
        //LOOP_OVER_OWNED(v, ind) {
          double dh = 0.0;
          if (have_m) {
            const double deriv_m = c*(f_m[ind]-f_m[ind-stride_m]);
            if (have_m_pml) {
              const double h_minus =
                the_f[ind] - the_f_pml[ind];
              dh += decay_m[ind]*
                (deriv_m - C_m[ind]*h_minus);
            } else {
              dh += deriv_m;
            }
          }
          if (have_p) {
            const double m_deriv_p = c*(f_p[ind-stride_p]-f_p[ind]);
            if (have_p_pml) {
              const double delta_p = decay_p[ind]*
                (m_deriv_p - C_p[ind]*the_f_pml[ind]);
              the_f_pml[ind] += delta_p;
              dh += delta_p;
            } else {
              if (have_m_pml) the_f_pml[ind] += m_deriv_p;
              dh += m_deriv_p;
            }
          }
          the_f[ind] += dh;
        }*/
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
              double hrz = f[Hr][cmp][z+ir] - f_pml[Hr][cmp][z+ir];
              f_pml[Hr][cmp][z+ir] += dhrp;
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
                                             - Czhp*f_pml[Hp][cmp][z+ir]);
              const double hpr = f[Hp][cmp][z+ir]-f_pml[Hp][cmp][z+ir];
              f_pml[Hp][cmp][z+ir] += dhpz;
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
                         - Crhz*f_pml[Hz][cmp][z+ir]);
            const double hzp = f[Hz][cmp][z+ir] - f_pml[Hz][cmp][z+ir];
            f_pml[Hz][cmp][z+ir] += dhzr;
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
            const double hrz = f[Hr][cmp][z] - f_pml[Hr][cmp][z];
            f_pml[Hr][cmp][z] += dhrp;
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
          if (f_pml[Hr][cmp])
            for (int z=0;z<=v.nz();z++) f_pml[Hr][cmp][z+ir] = 0;
        }
      }
    }
  } else {
    abort("Can't step H in these dimensions.");
  }
}

void fields::step_e() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_e();
}

void fields_chunk::step_e() {
  const volume v = this->v;
  if (v.dim == D1) {
    DOCMP {
      if (ma->C[Z][Ex])
        for (int z=1;z<=v.nz();z++) {
          const double C = ma->C[Z][Ex][z];
          const double inveps_plus_C = ma->Cdecay[Z][Ex][X][z];
          f[Ex][cmp][z] += inveps_plus_C*(-c*(f[Hy][cmp][z] - f[Hy][cmp][z-1])
                                          - C*f[Ex][cmp][z]);
        }
      else 
        for (int z=1;z<=v.nz();z++)
          f[Ex][cmp][z] += -c*ma->inveps[Ex][X][z]*
            (f[Hy][cmp][z] - f[Hy][cmp][z-1]);
    }
  } else if (v.dim == D2) {
    DOCMP {
      // Propogate Ez
      if ((ma->C[X][Ez] || ma->C[Y][Ez])&& f[Ez][cmp])
        for (int x=1;x<=v.nx();x++) {
          const int ix = x*(v.ny()+1);
          const int ixm1 = (x-1)*(v.ny()+1);
          for (int y=1;y<=v.ny();y++) {
            const double Cxez = (ma->C[X][Ez])?ma->C[X][Ez][y+ix]:0;
            const double Cyez = (ma->C[Y][Ez])?ma->C[Y][Ez][y+ix]:0;
            const double inveps_plus_Cxez = (ma->Cdecay[X][Ez][Z]) ?
              ma->Cdecay[X][Ez][Z][y+ix] : ma->inveps[Ez][Z][y+ix];
            const double inveps_plus_Cyez = (ma->Cdecay[Y][Ez][Z]) ?
              ma->Cdecay[Y][Ez][Z][y+ix] : ma->inveps[Ez][Z][y+ix];
            const double dezx = inveps_plus_Cxez*
              (-c*(f[Hy][cmp][y+ix] - f[Hy][cmp][y+ixm1]) - Cxez*f_pml[Ez][cmp][y+ix]);
            const double ezy = f[Ez][cmp][y+ix] - f_pml[Ez][cmp][y+ix];
            f_pml[Ez][cmp][y+ix] += dezx;
            f[Ez][cmp][y+ix] += dezx +
              inveps_plus_Cyez*(-c*(-(f[Hx][cmp][y+ix] - f[Hx][cmp][y+ix-1])) - Cyez*ezy);
          }
        }
      else if (f[Ez][cmp])
        for (int x=1;x<=v.nx();x++) {
          const int ix = x*(v.ny()+1);
          const int ixm1 = (x-1)*(v.ny()+1);
          for (int y=1;y<=v.ny();y++) {
            f[Ez][cmp][y+ix] += -c*ma->inveps[Ez][Z][y+ix]*
              ((f[Hy][cmp][y+ix] - f[Hy][cmp][y+ixm1]) -
               (f[Hx][cmp][y+ix] - f[Hx][cmp][y+ix-1]));
          }
        }
      // Propogate Ex
      if (ma->C[Y][Ex] && f[Ex][cmp])
        for (int x=0;x<v.nx();x++) {
          const int ix = x*(v.ny()+1);
          for (int y=1;y<=v.ny();y++) {
            const double Cyex = ma->C[Y][Ex][y+ix];
            const double inveps_plus_Cyex = (ma->Cdecay[Y][Ex][X]) ?
              ma->Cdecay[Y][Ex][X][y+ix] : ma->inveps[Ex][X][y+ix];
            f[Ex][cmp][y+ix] += inveps_plus_Cyex*
              (-c*(f[Hz][cmp][y+ix] - f[Hz][cmp][y+ix-1]) - Cyex*f[Ex][cmp][y+ix]);
          }
        }
      else if (f[Ex][cmp])
        for (int x=0;x<v.nx();x++) {
          const int ix = x*(v.ny()+1);
          for (int y=1;y<=v.ny();y++)
            f[Ex][cmp][y+ix] += -c*ma->inveps[Ex][X][y+ix]*
              (f[Hz][cmp][y+ix] - f[Hz][cmp][y+ix-1]);
        }
      // Propogate Ey
      if (ma->C[X][Ey] && f[Ey][cmp])
        for (int x=1;x<=v.nx();x++) {
          const int ix = x*(v.ny()+1);
          const int ixm1 = (x-1)*(v.ny()+1);
          for (int y=0;y<v.ny();y++) {
            const double Cxey = ma->C[X][Ey][y+ix];
            const double inveps_plus_Cxey = (ma->Cdecay[X][Ey][Y]) ?
              ma->Cdecay[X][Ey][Y][y+ix] : ma->inveps[Ey][Y][y+ix];
            f[Ey][cmp][y+ix] += inveps_plus_Cxey*
              (-c*(-(f[Hz][cmp][y+ix] - f[Hz][cmp][y+ixm1])) - Cxey*f[Ey][cmp][y+ix]);
          }
        }
      else if (f[Ey][cmp])
        for (int x=1;x<=v.nx();x++) {
          const int ix = x*(v.ny()+1);
          const int ixm1 = (x-1)*(v.ny()+1);
          for (int y=0;y<v.ny();y++)
            f[Ey][cmp][y+ix] += -c*ma->inveps[Ey][Y][y+ix]*
              (-(f[Hz][cmp][y+ix] - f[Hz][cmp][y+ixm1]));
        }
    }
  } else if (v.dim == Dcyl) {
    DOCMP {
      // Propogate Ep
      if (ma->C[Z][Ep] || ma->C[R][Ep])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
            const double mor = m*oor;
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czep = (ma->C[Z][Ep])?ma->C[Z][Ep][z+ir]:0;
              const double Crep = (ma->C[R][Ep])?ma->C[R][Ep][z+ir]:0;
              const double inveps_plus_Czep = (ma->Cdecay[Z][Ep][P]) ?
                ma->Cdecay[Z][Ep][P][z+ir] : ma->inveps[Ep][P][z+ir];
              const double inveps_plus_Crep = (ma->Cdecay[R][Ep][P]) ?
                ma->Cdecay[R][Ep][P][z+ir] : ma->inveps[Ep][P][z+ir];
              const double depz = inveps_plus_Czep*(c*(f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                                                    - Czep*f_pml[Ep][cmp][z+ir]);
              const double epr = f[Ep][cmp][z+ir] - f_pml[Ep][cmp][z+ir];
              f_pml[Ep][cmp][z+ir] += depz;
              f[Ep][cmp][z+ir] += depz +
                inveps_plus_Crep*(c*(-(f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1])) - Crep*epr);
            }
          }
      else
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Ep][cmp][z+ir] += c*ma->inveps[Ep][P][z+ir]*
                ((f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                 - (f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1]));
          }
      // Propogate Ez
      if (ma->C[R][Ez])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
          double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
          double mor = m*oor;
          const int ir = r*(v.nz()+1);
          const int irm1 = (r-1)*(v.nz()+1);
          for (int z=0;z<v.nz();z++) {
            const double Crez = ma->C[R][Ez][z+ir];
            const double inveps_plus_Crez = (ma->Cdecay[R][Ez][Z]) ?
              ma->Cdecay[R][Ez][Z][z+ir] : ma->inveps[Ez][Z][z+ir];
            const double dezr = inveps_plus_Crez*
              (c*(f[Hp][cmp][z+ir]*((int)(v.origin.r()*v.a+0.5) + r+0.5)-
                  f[Hp][cmp][z+irm1]*((int)(v.origin.r()*v.a+0.5) + r-0.5))*oor
               - Crez*f_pml[Ez][cmp][z+ir]);
            const double ezp = f[Ez][cmp][z+ir]-f_pml[Ez][cmp][z+ir];
            f_pml[Ez][cmp][z+ir] += dezr;
            f[Ez][cmp][z+ir] += dezr +
              c*ma->inveps[Ez][Z][z+ir]*(-it(cmp,f[Hr],z+ir)*mor);
          }
        }
      else
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
          double oor = 1.0/((int)(v.origin.r()*v.a + 0.5) + r);
          double mor = m*oor;
          const int ir = r*(v.nz()+1);
          const int irm1 = (r-1)*(v.nz()+1);
          for (int z=0;z<v.nz();z++)
            f[Ez][cmp][z+ir] += c*ma->inveps[Ez][Z][z+ir]*
              ((f[Hp][cmp][z+ir]*((int)(v.origin.r()*v.a+0.5) + r+0.5)-
                f[Hp][cmp][z+irm1]*((int)(v.origin.r()*v.a+0.5) + r-0.5))*oor
               - it(cmp,f[Hr],z+ir)*mor);
        }
      // Propogate Er
      if (ma->C[Z][Er])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1.0/((int)(v.origin.r()*v.a+0.5) + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czer = ma->C[Z][Er][z+ir];
              const double inveps_plus_Czer = (ma->Cdecay[Z][Er][R]) ?
                ma->Cdecay[Z][Er][R][z+ir] : ma->inveps[Er][R][z+ir];
              double derp = c*ma->inveps[Er][R][z+ir]*(it(cmp,f[Hz],z+ir)*morph);
              double erz = f[Er][cmp][z+ir] - f_pml[Er][cmp][z+ir];
              f_pml[Er][cmp][z+ir] += derp;
              f[Er][cmp][z+ir] += derp + inveps_plus_Czer*
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
              f[Er][cmp][z+ir] += c*ma->inveps[Er][R][z+ir]*
                (it(cmp,f[Hz],z+ir)*morph - (f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]));
          }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0 && v.origin.r() == 0.0) {
        for (int z=0;z<=v.nz();z++)
          f[Ez][cmp][z] += c*ma->inveps[Ez][Z][z]*
            (f[Hp][cmp][z] + it(cmp,f[Hr],z)*m);
      } else if (m == 1 && v.origin.r() == 0.0) {
        if (ma->C[Z][Ep])
          for (int z=1;z<=v.nz();z++) {
            const double Czep = ma->C[Z][Ep][z];
            const double inveps_plus_Czep = (ma->Cdecay[Z][Ep][P]) ?
              ma->Cdecay[Z][Ep][P][z] : ma->inveps[Ep][P][z];
            const double depz = inveps_plus_Czep*(c*(f[Hr][cmp][z]-f[Hr][cmp][z-1])
                                                  - Czep*f_pml[Ep][cmp][z]);
            const double epr = f[Ep][cmp][z] - f_pml[Ep][cmp][z];
            f_pml[Ep][cmp][z] += depz;
            f[Ep][cmp][z] += depz +
              c*ma->inveps[Ep][P][z]*(-f[Hz][cmp][z]*2.0);
          }
        else
          for (int z=1;z<=v.nz();z++)
            f[Ep][cmp][z] += c*ma->inveps[Ep][P][z]*
              ((f[Hr][cmp][z]-f[Hr][cmp][z-1]) - f[Hz][cmp][z]*2.0);
      } else {
        for (int r=0;r<=v.nr() && (int)(v.origin.r()*v.a+0.5) + r < m;r++) {
          const int ir = r*(v.nz()+1);
          for (int z=0;z<=v.nz();z++) f[Ep][cmp][z+ir] = 0;
          if (f_pml[Ep][cmp])
            for (int z=0;z<=v.nz();z++) f_pml[Ep][cmp][z+ir] = 0;
          for (int z=0;z<=v.nz();z++) f[Ez][cmp][z+ir] = 0;
          if (f_pml[Ez][cmp])
            for (int z=0;z<=v.nz();z++) f_pml[Ez][cmp][z+ir] = 0;
        }
      }
    }
  } else {
    abort("Unsupported dimension.\n");
  }
}

#include "config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

void fields::step_boundaries(field_type ft) {
  // First copy outgoing data to buffers...
  int *wh = new int[num_chunks];
  for (int i=0;i<num_chunks;i++) wh[i] = 0;
  for (int i=0;i<num_chunks;i++)
    for (int j=0;j<num_chunks;j++)
      if (chunks[j]->is_mine()) {
        const int pair = j+i*num_chunks;
        for (int n=0;n<comm_sizes[ft][pair];n++) {
          DOCMP {
            double stupid = *(chunks[j]->connections[ft][Outgoing][cmp][wh[j]]);
            comm_blocks[ft][pair][n*2+cmp] = stupid;
          }
          wh[j]++;
        }
      }
  // Communicate the data around!
#if 0 // This is the blocking version, which should always be safe!
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      const int pair = j+i*num_chunks;
      DOCMP {
        send(chunks[j]->n_proc(), chunks[i]->n_proc(),
             comm_blocks[ft][pair], comm_sizes[ft][pair]*2);
      }
    }
#endif
#ifdef HAVE_MPI
  am_now_working_on(MpiTime);
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
            MPI_Isend(comm_blocks[ft][pair], comm_sizes[ft][pair]*2,
                      MPI_DOUBLE, chunks[i]->n_proc(),
                      tagto[chunks[i]->n_proc()]++,
                      MPI_COMM_WORLD, &reqs[reqnum++]);
          if (chunks[i]->is_mine())
            MPI_Irecv(comm_blocks[ft][pair], comm_sizes[ft][pair]*2,
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
  finished_working();
#endif
  
  // Finally, copy incoming data to the fields themselves!
  for (int i=0;i<num_chunks;i++) wh[i] = 0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      for (int j=0;j<num_chunks;j++) {
        const int pair = j+i*num_chunks;
        for (int n=0;n<comm_sizes[ft][pair];n++) {
          complex<double> val = chunks[i]->connection_phases[ft][wh[i]]*
            complex<double>(comm_blocks[ft][pair][n*2],
                            comm_blocks[ft][pair][n*2+1]);
          *(chunks[i]->connections[ft][Incoming][0][wh[i]]) = real(val);
          *(chunks[i]->connections[ft][Incoming][1][wh[i]]) = imag(val);
          wh[i]++;
        }
      }
  delete[] wh;
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
  for (int c=0;c<10;c++)
    if (f[c][0] && is_magnetic((component)c)) {
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
  for (int c=0;c<10;c++)
    if (f[c][0] && is_electric((component)c)) {
      f[c][0][s->i] += real(A*s->A[c]);
      if (!is_real) f[c][1][s->i] += imag(A*s->A[c]);
    }
  step_e_source(s->next, time);
}

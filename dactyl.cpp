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
#include <complex>

#include "dactyl.h"
#include "dactyl_internals.h"

fields::~fields() {
  delete ma;
  DOCMP {
    for (int i=0;i<10;i++) delete[] f[i][cmp];
    for (int i=0;i<10;i++) delete[] f_backup[i][cmp];
    for (int i=0;i<10;i++) delete[] f_pml[i][cmp];
    for (int i=0;i<10;i++) delete[] f_backup_pml[i][cmp];
    delete[] e_connection_sinks[cmp];
    delete[] e_connection_sources[cmp];
    delete[] h_connection_sinks[cmp];
    delete[] h_connection_sources[cmp];
  }
  delete[] e_phases;
  delete[] h_phases;
  delete h_sources;
  delete e_sources;
  delete bands;
  delete pol;
  delete olpol;
}

void fields::use_bloch(double tk) {
  k = tk;
  cosknz = cos(k*2*pi*inva*v.nz());
  sinknz = sin(k*2*pi*inva*v.nz());
  eiknz = complex<double>(cosknz, sinknz);
  complex<double> emiknz = complex<double>(cosknz, -sinknz);
  if (v.dim == d1) {
    delete[] e_phases;
    delete[] h_phases;
    num_h_connections = 0;
    num_e_connections = 1;
    e_phases = new complex<double>[num_e_connections];
    e_phases[0] = eiknz;
    DOCMP {
      delete[] e_connection_sinks[cmp];
      delete[] e_connection_sources[cmp];
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp][0] = &f[Ex][cmp][0];
      e_connection_sinks[cmp][0] = &f[Ex][cmp][v.nz()];
    }
  } else if (v.dim == dcyl) {
    num_e_connections = 3*(v.nr()+1);
    num_h_connections = 3*(v.nr()+1);
    delete[] e_phases;
    delete[] h_phases;
    e_phases = new complex<double>[num_e_connections];
    h_phases = new complex<double>[num_h_connections];
    DOCMP {
      delete[] e_connection_sinks[cmp];
      delete[] e_connection_sources[cmp];
      delete[] h_connection_sinks[cmp];
      delete[] h_connection_sources[cmp];
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      h_connection_sinks[cmp] = new (double *)[num_h_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      h_connection_sources[cmp] = new (double *)[num_h_connections];
      int econ = 0;
      int hcon = 0;
      for (int r=0;r<=v.nr();r++) {
        int right = v.nz() + r*(v.nz()+1);
        int left = 0 + r*(v.nz()+1);
        e_connection_sources[cmp][econ] = &f[Ep][cmp][right];
        e_connection_sinks[cmp][econ] = &f[Ep][cmp][left];
        e_phases[econ] = emiknz;
        econ++;
        e_connection_sources[cmp][econ] = &f[Er][cmp][right];
        e_connection_sinks[cmp][econ] = &f[Er][cmp][left];
        e_phases[econ] = emiknz;
        econ++;
        e_connection_sources[cmp][econ] = &f[Ez][cmp][left];
        e_connection_sinks[cmp][econ] = &f[Ez][cmp][right];
        e_phases[econ] = eiknz;
        econ++;

        h_connection_sources[cmp][hcon] = &f[Hz][cmp][right];
        h_connection_sinks[cmp][hcon] = &f[Hz][cmp][left];
        h_phases[hcon] = emiknz;
        hcon++;
        h_connection_sources[cmp][hcon] = &f[Hr][cmp][left];
        h_connection_sinks[cmp][hcon] = &f[Hr][cmp][right];
        h_phases[hcon] = eiknz;
        hcon++;
        h_connection_sources[cmp][hcon] = &f[Hp][cmp][left];
        h_connection_sinks[cmp][hcon] = &f[Hp][cmp][right];
        h_phases[hcon] = eiknz;
        hcon++;
        // FIXME: Need to add PML connections here!
      }
    }
  } else {
    printf("Unsupported dimension?!\n");
    exit(1);
  }
}

double zero = 0.0;

fields::fields(const mat *the_ma, int tm) {
  int r, z;
  ma = new mat(the_ma);
  verbosity = 0;
  v = ma->v;
  outdir = ma->outdir;
  m = tm;
  phasein_time = 0;
  new_ma = NULL;
  bands = NULL;
  k = -1;
  a = ma->a;
  inva = 1.0/a;
  preferred_fmax = 2.5; // Some sort of reasonable maximum
                        // frequency... (assuming a has a value on the
                        // order of your frequency).
  t = 0;
  pol = polarization::set_up_polarizations(ma);
  olpol = polarization::set_up_polarizations(ma);
  h_sources = e_sources = NULL;
  DOCMP {
    for (int i=0;i<10;i++) f[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_backup[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_pml[i][cmp] = NULL;
    for (int i=0;i<10;i++) f_backup_pml[i][cmp] = NULL;

    for (int i=0;i<10;i++) if (v.has_field((component)i))
      f[i][cmp] = new double[v.ntot()];
    for (int i=0;i<10;i++) if (v.has_field((component)i)) {
      f_pml[i][cmp] = new double[v.ntot()];
      if (f_pml[i][cmp] == NULL) {
        printf("Out of memory!\n");
        exit(1);
      }
    }
  }
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c))
        for (int i=0;i<v.ntot();i++)
          f[c][cmp][i] = 0.0;
    // Now for pml extra fields...
    for (int c=0;c<10;c++)
      if (v.has_field((component)c))
        for (int i=0;i<v.ntot();i++)
          f_pml[c][cmp][i] = 0.0;
  }
  if (v.dim == d1) {
    // Set up by default with metallic boundary conditions.
    num_h_connections = 0;
    num_e_connections = 1;
    e_phases = new complex<double>[num_e_connections];
    e_phases[0] = 0.0;
    h_phases = NULL;
    DOCMP {
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      e_connection_sources[cmp][0] = &zero;
      e_connection_sinks[cmp][0] = &f[Ex][cmp][v.nz()];
      h_connection_sinks[cmp] = NULL;
      h_connection_sources[cmp] = NULL;
    }
  } else if (v.dim == dcyl) {
    // Set up by default with metallic boundary conditions.
    num_e_connections = 2*(v.nr()+1);
    num_h_connections =   (v.nr()+1);
    if (f_pml[Ep][0]) num_e_connections += v.nr()+1;
    if (f_pml[Er][0]) num_e_connections += v.nr()+1;
    if (f_pml[Hz][0]) num_h_connections += v.nr()+1;
    e_phases = new complex<double>[num_e_connections];
    h_phases = new complex<double>[num_h_connections];
    DOCMP {
      e_connection_sinks[cmp] = new (double *)[num_e_connections];
      h_connection_sinks[cmp] = new (double *)[num_h_connections];
      e_connection_sources[cmp] = new (double *)[num_e_connections];
      h_connection_sources[cmp] = new (double *)[num_h_connections];
      int econ = 0;
      int hcon = 0;
      for (int i=0;i<num_e_connections;i++) e_phases[i] = 1.0;
      for (int i=0;i<num_h_connections;i++) h_phases[i] = 1.0;
      for (int r=0;r<=v.nr();r++) {
        int i = v.nz() + r*(v.nz()+1);
        e_connection_sources[cmp][econ] = &zero;
        e_connection_sinks[cmp][econ] = &f[Ep][cmp][i];
        econ++;
        e_connection_sources[cmp][econ] = &zero;
        e_connection_sinks[cmp][econ] = &f[Er][cmp][i];
        econ++;
        h_connection_sources[cmp][hcon] = &zero;
        h_connection_sinks[cmp][hcon] = &f[Hz][cmp][i];
        hcon++;
        // I have to add PML connections here too:
        if (f_pml[Ep][cmp]) {
          e_connection_sources[cmp][econ] = &zero;
          e_connection_sinks[cmp][econ] = &f_pml[Ep][cmp][i];
          econ++;
        }
        if (f_pml[Er][cmp]) {
          e_connection_sources[cmp][econ] = &zero;
          e_connection_sinks[cmp][econ] = &f_pml[Er][cmp][i];
          econ++;
        }
        if (f_pml[Hz][0]) {
          h_connection_sources[cmp][hcon] = &zero;
          h_connection_sinks[cmp][hcon] = &f_pml[Hz][cmp][i];
          hcon++;
        }
      }
    }
  }
}

void fields::use_real_sources() {
  if (e_sources) e_sources->use_real_sources();
  if (h_sources) h_sources->use_real_sources();
}

void src::use_real_sources() {
  is_real = 1;
  if (next) next->use_real_sources();
}

int fields::phase_in_material(const mat *newma, double time) {
  new_ma = newma;
  phasein_time = (int) (time*a/c);
  if (phasein_time == 0) phasein_time = 1;
  printf("I'm going to take %d time steps to phase in the material.\n", phasein_time);
  return phasein_time;
}

int fields::is_phasing() {
  return phasein_time > 0;
}

complex<double> src::get_amplitude_at_time(int t) const {
  double envelope = get_envelope_at_time(t);
  if (envelope == 0.0)
    return 0.0;
  double tt = t - peaktime;
  if (is_real) return real( (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope );
  else return (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope;
}

double src::get_envelope_at_time(int t) const {
  double tt = t - peaktime;
  if (is_continuous && tt > 0) {
    return 1.0;
  } else if (fabs(tt) > cutoff) {
    return 0.0;
  } else {
    return exp(-tt*tt/(2*width*width));
  }
}

static double integrate_envelope(const src *s) {
  if (s == NULL) {
    printf("Bad arg to integrate_envelope!\n");
    exit(1);
  }
  double sofar = 0.0;
  for (int t=(int)s->peaktime-s->cutoff;t<(1<<30);t++) {
    double e = s->get_envelope_at_time(t);
    sofar += e;
    if (e == 0) break; // But here if there is a source that starts late,
                       // or a source that never stops.
  }
  return sofar;
}

static complex<double> integrate_source(const src *s) {
  if (s == NULL) {
    printf("Bad arg to integrate_source!\n");
    exit(1);
  }
  complex<double> sofar = 0.0;
  for (int t=0;1<<30;t++) {
    complex<double> A = s->get_amplitude_at_time(t);
    sofar += A;
    if (A == 0) break; // But here if there is a source that starts late,
                       // or a source that never stops.
  }
  return sofar;
}

void fields::add_point_source(component whichf, double freq, double width, double peaktime,
                              double cutoff, const vec &p, complex<double> amp,
                              int is_c) {
  // FIXME this really should call an interpolation routine...
  if (p.dim != v.dim) {
    printf("Error:  source doesn't have right dimensions!\n");
    exit(1);
  } else if (!v.has_field(whichf)) {
    printf("Error:  source component %s is invalid.\n", component_name(whichf));
    exit(1);
  }
  int theindex = v.index(whichf, p);
  add_indexed_source(whichf, freq, width, peaktime, cutoff, theindex, amp, is_c);
}

void fields::add_plane_source(double freq, double width, double peaktime,
                              double cutoff, double envelope (const vec &),
                              const vec &p, const vec &norm,
                              int is_c) {
  if (v.dim == dcyl) {
    // We ignore norm in this case...
    if (m != 1) {
      printf("Can only use plane source with m == 1!\n");
      exit(1);
    }
    const complex<double> I = complex<double>(0,1);
    const double z = p.z();
    const double eps = sqrt(ma->eps[(int)(z+0.5)]);
    for (int ir=0;ir<v.nr();ir++) {
      {
        const double r = ir*inva;
        // E_phi
        add_point_source(Ep, freq, width, peaktime, cutoff, vec(r,z),
                         envelope(vec(r,z)), is_c);        
        // iH_r = d(rH_phi)/dr
        const double slope = ((r+0.5)*envelope(vec(r+0.5*inva,z)) -
                              (r-0.5)*envelope(vec(r-0.5*inva,z)))*a;
        add_point_source(Hr, freq, width, peaktime, cutoff, vec(r,z), -eps*slope, is_c);
      }
      {
        const double r = (ir+0.5)*inva;
        const double sc = (ir == 0)?0.5:1.0;
        // iE_r = d(rE_phi)/dr
        const double slope = ((r+0.5)*envelope(vec(r+0.5*inva,z)) -
                              (r-0.5)*envelope(vec(r-0.5*inva,z)))*a;
        add_point_source(Er, freq, width, peaktime, cutoff, vec(r,z), -I*sc*slope, is_c);
        // H_phi
        add_point_source(Hp, freq, width, peaktime, cutoff, vec(r,z),
                         -I*eps*sc*envelope(vec(r,z)), is_c);
      }
    }
  } else if (v.dim == d1) {
    const double z = p.z();
    const double eps = sqrt(ma->eps[(int)(z+0.5)]);
    add_point_source(Ex, freq, width, peaktime, cutoff, vec(z), envelope(vec(z)), is_c);
    add_point_source(Hy, freq, width, peaktime, cutoff, vec(z), envelope(vec(z))*eps, is_c);
  } else {
    printf("Can't use plane source in this number of dimensions.\n");
    exit(1);
  }
}

void fields::add_indexed_source(component whichf, double freq, double width,
                                double peaktime, int cutoff, int theindex, 
                                complex<double> amp, int is_c) {
  if (theindex >= v.ntot() || theindex < 0) {
    printf("Error:  source is outside of cell! (%d)\n", theindex);
    exit(1);
  }
  src *tmp = new src;
  tmp->freq = freq*c*inva;
  tmp->width = width/tmp->freq; // this is now time width
  for (int com=0;com<10;com++) tmp->A[com] = 0;
  tmp->A[whichf] = amp;
  tmp->i = theindex;
  tmp->amp_shift = 0.0;
  tmp->is_real = 0;
  tmp->is_continuous = is_c;
  if (whichf == Hz || whichf == Hr || whichf == Hp) {
    tmp->next = h_sources;
    h_sources = tmp;
  } else {
    tmp->next = e_sources;
    e_sources = tmp;
  }
  tmp->cutoff = 1+ (int)(cutoff*tmp->width);
  tmp->peaktime = peaktime*a/c;
  if (peaktime <= 0.0) tmp->peaktime = t+tmp->cutoff;
  // Apply a shift so that we won't end up with a static polarization when
  // the source is gone:
  if (is_c) tmp->amp_shift = 0.0;
  else tmp->amp_shift = integrate_source(tmp)/integrate_envelope(tmp);
}

void fields::step_right() {
  t += 1;

  phase_material();

  step_h_right();
  step_h_source(h_sources);
  step_h_boundaries();

  step_e_right();
  step_e_source(e_sources);
  step_e_boundaries();

  step_e_polarization();
  step_polarization_itself();
}

void fields::step() {
  t += 1;

  phase_material();

  step_h();
  step_h_source(h_sources);
  step_h_boundaries();

  prepare_step_polarization_energy();
  half_step_polarization_energy();
  step_e();
  step_e_source(e_sources);
  step_e_boundaries();
  half_step_polarization_energy();

  step_e_polarization();
  step_polarization_itself();
}

void fields::phase_material() {
  if (new_ma && phasein_time > 0) {
    // First multiply the electric field by epsilon...
    DOCMP {
      for (int c=0;c<10;c++)
        if (v.has_field((component)c) && is_electric((component)c))
          for (int i=0;i<v.ntot();i++)
            f[c][cmp][i] /= ma->inveps[c][i];
      for (int c=0;c<10;c++)
        if (v.has_field((component)c) && is_electric((component)c))
          for (int i=0;i<v.ntot();i++)
            f_pml[c][cmp][i] /= ma->inveps[c][i];
    }
    // Then change epsilon...
    ma->mix_with(new_ma, 1.0/phasein_time);
    phasein_time--;
    // Finally divide the electric field by the new epsilon...
    DOCMP {
      for (int c=0;c<10;c++)
        if (v.has_field((component)c) && is_electric((component)c))
          for (int i=0;i<v.ntot();i++)
            f[c][cmp][i] *= ma->inveps[c][i];
      for (int c=0;c<10;c++)
        if (v.has_field((component)c) && is_electric((component)c))
          for (int i=0;i<v.ntot();i++)
            f_pml[c][cmp][i] *= ma->inveps[c][i];
    }
  }
}

inline double it(int cmp, double *(f[2]), int ind) { return (1-2*cmp)*f[1-cmp][ind]; }

inline int rstart_0(const volume &v, int m) {
  return (int) max(0.0, m - v.origin.r()*v.a - 1);
}
inline int rstart_1(const volume &v, int m) {
  return (int) max(1.0, m - v.origin.r()*v.a);
}

void fields::step_h_right() {
  const volume v = this->v;
  if (v.dim == d1) {
    DOCMP {
      for (int z=0;z<v.nz();z++)
        f[Hy][cmp][z] = f[Ex][cmp][z];
    }
  } else if (v.dim == dcyl) {
    DOCMP {
      for (int r=rstart_1(v,m);r<=v.nr();r++) {
        const int ir = r*(v.nz()+1);
        for (int z=0;z<v.nz();z++) f[Hr][cmp][z+ir] = f[Ep][cmp][z+ir];
      }
      for (int r=rstart_0(v,m);r<v.nr();r++) {
        const int ir = r*(v.nz()+1);
        const int irp1 = (r+1)*(v.nz()+1);
        for (int z=0;z<v.nz();z++) f[Hp][cmp][z+ir] = f[Ez][cmp][z+irp1];
      }
      for (int r=rstart_0(v,m);r<v.nr();r++) {
        const int ir = r*(v.nz()+1);
        for (int z=1;z<=v.nz();z++) f[Hz][cmp][z+ir] = f[Er][cmp][z+ir];
      }
    }
  }
}

void fields::step_e_right() {
  const volume v = this->v;
  if (v.dim == d1) {
    DOCMP {
      for (int z=0;z<v.nz();z++)
        f[Hy][cmp][z] = f[Ex][cmp][z];
    }
  } else if (v.dim == dcyl) {
    DOCMP {
      for (int r=rstart_1(v,m);r<=v.nr();r++) {
        const int ir = r*(v.nz()+1);
        const int irm1 = (r-1)*(v.nz()+1);
        for (int z=1;z<=v.nz();z++) f[Ep][cmp][z+ir] = f[Hz][cmp][z+irm1];
      }
      for (int r=rstart_0(v,m);r<v.nr();r++) {
        const int ir = r*(v.nz()+1);
        for (int z=1;z<=v.nz();z++) f[Er][cmp][z+ir] = f[Hp][cmp][z+ir-1];
      }
      for (int r=rstart_1(v,m);r<=v.nr();r++) {
        const int ir = r*(v.nz()+1);
        for (int z=0;z<v.nz();z++) f[Ez][cmp][z+ir] = f[Hr][cmp][z+ir];
      }
    }
  }
}

void fields::step_h() {
  const volume v = this->v;
  if (v.dim == d1) {
    DOCMP {
      if (ma->Cmain[Hy])
        for (int z=0;z<v.nz();z++) {
          const double C = ma->Cmain[Hy][z];
          const double ooop_C = 1.0/(1.0+0.5*C);
          f[Hy][cmp][z] += ooop_C*(-c*(f[Ex][cmp][z+1] - f[Ex][cmp][z])
                                   - C*f[Hy][cmp][z]);
        }
      else
        for (int z=0;z<v.nz();z++)
          f[Hy][cmp][z] += -c*(f[Ex][cmp][z+1] - f[Ex][cmp][z]);
    }
  } else if (v.dim == dcyl) {
    for (int cmp=0;cmp<=1;cmp++) {
      // Propogate Hr
      if (ma->Cother[Hr])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/(v.origin.r() + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Czhr = ma->Cother[Hr][z+ir];
              const double ooop_Czhr = 1.0/(1.0+0.5*Czhr);
              double dhrp = c*(-it(cmp,f[Ez],z+ir)*mor);
              double hrz = f[Hr][cmp][z+ir] - f_pml[Hr][cmp][z+ir];
              f_pml[Hr][cmp][z+ir] += dhrp;
              f[Hr][cmp][z+ir] += dhrp +
                ooop_Czhr*(c*(f[Ep][cmp][z+ir+1]-f[Ep][cmp][z+ir]) - Czhr*hrz);
            }
          }
      else
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            double oor = 1.0/(v.origin.r() + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Hr][cmp][z+ir] += c*
                ((f[Ep][cmp][z+ir+1]-f[Ep][cmp][z+ir]) - it(cmp,f[Ez],z+ir)*mor);
          }
      // Propogate Hp
      if (ma->Cmain[Hp] || ma->Cother[Hp])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1/(v.origin.r() + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++) {
              const double Czhp = (ma->Cmain[Hp])?ma->Cmain[Hp][z+ir]:0;
              const double Crhp = (ma->Cother[Hp])?ma->Cother[Hp][z+ir]:0;
              const double ooop_Czhp = 1.0/(1.0+0.5*Czhp);
              const double ooop_Crhp = 1.0/(1.0+0.5*Crhp);
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
            double oorph = 1/(v.origin.r() + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=0;z<v.nz();z++)
              f[Hp][cmp][z+ir] += c*
                ((f[Ez][cmp][z+irp1]-f[Ez][cmp][z+ir])
                 - (f[Er][cmp][z+ir+1]-f[Er][cmp][z+ir]));
          }
      // Propogate Hz
      if (ma->Cother[Hz])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
          double oorph = 1/(v.origin.r() + r+0.5);
          double morph = m*oorph;
          const int ir = r*(v.nz()+1);
          const int irp1 = (r+1)*(v.nz()+1);
          for (int z=1;z<=v.nz();z++) {
            const double Crhz = ma->Cother[Hz][z+ir];
            const double ooop_Crhz = 1.0/(1.0+0.5*Crhz);
            const double dhzr =
              ooop_Crhz*(-c*(f[Ep][cmp][z+irp1]*(r+1.)-f[Ep][cmp][z+ir]*r)*oorph
                         - Crhz*f_pml[Hz][cmp][z+ir]);
            const double hzp = f[Hz][cmp][z+ir] - f_pml[Hz][cmp][z+ir];
            f_pml[Hz][cmp][z+ir] += dhzr;
            f[Hz][cmp][z+ir] += dhzr + c*(it(cmp,f[Er],z+ir)*morph);
          }
        }
      else
        for (int r=rstart_0(v,m);r<v.nr();r++) {
          double oorph = 1/(v.origin.r() + r+0.5);
          double morph = m*oorph;
          const int ir = r*(v.nz()+1);
          const int irp1 = (r+1)*(v.nz()+1);
          for (int z=1;z<=v.nz();z++)
            f[Hz][cmp][z+ir] += c*
              (it(cmp,f[Er],z+ir)*morph
               - (f[Ep][cmp][z+irp1]*(r+1.)-f[Ep][cmp][z+ir]*r)*oorph);
        }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0) {
        // Nothing needed for H.
      } else if (m == 1) {
        if (ma->Cmain[Hr])
          for (int z=0;z<v.nz();z++) {
            const double Czhr = ma->Cmain[Hr][z];
            const double ooop_Czhr = 1.0/(1.0+0.5*Czhr);
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
        for (int r=0;r<=v.nr() && v.origin.r() + r < m;r++) {
          const int ir = r*(v.nz()+1);
          for (int z=0;z<=v.nz();z++) f[Hr][cmp][z+ir] = 0;
          if (f_pml[Hr][cmp])
            for (int z=0;z<=v.nz();z++) f_pml[Hr][cmp][z+ir] = 0;
        }
      }
    }
  }
}

void fields::step_h_boundaries() {
  for (int i=0;i<num_h_connections;i++) {
    complex<double> val = h_phases[i]*
      complex<double>(*h_connection_sources[0][i],*h_connection_sources[1][i]);
    *h_connection_sinks[0][i] = real(val);
    *h_connection_sinks[1][i] = imag(val);
  }
}

void fields::step_e() {
  const volume v = this->v;
  if (v.dim == d1) {
    DOCMP {
      if (ma->Cmain[Ex])
        for (int z=1;z<=v.nz();z++) {
          const double C = ma->Cmain[Ex][z];
          const double inveps_plus_C = ma->inveps[Ex][z]/(1.0+0.5*ma->inveps[Ex][z]*C);
          f[Ex][cmp][z] += inveps_plus_C*(-c*(f[Hy][cmp][z] - f[Hy][cmp][z-1])
                                          - C*f[Ex][cmp][z]);
        }
      else 
        for (int z=1;z<=v.nz();z++)
          f[Ex][cmp][z] += -c*ma->inveps[Ex][z]*(f[Hy][cmp][z] - f[Hy][cmp][z-1]);
    }
  } else if (v.dim == dcyl) {
    for (int cmp=0;cmp<=1;cmp++) {
      // Propogate Ep
      if (ma->Cmain[Ep] || ma->Cother[Ep])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
            const double oor = 1.0/(v.origin.r() + r);
            const double mor = m*oor;
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czep = (ma->Cmain[Ep])?ma->Cmain[Ep][z+ir]:0;
              const double Crep = (ma->Cmain[Er])?ma->Cmain[Er][z+ir]:0;
              const double inveps_plus_Czep = ma->inveps[Ep][z+ir]/
                (1+.5*ma->inveps[Ep][z+ir]*Czep);
              const double inveps_plus_Crep = ma->inveps[Ep][z+ir]/
                (1+.5*ma->inveps[Ep][z+ir]*Crep);
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
            double oor = 1.0/(v.origin.r() + r);
            double mor = m*oor;
            const int ir = r*(v.nz()+1);
            const int irm1 = (r-1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Ep][cmp][z+ir] += c*ma->inveps[Ep][z+ir]*
                ((f[Hr][cmp][z+ir]-f[Hr][cmp][z+ir-1])
                 - (f[Hz][cmp][z+ir]-f[Hz][cmp][z+irm1]));
          }
      // Propogate Ez
      if (ma->Cother[Ez])
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
          double oor = 1.0/(v.origin.r() + r);
          double mor = m*oor;
          const int ir = r*(v.nz()+1);
          const int irm1 = (r-1)*(v.nz()+1);
          for (int z=0;z<v.nz();z++) {
            const double Crez = ma->Cother[Ez][z+ir];
            const double inveps_plus_Crez = ma->inveps[Ez][z+ir]/
              (1+.5*ma->inveps[Ez][z+ir]*Crez);

            const double dezr = inveps_plus_Crez*
              (c*(f[Hp][cmp][z+ir]*(r+0.5)-f[Hp][cmp][z+irm1]*(r-0.5))*oor
               - Crez*f_pml[Ez][cmp][z+ir]);
            const double ezp = f[Ez][cmp][z+ir]-f_pml[Ez][cmp][z+ir];
            f_pml[Ez][cmp][z+ir] += dezr;
            f[Ez][cmp][z+ir] += dezr +
              c*ma->inveps[Ez][z+ir]*(-it(cmp,f[Hr],z+ir)*mor);
          }
        }
      else
        for (int r=rstart_1(v,m);r<=v.nr();r++) {
          double oor = 1.0/(v.origin.r() + r);
          double mor = m*oor;
          const int ir = r*(v.nz()+1);
          const int irm1 = (r-1)*(v.nz()+1);
          for (int z=0;z<v.nz();z++)
            f[Ez][cmp][z+ir] += c*ma->inveps[Ez][z+ir]*
              ((f[Hp][cmp][z+ir]*(r+0.5)-f[Hp][cmp][z+irm1]*(r-0.5))*oor
               - it(cmp,f[Hr],z+ir)*mor);
        }
      // Propogate Er
      if (ma->Cother[Er])
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1/(v.origin.r() + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++) {
              const double Czer = ma->Cother[Er][z+ir];
              const double inveps_plus_Czer = ma->inveps[Er][z+ir]/
                (1+.5*ma->inveps[Er][z+ir]*Czer);
              double derp = c*ma->inveps[Er][z+ir]*(it(cmp,f[Hz],z+ir)*morph);
              double erz = f[Er][cmp][z+ir] - f_pml[Er][cmp][z+ir];
              f_pml[Er][cmp][z+ir] += derp;
              f[Er][cmp][z+ir] += derp + inveps_plus_Czer*
                (-c*(f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]) - Czer*erz);
            }
          }
      else
        for (int r=rstart_0(v,m);r<v.nr();r++) {
            double oorph = 1/(v.origin.r() + r+0.5);
            double morph = m*oorph;
            const int ir = r*(v.nz()+1);
            const int irp1 = (r+1)*(v.nz()+1);
            for (int z=1;z<=v.nz();z++)
              f[Er][cmp][z+ir] += c*ma->inveps[Er][z+ir]*
                (it(cmp,f[Hz],z+ir)*morph - (f[Hp][cmp][z+ir]-f[Hp][cmp][z+ir-1]));
          }
      // Deal with annoying r==0 boundary conditions...
      if (m == 0) {
        for (int z=1;z<=v.nz();z++)
          f[Ez][cmp][z] += c*ma->inveps[Ez][z]*(f[Hp][cmp][z] + it(cmp,f[Hr],z)*m);
      } else if (m == 1) {
        if (ma->Cmain[Ep])
          for (int z=1;z<=v.nz();z++) {
            const double Czep = ma->Cmain[Ep][z];
            const double inveps_plus_Czep = ma->inveps[Ep][z]/
              (1+.5*ma->inveps[Ep][z]*Czep);
            const double depz = inveps_plus_Czep*(c*(f[Hr][cmp][z]-f[Hr][cmp][z-1])
                                                  - Czep*f_pml[Ep][cmp][z]);
            const double epr = f[Ep][cmp][z] - f_pml[Ep][cmp][z];
            f_pml[Ep][cmp][z] += depz;
            f[Ep][cmp][z] += depz +
              c*ma->inveps[Ep][z]*(-f[Hz][cmp][z]*2.0);
          }
        else
          for (int z=1;z<=v.nz();z++)
            f[Ep][cmp][z] += c*ma->inveps[Ep][z]*
              ((f[Hr][cmp][z]-f[Hr][cmp][z-1]) - f[Hz][cmp][z]*2.0);
      } else {
        for (int r=0;r<=v.nr() && v.origin.r() + r < m;r++) {
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
    printf("Unsupported dimension.\n");
    exit(1);
  }
}

void fields::step_e_boundaries() {
  for (int i=0;i<num_e_connections;i++) {
    complex<double> val = e_phases[i]*
      complex<double>(*e_connection_sources[0][i],*e_connection_sources[1][i]);
    *e_connection_sinks[0][i] = real(val);
    *e_connection_sinks[1][i] = imag(val);
  }
}

void fields::step_h_source(const src *s) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(t);
  if (A == 0.0) {
    step_h_source(s->next);
    return;
  }
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_magnetic((component)c)) {
      f[c][0][s->i] += imag(A*s->A[c]);
      f[c][1][s->i] += real(A*s->A[c]);
    }
  step_h_source(s->next);
}

void fields::step_e_source(const src *s) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(t);
  if (A == 0.0) {
    step_e_source(s->next);
    return;
  }
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c)) {
      f[c][0][s->i] += imag(A*s->A[c]);
      f[c][1][s->i] += real(A*s->A[c]);
    }
  step_e_source(s->next);
}

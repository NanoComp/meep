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

#include "tidod.h"
#include "tidod_internals.h"

fields_1d::~fields_1d() {
  delete ma;
  DOCMP {
    delete[] hy[cmp];
    delete[] ex[cmp];
    delete[] backup_hy[cmp];
  }
  if (bands) delete bands;
  delete pol;
  delete olpol;
}

void fields_1d::use_bloch(double tk) {
  if (npmlz) npmlz = 0;
  k = tk;
  cosknz = cos(k*2*pi*inva*nz);
  sinknz = sin(k*2*pi*inva*nz);
  eiknz = complex<double>(cosknz, sinknz);
}

fields_1d::fields_1d(const mat_1d *the_ma) {
  int z;
  ma = new mat_1d(the_ma);
  verbosity = 0;
  nz = ma->nz;
  outdir = ma->outdir;
  phasein_time = 0;
  new_ma = NULL;
  bands = NULL;
  k = -1;
  a = ma->a;
  inva = 1.0/a;
  t = 0;
  pol = polarization_1d::set_up_polarizations(ma);
  olpol = polarization_1d::set_up_polarizations(ma);
  h_sources = e_sources = NULL;
  DOCMP {
    hy[cmp] = new double[nz+1];
    ex[cmp] = new double[nz+1];
    backup_hy[cmp] = NULL;
  }
  if (ex[1] == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  npmlz = ma->npmlz;
  DOCMP {
    for (z=0;z<nz+1;z++) CM(hy,z) = 0.0;
    for (z=0;z<nz+1;z++) CM(ex,z) = 0.0;
  }
}

void fields_1d::use_real_sources() {
  if (e_sources) e_sources->use_real_sources();
  if (h_sources) h_sources->use_real_sources();
}

void src_1d::use_real_sources() {
  is_real = 1;
  if (next) next->use_real_sources();
}

void fields_1d::initialize_with_n(int ntot) {
  for (int n=0;n<ntot;n++) initialize_with_nth(n);
}

inline double expi(int cmp, double x) {
  return (cmp) ? cos(x) : sin(x);
}

void fields_1d::initialize_with_nth(int np1) {
  const int n = np1 - 1;
  double kk = k*2*pi*inva + n*2*pi/nz;
  DOCMP {
    for (int z=0;z<=nz;z++) {
      CM(ex,z) += expi(cmp, kk*z);
      CM(ex,z) += expi(cmp, kk*z+pi/2);
    }
  }
}

int fields_1d::phase_in_material(const mat_1d *newma, double time) {
  new_ma = newma;
  phasein_time = (int) (time*a/c);
  if (phasein_time == 0) phasein_time = 1;
  printf("I'm going to take %d time steps to phase in the material.\n", phasein_time);
  return phasein_time;
}

int fields_1d::is_phasing() {
  return phasein_time > 0;
}

complex<double> src_1d::get_amplitude_at_time(int t) const {
  double envelope = get_envelope_at_time(t);
  if (envelope == 0.0)
    return 0.0;
  double tt = t - peaktime;
  if (is_real) return real( (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope );
  else return (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope;
}

double src_1d::get_envelope_at_time(int t) const {
  double tt = t - peaktime;
  if (is_continuous && tt > 0) {
    return 1.0;
  } else if (fabs(tt) > cutoff) {
    return 0.0;
  } else {
    return exp(-tt*tt/(2*width*width));
  }
}

static double integrate_envelope(const src_1d *s) {
  if (s == NULL) {
    printf("Bad arg to integrate_envelope!\n");
    exit(1);
  }
  double sofar = 0.0;
  for (int t=1+s->peaktime-s->cutoff;t<(1<<30);t++) {
    double e = s->get_envelope_at_time(t);
    sofar += e;
    if (e == 0) break; // But here if there is a source that starts late,
                       // or a source that never stops.
  }
  return sofar;
}

static complex<double> integrate_source(const src_1d *s) {
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

void fields_1d::add_src_pt(int z, complex<double> amp,
                           double freq, double width, double peaktime,
                           double cutoff, int is_h, int is_contin) {
  if (amp == 0) return;
  if (z >= nz || z < 0) {
    printf("Error:  source is outside of cell!\n");
    exit(1);
  }
  src_1d *tmp = new src_1d;
  tmp->freq = freq*c*inva;
  tmp->width = width/tmp->freq; // this is now time width
  tmp->amp = amp;
  tmp->z = z;
  tmp->amp_shift = 0.0;
  tmp->is_real = 0;
  tmp->is_continuous = is_contin;
  if (is_h) {
    tmp->next = h_sources;
    h_sources = tmp;
  } else {
    tmp->next = e_sources;
    e_sources = tmp;
  }
  tmp->cutoff = 1+ (int)(cutoff*tmp->width);
  tmp->peaktime = peaktime*a;
  if (peaktime <= 0.0) tmp->peaktime = t+tmp->cutoff;
  // Apply a shift so that we won't end up with a static polarization when
  // the source is gone:
  if (is_contin) tmp->amp_shift = 0.0;
  else tmp->amp_shift = integrate_source(tmp)/integrate_envelope(tmp);
}

void fields_1d::find_source_z_position(double z, double shift, int *z1, int *z2,
                                       complex<double> *amp1, complex<double> *amp2) {
    
  int z_floor = (int)floor(z - shift);
  int z_ceil = z_floor + 1;
  
  // The bulk case:
  *z1 = z_floor;
  *amp1 = (z_ceil + shift) - z;
  *z2 = z_ceil;
  *amp2 = z - (z_floor + shift);

  if (npmlz > 0) { 
    if (z == 0.0) { // place source on first lattice point that's not in the PML
      *z1 = npmlz + 1 - (int) (2.0*shift);
      *amp1 = 1.0;
      *z2 = 0;
      *amp2 = 0.0;
    } else if (z == nz) {
      *z1 = nz - 1 - npmlz;
      *amp1 = 1.0;
      *z2 = 0;
      *amp2 = 0.0;
    }
  } else if (k >= 0.0) { // Periodic boundary conditions...
    while (*z1 < 0) {
      *amp1 *= eiknz;
      *z1 += nz;
    }
    while (*z2 < 0) {
      *amp2 *= eiknz;
      *z2 += nz;
    }
    while (*z1 >= nz) {
      *z1 -= nz;
      *amp1 /= eiknz;
    }
    while (*z2 >= nz) {
      *z2 -= nz;
      *amp2 /= eiknz;
    }
  } else if (k < 0) { // metal
    if (z <= (1 - shift)) { // if z_floor is on metal or out of the cell
      *z1 = z_ceil;
      *amp1 = 1.0;
      *z2 = 0;
      *amp2 = 0.0;
    } else if (z >= (nz - 1 + shift)) { // if z_ceil is on metal or out of the cell
      *z1 = z_floor;
      *amp1 = 1.0;
      *z2 = 0;
      *amp2 = 0.0;
    }
  }
}

void fields_1d::add_source(component_1d whichf, double freq, double width, double peaktime,
                           double cutoff, double z, int is_c) {
  double zshift = 0.0;
  if (whichf == Hy) zshift = 0.5;
  int is_h = whichf == Hy;
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, zshift, &z1, &z2, &amp1, &amp2);
  add_src_pt(z1, amp1, freq, width, peaktime, cutoff, is_h, is_c);
  add_src_pt(z2, amp2, freq, width, peaktime, cutoff, is_h, is_c);
}

void fields_1d::add_continuous_source(component_1d whichf, double freq, double width, double peaktime,
                                      double cutoff, double z) {
  add_source(whichf, freq, width, peaktime, cutoff, z, 1);
}


void fields_1d::step() {
  t += 1;

  phase_material();

  step_h();
  step_h_source(h_sources);

  step_e();
  step_e_source(e_sources);

  step_e_polarization();
  step_polarization_itself();
}

void fields_1d::phase_material() {
  if (new_ma && phasein_time > 0) {
    // First multiply the electric field by epsilon...
    DOCMP {
      for (int i=0;i<(nz+1);i++) ex[cmp][i] /= ma->inveps[i];
    }
    // Then change epsilon...
    ma->mix_with(new_ma, 1.0/phasein_time);
    phasein_time--;
    // Finally divide the electric field by the new epsilon...
    DOCMP {
      for (int i=0;i<(nz+1);i++) ex[cmp][i] *= ma->inveps[i];
    }
  }
}

void fields_1d::step_h() {
  DOCMP {
    int z = 0;
    if (npmlz) for (;z<npmlz;z++) {
      double Czhy = ma->Czhy[z];
      double ooop_Czhy = 1.0/(1.0+0.5*Czhy);
      CM(hy,z) += -c*(CM(ex,z+1)-CM(ex,z))*ooop_Czhy - Czhy*CM(hy,z);
    }
    for (;z<nz-npmlz;z++) {
      CM(hy,z)+= -c*(CM(ex,z+1)-CM(ex,z));
    }
    if (npmlz) for (;z<nz;z++) {
      double Czhy = ma->Czhy[nz-z-1];
      double ooop_Czhy = 1.0/(1.0+0.5*Czhy);
      CM(hy,z) += -c*(CM(ex,z+1)-CM(ex,z))*ooop_Czhy - Czhy*CM(hy,z);
    }
  }
}

void fields_1d::step_e() {
  DOCMP {
    if (k>=0.0) {
      CM(ex,0) += -c*MA(ma->inveps,0)*(CM(hy,0)-EMIKZ(hy,nz-1));
      CM(ex,nz) += -c*MA(ma->inveps,0)*(EIKZ(hy,0)-CM(hy,nz-1));
    }
    int z = 1;
    if (npmlz) for (;z<npmlz+1;z++) {
      double Czex = ma->Czex[z-1];
      double ooop_Czex = 1.0/(1.0+0.5*Czex);
      CM(ex,z) += -c*MA(ma->inveps,z)*(CM(hy,z)-CM(hy,z-1))*ooop_Czex - Czex*CM(ex,z);
    }
    for (;z<nz-npmlz;z++) {
      CM(ex,z) += -c*MA(ma->inveps,z)*(CM(hy,z)-CM(hy,z-1));
    }
    if (npmlz) for (;z<nz;z++) {
      double Czex = ma->Czex[nz-z-1];
      double ooop_Czex = 1.0/(1.0+0.5*Czex);
      CM(ex,z) += -c*MA(ma->inveps,z)*(CM(hy,z)-CM(hy,z-1))*ooop_Czex - Czex*CM(ex,z);
    }
  }
}

void fields_1d::step_h_source(const src_1d *s) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(t);
  if (A == 0.0) {
    step_h_source(s->next);
    return;
  }
  IM(hy,s->z) += imag(A*s->amp);
  RE(hy,s->z) += real(A*s->amp);
  step_h_source(s->next);
}

void fields_1d::step_e_source(const src_1d *s) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(t);
  if (A == 0.0) {
    step_e_source(s->next);
    return;
  }
  IM(ex,s->z) += imag(A*s->amp);
  RE(ex,s->z) += real(A*s->amp);
  step_e_source(s->next);
}

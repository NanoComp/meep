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
    delete[] hr[cmp];
    delete[] hp[cmp];
    delete[] hz[cmp];
    delete[] er[cmp];
    delete[] ep[cmp];
    delete[] ez[cmp];
    delete[] backup_hr[cmp];
    delete[] backup_hp[cmp];
    delete[] backup_hz[cmp];
    delete[] hrp[cmp];
    delete[] hpz[cmp];
    delete[] hzr[cmp];
    delete[] erp[cmp];
    delete[] epz[cmp];
    delete[] ezr[cmp];
    delete[] z_hrp[cmp][0];
    delete[] z_hrp[cmp][1];
    delete[] z_hpz[cmp][0];
    delete[] z_hpz[cmp][1];
    delete[] z_erp[cmp][0];
    delete[] z_erp[cmp][1];
    delete[] z_epz[cmp][0];
    delete[] z_epz[cmp][1];
    delete[] erw[cmp];
    delete[] epw[cmp];
    delete[] ezw[cmp];
    delete[] hrw[cmp];
    delete[] hpw[cmp];
    delete[] hzw[cmp];
  }
  if (nfreq) delete[] freqs;
  if (bands) delete bands;
  delete pol;
  delete olpol;
}

void fields::use_bloch(double tk) {
  if (npmlz) npmlz = 0;
  k = tk;
  cosknz = cos(k*2*pi*inva*nz);
  sinknz = sin(k*2*pi*inva*nz);
  eiknz = complex<double>(cosknz, sinknz);
}

fields::fields(const mat *the_ma, int tm) {
  int r, z;
  ma = new mat(the_ma);
  verbosity = 0;
  nr = ma->nr;
  nz = ma->nz;
  outdir = ma->outdir;
  m = tm;
  phasein_time = 0;
  new_ma = NULL;
  bands = NULL;
  freqs = NULL;
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
    hr[cmp] = new double[(nr+1)*(nz+1)];
    hp[cmp] = new double[(nr+1)*(nz+1)];
    hz[cmp] = new double[(nr+1)*(nz+1)];
    er[cmp] = new double[(nr+1)*(nz+1)];
    ep[cmp] = new double[(nr+1)*(nz+1)];
    ez[cmp] = new double[(nr+1)*(nz+1)];
    backup_hr[cmp] = NULL;
    backup_hp[cmp] = NULL;
    backup_hz[cmp] = NULL;
  }
  if (ez[1] == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  npmlz = ma->npmlz;
  npmlr = ma->npmlr;
  // Note: some of the following sizes are wrong...
  DOCMP {
    hrp[cmp] = new double[npmlr*(nz+1)];
    hpz[cmp] = new double[npmlr*(nz+1)];
    hzr[cmp] = new double[npmlr*(nz+1)];
    erp[cmp] = new double[npmlr*(nz+1)];
    epz[cmp] = new double[npmlr*(nz+1)];
    ezr[cmp] = new double[npmlr*(nz+1)];
    z_hrp[cmp][0] = new double[npmlz*(nr+1)];
    z_hrp[cmp][1] = new double[npmlz*(nr+1)];
    z_hpz[cmp][0] = new double[npmlz*(nr+1)];
    z_hpz[cmp][1] = new double[npmlz*(nr+1)];
    z_erp[cmp][0] = new double[npmlz*(nr+1)];
    z_erp[cmp][1] = new double[npmlz*(nr+1)];
    z_epz[cmp][0] = new double[npmlz*(nr+1)];
    z_epz[cmp][1] = new double[npmlz*(nr+1)];
  }
  DOCMP {
    for (r=0;r<nr+1;r++) for (z=0;z<nz+1;z++)
      CM(hr,r,z) = 0.0;
    for (r=0;r<nr+1;r++) for (z=0;z<nz+1;z++)
      CM(hp,r,z) = 0.0;
    for (r=0;r<nr+1;r++) for (z=0;z<nz+1;z++)
      CM(hz,r,z) = 0.0;
    for (r=0;r<nr+1;r++) for (z=0;z<nz+1;z++)
      CM(er,r,z) = 0.0;
    for (r=0;r<nr+1;r++) for (z=0;z<nz+1;z++)
      CM(ep,r,z) = 0.0;
    for (r=0;r<nr+1;r++) for (z=0;z<nz+1;z++)
      CM(ez,r,z) = 0.0;
    // Now for pml extra fields...
    if (npmlr) {
      for (r=nr-npmlr;r<nr;r++) for (z=0;z<nz+1;z++)
        PMLR(hrp,r,z) = 0.0;
      for (r=nr-npmlr;r<nr;r++) for (z=0;z<nz+1;z++)
        PMLR(hpz,r,z) = 0.0;
      for (r=nr-npmlr;r<nr;r++) for (z=0;z<nz+1;z++)
        PMLR(hzr,r,z) = 0.0;
      for (r=nr-npmlr;r<nr;r++) for (z=0;z<nz+1;z++)
        PMLR(erp,r,z) = 0.0;
      for (r=nr-npmlr;r<nr;r++) for (z=0;z<nz+1;z++)
        PMLR(epz,r,z) = 0.0;
      for (r=nr-npmlr;r<nr;r++) for (z=0;z<nz+1;z++)
        PMLR(ezr,r,z) = 0.0;
    }
    if (npmlz) {
      for (int lr=-1;lr<2;lr+=2) for (int r=0;r<nr+1;r++) for (int iz=0;iz<npmlz;iz++)
        PMLZ(z_hrp,r) = 0.0;
      for (int lr=-1;lr<2;lr+=2) for (int r=0;r<nr+1;r++) for (int iz=0;iz<npmlz;iz++)
        PMLZ(z_hpz,r) = 0.0;
      for (int lr=-1;lr<2;lr+=2) for (int r=0;r<nr+1;r++) for (int iz=0;iz<npmlz;iz++)
        PMLZ(z_erp,r) = 0.0;
      for (int lr=-1;lr<2;lr+=2) for (int r=0;r<nr+1;r++) for (int iz=0;iz<npmlz;iz++)
        PMLZ(z_epz,r) = 0.0;
    }
  }
  nfreq = 0;
  nzflux = 0;
  nrflux = 0;
  setifreqmax_and_iposmax(0, 0);
}

void fields::use_real_sources() {
  if (e_sources) e_sources->use_real_sources();
  if (h_sources) h_sources->use_real_sources();
}

void src::use_real_sources() {
  is_real = 1;
  if (next) next->use_real_sources();
}

void fields::initialize_with_n_te(int ntot) {
  for (int n=0;n<ntot;n++) initialize_with_nth_te(n+1);
}

void fields::initialize_with_n_tm(int ntot) {
  for (int n=0;n<ntot;n++) initialize_with_nth_tm(n+1);
}

#define EIKZR (cmp? coskz : sinkz)
#define IEIKZR (cmp? -sinkz : coskz)
#define EIKZRH (cmp? hcoskz : hsinkz)
#define IEIKZRH (cmp? -hsinkz : hcoskz)

#include <gsl/gsl_sf_bessel.h>

double J(int m, double kr) { return gsl_sf_bessel_Jn(m, kr); }
double Jprime(int m, double kr) { 
  if (m) return 0.5*(J(m-1,kr)-J(m+1,kr));
  else return -J(1,kr);
}
double Jroot(int m, int n) { return gsl_sf_bessel_zero_Jnu(m, n+1); }
double Jmax(int m, int n) {
  double rlow, rhigh = Jroot(m,n), rtry;
  if (n == 0) rlow = 0;
  else rlow = Jroot(m, n-1);
  double jplow = Jprime(m,rlow), jptry;
  do {
    rtry = rlow + (rhigh - rlow)*0.5;
    jptry = Jprime(m,rtry);
    if (jplow*jptry < 0) rhigh = rtry;
    else rlow = rtry;
  } while (rhigh - rlow > rhigh*1e-15);
  return rtry;
}

inline double expi(int cmp, double x) {
  return (cmp) ? cos(x) : sin(x);
}

void fields::initialize_with_nth_te(int np0) {
  const int n = np0 - 0;
  double rmax = Jmax(m,n);
  double ktrans = rmax/nr;
  double kk = k*2*pi*inva;
  double eps = ma->eps[1];
  double inveps = 1.0/eps;
  double om = c*sqrt(inveps*(ktrans*ktrans + kk*kk));
  double omt = om*0.5;
  if (om/(2*pi) > preferred_fmax) preferred_fmax = om/(2*pi);
  for (int r=rmin_bulk(m)-1;r<nr;r++) {
    double Jm = J(m,ktrans*(r+0.5));
    double Jmp_h = (Jm-J(m,ktrans*(r-0.5)))/ktrans;//Jprime(m,ktrans*r);//
    if (!r) Jmp_h = 2.0*Jm/ktrans;
    for (int z=0;z<nz;z++) {
      double kz = -kk*z;
      double kzph = -kk*(z+0.5);
      DOCMP {
        CM(hz,r,z) += Jm*expi(cmp, kz);
        CM(hr,r,z) += (-kk*c/om)*(-c*inveps*om/(ktrans*ktrans))*ktrans*Jmp_h
          *expi(cmp, kzph+pi/2);
        CM(hp,r,z) += (kk*c/om)*(-m*c*inveps/(r+0.5)*om/(ktrans*ktrans))
          *Jm*expi(cmp, kzph);
        
        if (r >= rmin_bulk(m) || m == 1)
          CM(ep,r,z) += (-c*inveps*om/(ktrans*ktrans))*ktrans*Jmp_h
            *expi(cmp, kz-omt+pi/2);
        if (r >= rmin_bulk(m) || m < 2)
          CM(er,r,z) += (-m*c*inveps/(r+0.5)*om/(ktrans*ktrans))*Jm*expi(cmp, kz-omt);
      }
    }
    DOCMP {
      const int z=nz;
      double kz = -kk*z;
      double kzph = -kk*(z+0.5);
      CM(hz,r,z) += Jm*expi(cmp, kz);

      if (r >= rmin_bulk(m) || m == 1)
        CM(ep,r,z) += (-c*inveps*om/(ktrans*ktrans))*ktrans*Jmp_h
          *expi(cmp, kz-omt+pi/2);
      if (r >= rmin_bulk(m) || m < 2)
        CM(er,r,z) += (-m*c*inveps/(r+0.5)*om/(ktrans*ktrans))*Jm*expi(cmp, kz-omt);
    }
  }
}

void fields::initialize_with_nth_tm(int np1) {
  const int n = np1 - 1;
  double rroot = Jroot(m,n);
  double ktrans = rroot/nr;
  double kk = k*2*pi*inva;
  double eps = ma->eps[1];
  double inveps = 1.0/eps;
  double om = c*sqrt(inveps*(ktrans*ktrans + kk*kk));
  double omt = om*0.5;
  if (abs(om/(2*pi)) > preferred_fmax) preferred_fmax = om/(2*pi);
  for (int r=rmin_bulk(m)-1;r<nr;r++) {
    double Jm = J(m,ktrans*r);
    double Jmp = Jprime(m,ktrans*r);
    double Jm_h = J(m,ktrans*(r+0.5));
    double Jmp_h = Jprime(m,ktrans*(r+0.5));
    for (int z=0;z<nz;z++) {
      double kz = kk*z;
      double kzmh = kk*(z-0.5);
      DOCMP {
        if (r >= rmin_bulk(m) || m == 0)
          CM(ez,r,z) += Jm*expi(cmp, kz);
        if (r >= rmin_bulk(m))
          CM(hr,r,z) += (m*c/r*om/(ktrans*ktrans))*Jm*expi(cmp, kz+omt);
        if (r == 0 && m == 1)
          CM(hr,r,z) += (-m*c*om/(ktrans*ktrans))*ktrans*Jmp*expi(cmp, kz+omt);
        CM(hp,r,z) += (c*om/(ktrans*ktrans))*(ktrans*Jmp_h)
                      *expi(cmp, kz+omt+pi/2);
        
        if (r >= rmin_bulk(m))
          CM(ep,r,z) += (-kk*c*inveps/om)*(m*c/r*om/(ktrans*ktrans))*Jm
            *expi(cmp, kzmh);
        CM(er,r,z) += (kk*c*inveps/om)*(c*om/(ktrans*ktrans))*(ktrans*Jmp_h)
                      *expi(cmp, kzmh+pi/2);
      }
    }
    DOCMP {
      const int z=nz;
      double kz = kk*z;
      double kzmh = kk*(z-0.5);
      if (r >= rmin_bulk(m))
        CM(ep,r,z) += (-kk*c*inveps/om)*(m*c/r*om/(ktrans*ktrans))*Jm
          *expi(cmp, kzmh);
      CM(er,r,z) += (kk*c*inveps/om)*(c*om/(ktrans*ktrans))*(ktrans*Jmp_h)
        *expi(cmp, kzmh+pi/2);
    }
  }
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
  double tt = t - peaktime;
  if (fabs(tt) > cutoff) {
    return 0.0;
  }
  double envelope = exp(-tt*tt/(2*width*width));
  if (is_real) return real( (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope );
  else return (polar(1.0,-2*pi*freq*tt) - amp_shift)*envelope;
}

double src::get_envelope_at_time(int t) const {
  double tt = t - peaktime;
  if (fabs(tt) > cutoff) {
    return 0.0;
  }
  return exp(-tt*tt/(2*width*width));
}

static double integrate_envelope(const src *s) {
  if (s == NULL) {
    printf("Bad arg to integrate_envelope!\n");
    exit(1);
  }
  double sofar = 0.0;
  for (int t=0;1<<30;t++) {
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

void fields::add_src_pt(int r, int z,
                        complex<double> Pr, complex<double> Pp, complex<double> Pz,
                        double freq, double width, double peaktime,
                        double cutoff, int is_h) {
  const double pi=3.14159265;
  if (m!=0 && r < rmin_bulk(m)-1) return;
  if (r >= nr - npmlr) return;
  if (Pr == 0 && Pp == 0 && Pz == 0) return;
  if (m == 0 && r == 0 && (Pp != 0 || Pr != 0)) {
    return; // Can't have these sources at r == 0 if m == 0.
  } else if (m ==1 && r == 0 && Pz != 0) {
    return;
  } else if (m > 1 && r < rmin_bulk(m) - 1) {
    return;
  } else if (m > 1 && r == rmin_bulk(m) - 1 &&
             ((is_h && Pr != 0) || (!is_h && Pp != 0) || (!is_h && Pz != 0))) {
    return;
  }
  if (z >= nz || z < 0) {
    printf("Error:  source is outside of cell!\n");
    exit(1);
  }
  src *tmp = new src;
  tmp->freq = freq*c*inva;
  tmp->width = width/tmp->freq; // this is now time width
  tmp->ez = Pz;
  tmp->er = Pr;
  tmp->ep = Pp;
  tmp->r = r;
  tmp->z = z;
  tmp->amp_shift = 0.0;
  tmp->is_real = 0;
  if (is_h) {
    tmp->next = h_sources;
    h_sources = tmp;
  } else {
    tmp->next = e_sources;
    e_sources = tmp;
  }
  tmp->cutoff = 1+ (int)(cutoff*tmp->width);
  tmp->peaktime = peaktime*a;
  if (peaktime <= 0.0) tmp->peaktime = tmp->cutoff;
  // Apply a shift so that we won't end up with a static polarization when
  // the source is gone:
  tmp->amp_shift = integrate_source(tmp)/integrate_envelope(tmp);
}

void fields::find_source_z_position(double z, double shift, int *z1, int *z2,
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

void fields::add_plane_source(double freq, double width, double peaktime,
                              double cutoff, double z, complex<double> amp(double r)) {
  int thez = (int) (z*a + 0.5);
  const complex<double> I = complex<double>(0,1);
  complex<double> eiomt = exp(-I*freq*2*pi*inva/c*0.5);
  if (z == 0.0 && npmlz) thez = npmlz+1;
  for (int r=0;r<nr;r++) {
    double rh = r+0.5, rmh = r-0.5, rp = r+1;
    complex<double> A = amp(r*inva), Ah = amp(rh*inva),
      Ap = amp(rp*inva), Amh = amp(rmh*inva);;

    double eps = sqrt(MA(ma->eps, r, thez));
    if (A != 0.0) {
      // E_phi
      add_src_pt(r, thez, 0.0, A, 0.0, freq, width, peaktime, cutoff, 0);
      // iH_r = d(rH_phi)/dr
      complex<double> A_Hr = -(rh*Ah - rmh*Amh)*eps;
      if (r == 0) A_Hr = Ah*eps;
      add_src_pt(r, thez, -A*eps, 0.0, 0.0, freq, width, peaktime, cutoff, 1);
    }
    if (Ah != 0.0) {
      // iE_r = d(rE_phi)/dr
      complex<double> A_Er = -I*0.5*(Ap + A);
      add_src_pt(r, thez, A_Er, 0.0, 0.0, freq, width, peaktime, cutoff, 0);
      // H_phi
      add_src_pt(r, thez, 0.0, -(Ap+A)*0.5*I*eps, 0.0, freq, width, peaktime, cutoff, 1);
    }
  }
}

void fields::add_hr_source(double freq, double width, double peaktime,
                           double cutoff, double z, complex<double> amp(double r)) {
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, 0.5, &z1, &z2, &amp1, &amp2);
  for (int r=0;r<nr;r++)
    if (amp(r*inva) != 0.0) {
      add_src_pt(r, z1, amp1*amp(r*inva), 0.0 , 0.0, freq, width, peaktime, cutoff, 1);
      add_src_pt(r, z2, amp2*amp(r*inva), 0.0 , 0.0, freq, width, peaktime, cutoff, 1);
    }
}

void fields::add_hp_source(double freq, double width, double peaktime,
                           double cutoff, double z, complex<double> amp(double r)) {
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, 0.5, &z1, &z2, &amp1, &amp2);
  for (int r=0;r<nr;r++)
    if (amp((r+0.5)*inva) != 0.0) {
      add_src_pt(r, z1, 0.0 , amp1*amp((r+0.5)*inva), 0.0, freq, width, peaktime, cutoff, 1);
      add_src_pt(r, z2, 0.0 , amp2*amp((r+0.5)*inva), 0.0, freq, width, peaktime, cutoff, 1);
    }
}

void fields::add_hz_source(double freq, double width, double peaktime,
                           double cutoff, double z, complex<double> amp(double r)) {
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, 0.0, &z1, &z2, &amp1, &amp2);
  for (int r=0;r<nr;r++)
    if (amp((r+0.5)*inva) != 0.0) {
      add_src_pt(r, z1, 0.0 , 0.0, amp1*amp((r+0.5)*inva), freq, width, peaktime, cutoff, 1);
      add_src_pt(r, z2, 0.0 , 0.0, amp2*amp((r+0.5)*inva), freq, width, peaktime, cutoff, 1);
    }
}

void fields::add_er_source(double freq, double width, double peaktime,
                           double cutoff, double z, complex<double> amp(double r)) {
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, 0.0, &z1, &z2, &amp1, &amp2);
  for (int r=0;r<nr;r++)
    if (amp((r+0.5)*inva) != 0.0) {
      add_src_pt(r, z1, amp1*amp((r+0.5)*inva), 0.0 , 0.0, freq, width, peaktime, cutoff);
      add_src_pt(r, z2, amp2*amp((r+0.5)*inva), 0.0 , 0.0, freq, width, peaktime, cutoff);
    }
}

void fields::add_ep_source(double freq, double width, double peaktime,
                           double cutoff, double z, complex<double> amp(double r)) {
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, 0.0, &z1, &z2, &amp1, &amp2);
  for (int r=0;r<nr;r++)
    if (amp(r*inva) != 0.0) {
      add_src_pt(r, z1, 0.0, amp1*amp(r*inva), 0.0, freq, width, peaktime, cutoff);
      add_src_pt(r, z2, 0.0, amp2*amp(r*inva), 0.0, freq, width, peaktime, cutoff);
    }
}

void fields::add_ez_source(double freq, double width, double peaktime,
                           double cutoff, double z, complex<double> amp(double r)) {
  int z1, z2;
  complex<double> amp1, amp2;
  find_source_z_position(z*a, 0.5, &z1, &z2, &amp1, &amp2);
  for (int r=0;r<nr;r++)
    if (amp(r*inva) != 0.0) {
      add_src_pt(r, z1, 0.0 , 0.0, amp1*amp(r*inva), freq, width, peaktime, cutoff);
      add_src_pt(r, z2, 0.0 , 0.0, amp2*amp(r*inva), freq, width, peaktime, cutoff);
    }
}

void fields::step() {
  t += 1;

  phase_material();

  step_h_bulk();
  step_h_pml();
  step_h_boundaries();
  step_h_source(h_sources);

  step_e_bulk();
  step_e_pml();
  step_e_boundaries();
  step_e_source(e_sources);

  step_e_polarization();
  step_e_pml_polarization();
  step_polarization_itself();
  step_polarization_pml();
}

void fields::phase_material() {
  if (new_ma && phasein_time > 0) {
    // First multiply the electric field by epsilon...
    DOCMP {
      for (int i=0;i<(nr+1)*(nz+1);i++) er[cmp][i] /= ma->invepser[i];
      for (int i=0;i<(nr+1)*(nz+1);i++) ep[cmp][i] /= ma->invepsep[i];
      for (int i=0;i<(nr+1)*(nz+1);i++) ez[cmp][i] /= ma->invepsez[i];
      if (npmlr) {
        for (int r=nr-npmlr;r<nr;r++) for (int z=0;z<nz+1;z++)
          PMLR(erp,r,z) /= MA(ma->invepser,r,z);
        for (int r=nr-npmlr;r<nr;r++) for (int z=0;z<nz+1;z++)
          PMLR(epz,r,z) /= MA(ma->invepsep,r,z);
        for (int r=nr-npmlr;r<nr;r++) for (int z=0;z<nz+1;z++)
          PMLR(ezr,r,z) /= MA(ma->invepsez,r,z);
      }
      if (npmlz) {
        for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
          int z0 = npmlz;
          for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
            int z = z0;
            for (int iz=0;iz<npmlz;iz++,z+=lr) {
              PMLZ(z_erp,r) /= MA(ma->invepser,r,z);
              PMLZ(z_epz,r) /= MA(ma->invepsep,r,z);
            }
          }
        }
      }
    }
    // Then change epsilon...
    ma->mix_with(new_ma, 1.0/phasein_time);
    phasein_time--;
    // Finally divide the electric field by the new epsilon...
    DOCMP {
      for (int i=0;i<(nr+1)*(nz+1);i++) er[cmp][i] *= ma->invepser[i];
      for (int i=0;i<(nr+1)*(nz+1);i++) ep[cmp][i] *= ma->invepsep[i];
      for (int i=0;i<(nr+1)*(nz+1);i++) ez[cmp][i] *= ma->invepsez[i];
      if (npmlr) {
        for (int r=nr-npmlr;r<nr;r++) for (int z=0;z<nz+1;z++)
          PMLR(erp,r,z) *= MA(ma->invepser,r,z);
        for (int r=nr-npmlr;r<nr;r++) for (int z=0;z<nz+1;z++)
          PMLR(epz,r,z) *= MA(ma->invepsep,r,z);
        for (int r=nr-npmlr;r<nr;r++) for (int z=0;z<nz+1;z++)
          PMLR(ezr,r,z) *= MA(ma->invepsez,r,z);
      }
      if (npmlz) {
        for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
          int z0 = npmlz;
          for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
            int z = z0;
            for (int iz=0;iz<npmlz;iz++,z+=lr) {
              PMLZ(z_erp,r) *= MA(ma->invepser,r,z);
              PMLZ(z_epz,r) *= MA(ma->invepsep,r,z);
            }
          }
        }
      }
    }
  }
}

void fields::step_h_bulk() {
  DOCMP {
    for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
      double oorph = 1/(r+0.5);
      double oor = 1.0/(double)r;
      double morph = m*oorph;
      double mor = m*oor;
      for (int z=1+npmlz;z<nz-npmlz;z++) {
        CM(hr,r,z)+= c*
          ((CM(ep,r,z+1)-CM(ep,r,z)) - IT(ez,r,z)*mor);
        CM(hp,r,z)+= c*
          ((CM(ez,r+1,z)-CM(ez,r,z)) - (CM(er,r,z+1)-CM(er,r,z)));
        CM(hz,r,z)+= c*
          (IT(er,r,z)*morph - (CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph);
      }
    }
    for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
      const int z=npmlz; // False boundary layer.
      double mor = m/(double)r;
      double oorph = 1/(r+0.5);
      CM(hr,r,z)+= c*
        ((CM(ep,r,z+1)-CM(ep,r,z)) - IT(ez,r,z)*mor);
      CM(hp,r,z)+= c*
        ((CM(ez,r+1,z)-CM(ez,r,z)) - (CM(er,r,z+1)-CM(er,r,z)));
    }
    {
      int r=rmin_bulk(m)-1, z=npmlz;
      double oorph = 1/(r+0.5);
      CM(hp,r,z)+= c* /* false boundary layer */
        ((CM(ez,r+1,z)-CM(ez,r,z)) - (CM(er,r,z+1)-CM(er,r,z)));
    }
    for (int z=1+npmlz;z<nz-npmlz;z++) {
      int r=rmin_bulk(m)-1;
      double oorph = 1.0/(r+0.5);
      double morph = m*oorph;
      CM(hp,r,z)+= c* /* false boundary layer */
        ((CM(ez,r+1,z)-CM(ez,r,z)) - (CM(er,r,z+1)-CM(er,r,z)));
      CM(hz,r,z)+= c* /* false boundary layer */
        (IT(er,r,z)*morph - (CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph);
    }
  }
}

void fields::step_h_pml() {
  DOCMP {
    for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
      double oorph = 1/(r+0.5);
      double oor = 1.0/(double)r;
      double morph = m*oorph;
      double mor = m*oor;
      if (npmlz) {
        { // Do update in low r low z pml region.
          int z0 = npmlz;
          for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
            int z = z0;
            for (int iz=0;iz<npmlz;iz++,z+=lr) {
              CM(hz,r,z)+= c*(IT(er,r,z)*morph -
                              (CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph);
            }
          }
        }
        { // Do update in low r high z pml region.
          int z0 = npmlz-1;
          for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
            int z = z0;
            for (int iz=0;iz<npmlz;iz++,z+=lr) {
              double Czhr = ma->Czhr[iz];
              double Czhp = ma->Czhp[iz];
              double ooop_Czhr = 1.0/(1.0+0.5*Czhr);
              double ooop_Czhp = 1.0/(1.0+0.5*Czhp);
              double dhrp = c*(-IT(ez,r,z)*mor);
              double hrz = CM(hr,r,z) - PMLZ(z_hrp,r);
              PMLZ(z_hrp,r) += dhrp;
              CM(hr,r,z)+= dhrp + ooop_Czhr*(c*(CM(ep,r,z+1)-CM(ep,r,z)) - Czhr*hrz);
              
              double dhpz = ooop_Czhp*(-c*(CM(er,r,z+1)-CM(er,r,z))-Czhp*PMLZ(z_hpz,r));
              PMLZ(z_hpz,r) += dhpz;
              CM(hp,r,z)+= dhpz + c*(CM(ez,r+1,z)-CM(ez,r,z));
            }
          }
        }
      }
    }
    if (npmlz) { // update r minimum and z high and low pml regions.
      int r=rmin_bulk(m)-1;
      double oorph = 1/(r+0.5);
      double morph = m*oorph;
      int z0 = npmlz;
      for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
        int z = z0;
        for (int iz=0;iz<npmlz;iz++,z+=lr) { // False boundary layer!
          CM(hz,r,z)+= c*(IT(er,r,z)*morph -(CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph);
        }
      }
    }
    if (npmlz) { // update r minimum and z high and low pml regions.
      int r=rmin_bulk(m)-1;
      double oorph = 1/(r+0.5);
      double morph = m*oorph;
      int z0 = npmlz-1;
      for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
        int z = z0;
        for (int iz=0;iz<npmlz;iz++,z+=lr) { // False boundary layer!
          double Czhp = ma->Czhp[iz];
          double ooop_Czhp = 1.0/(1.0+0.5*Czhp);
          double dhpz = ooop_Czhp*(-c*(CM(er,r,z+1)-CM(er,r,z))-Czhp*PMLZ(z_hpz,r));
          PMLZ(z_hpz,r) += dhpz;
          CM(hp,r,z)+= dhpz + c*(CM(ez,r+1,z)-CM(ez,r,z));
        }
      }
    }
    if (m==1 && npmlz) {
      int z0 = npmlz-1;
      for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
        int z = z0;
        for (int iz=0;iz<npmlz;iz++,z+=lr) {
          const int r=0;
          double Czhr = ma->Czhr[iz];
          double ooop_Czhr = 1.0/(1.0+0.5*Czhr);
          double dhrp = c*(-IT(ez,r+1,z)*m/* /1.0 */);
          double hrz = CM(hr,r,z) - PMLZ(z_hrp,r);
          PMLZ(z_hrp,r) += dhrp;
          CM(hr,r,z)+= dhrp + ooop_Czhr*(c*(CM(ep,r,z+1)-CM(ep,r,z)) - Czhr*hrz);
        }
      }
    }
    // Now update all the large r pml region (except boundaries)...
    if (npmlr) {
      for (int r=nr-npmlr;r<nr;r++) {
        double oorph = 1/(r+0.5);
        double oor = 1.0/(double)r;
        double morph = m*oorph;
        double mor = m*oor;
        
        double Cphr = ma->Cphr[r-nr+npmlr];
        double Crhp = ma->Crhp[r-nr+npmlr];
        double Crhz = ma->Crhz[r-nr+npmlr];
        double Cphz = ma->Cphz[r-nr+npmlr];
        double ooop_Cphr = 1.0/(1.0+0.5*Cphr);
        double ooop_Crhp = 1.0/(1.0+0.5*Crhp);
        double ooop_Crhz = 1.0/(1.0+0.5*Crhz);
        double ooop_Cphz = 1.0/(1.0+0.5*Cphz);
        for (int z=1;z<nz;z++) {
          double Czhp, Czhr;
          if (z < npmlz) {
            Czhr = ma->Czhr[npmlz - z - 1];
            Czhp = ma->Czhp[npmlz - z - 1];
          } else if (z >= nz - npmlz) {
            Czhr = ma->Czhr[z+npmlz-nz];
            Czhp = ma->Czhp[z+npmlz-nz];
          } else {
            Czhr = 0;
            Czhp = 0;
          }
          double ooop_Czhr = 1.0/(1.0+0.5*Czhr);
          double ooop_Czhp = 1.0/(1.0+0.5*Czhp);

          double dhrp = ooop_Cphr*(- c*IT(ez,r,z)*mor-Cphr*PMLR(hrp,r,z) );
          double hrz = CM(hr,r,z) - PMLR(hrp,r,z);
          PMLR(hrp,r,z) += dhrp;
          CM(hr,r,z)+= dhrp + ooop_Czhr*(c*(CM(ep,r,z+1)-CM(ep,r,z)) - Czhr*hrz);
          double dhpz = ooop_Czhp*(-c*(CM(er,r,z+1)-CM(er,r,z))-Czhp*PMLR(hpz,r,z));
          double hpr = CM(hp,r,z)-PMLR(hpz,r,z);
          PMLR(hpz,r,z) += dhpz;
          CM(hp,r,z)+= dhpz + ooop_Crhp*(c*(CM(ez,r+1,z)-CM(ez,r,z)) - Crhp*hpr);
          double dhzr = ooop_Crhz*(-c*(CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph
                                   - Crhz*PMLR(hzr,r,z));
          double hzp = CM(hz,r,z) - PMLR(hzr,r,z);
          PMLR(hzr,r,z) += dhzr;
          CM(hz,r,z)+= dhzr + ooop_Cphz*(c*IT(er,r,z)*morph - Cphz*hzp);
        }
      }
      for (int r=nr-npmlr;r<nr;r++) {
        const int z=0, iz=npmlz-1, lr=-1; // False boundary layer.
        double Czhr;
        double Czhp;
        if (npmlz) {
          Czhr = ma->Czhr[iz];
          Czhp = ma->Czhp[iz];
        } else {
          Czhr = 0;
          Czhp = 0;
        }
        double mor = m/(double)r;
        double oorph = 1/(r+0.5);
        double Cphr = ma->Cphr[r-nr+npmlr];
        double Crhp = ma->Crhp[r-nr+npmlr];
        double ooop_Cphr = 1.0/(1.0+0.5*Cphr);
        double ooop_Crhp = 1.0/(1.0+0.5*Crhp);
        double ooop_Czhr = 1.0/(1.0+0.5*Czhr);
        double ooop_Czhp = 1.0/(1.0+0.5*Czhp);
        double dhrp = ooop_Cphr*(- c*IT(ez,r,z)*mor-Cphr*PMLR(hrp,r,z) );
        double hrz = CM(hr,r,z) - PMLR(hrp,r,z);
        PMLR(hrp,r,z) += dhrp;
        CM(hr,r,z)+= dhrp + ooop_Czhr*(c*(CM(ep,r,z+1)-CM(ep,r,z)) - Czhr*hrz);
        double dhpz = ooop_Czhp*(-c*(CM(er,r,z+1)-CM(er,r,z))-Czhp*PMLR(hpz,r,z));
        double hpr = CM(hp,r,z)-PMLR(hpz,r,z);
        PMLR(hpz,r,z) += dhpz;
        CM(hp,r,z)+= dhpz + ooop_Crhp*(c*(CM(ez,r+1,z)-CM(ez,r,z)) - Crhp*hpr);
      }
    }
    if (k >= 0.0 && npmlr) { // k < 0 indicates metallic axial boundaries.
      for (int z=0;z<=nz;z+=nz) { // z=0,nz
        for (int r=nr-npmlr;r<nr;r++) {
          double oorph = 1.0/(r+0.5);
          double morph = m*oorph;
          double Crhz = ma->Crhz[r-nr+npmlr];
          double Cphz = ma->Cphz[r-nr+npmlr];
          double ooop_Crhz = 1.0/(1.0+0.5*Crhz);
          double ooop_Cphz = 1.0/(1.0+0.5*Cphz);
          double dhzr = ooop_Crhz*(-c*(CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph
                                   - Crhz*PMLR(hzr,r,z));
          double hzp = CM(hz,r,z) - PMLR(hzr,r,z);
          PMLR(hzr,r,z) += dhzr;
          CM(hz,r,z)+= dhzr + ooop_Cphz*(c*IT(er,r,z)*morph - Cphz*hzp);
        }
      }
    }
  }
}

void fields::step_h_boundaries() {
  DOCMP {
    if (k >= 0.0) { // k < 0 indicates metallic axial boundaries.
      for (int z=0;z<=nz;z+=nz) { // z=0,nz
        for (int r=rmin_bulk(m)-1;r<nr-npmlr;r++) {
          double oorph = 1.0/(r+0.5);
          double morph = m*oorph;
          CM(hz,r,z)+= c*
            (IT(er,r,z)*morph - (CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph);
        }
      }
    }
    if (m==1) {
      for (int z=npmlz;z<nz-npmlz;z++) {
        const int r=0;
        CM(hr,r,z)+= c*
          ((CM(ep,r,z+1)-CM(ep,r,z)) - IT(ez,r+1,z)*m/* /1.0 */);
      }
    }
  }
}

void fields::step_e_bulk() {
  DOCMP {
    for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
      double oorph = 1/(r+0.5);
      double oor = 1.0/(double)r;
      double morph = m*oorph;
      double mor = m*oor;
      for (int z=1+npmlz;z<nz-npmlz;z++) {
        CM(er,r,z)+= c*MA(ma->invepser,r,z)*
          (IT(hz,r,z)*morph - (CM(hp,r,z)-CM(hp,r,z-1)));
        CM(ep,r,z)+= c*MA(ma->invepsep,r,z)*
          ((CM(hr,r,z)-CM(hr,r,z-1)) - (CM(hz,r,z)-CM(hz,r-1,z)));
        CM(ez,r,z)+= c*MA(ma->invepsez,r,z)*
          ((CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor - IT(hr,r,z)*mor);
      }
    }
    if (npmlz==0) {
      for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
        const int z=npmlz;
        double oor = 1.0/(double)r;
        double mor = m*oor;
        CM(ez,r,z)+= c*MA(ma->invepsez,r,z)* /* false boundary layer */
          ((CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor - IT(hr,r,z)*mor);
      }
    }
    for (int z=1+npmlz;z<nz-npmlz;z++) {
      const int r=rmin_bulk(m)-1;
      double oorph = 1/(r+0.5);
      double morph = m*oorph;
      CM(er,r,z)+= c*MA(ma->invepser,r,z)* /* false boundary layer */
        (IT(hz,r,z)*morph - (CM(hp,r,z)-CM(hp,r,z-1)));
    }
  }
}

void fields::step_e_pml() {
  DOCMP {
    for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
      double oorph = 1/(r+0.5);
      double oor = 1.0/(double)r;
      double morph = m*oorph;
      double mor = m*oor;
      if (npmlz) {
        int z0 = npmlz;
        for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
          int z = z0;
          for (int iz=0;iz<npmlz;iz++,z+=lr) {
            double Czer = ma->Czer[iz];
            double Czep = ma->Czep[iz];
            double inveps_plus_Czer = MA(ma->invepser,r,z)/(1+.5*MA(ma->invepser,r,z)*Czer);
            double inveps_plus_Czep = MA(ma->invepsep,r,z)/(1+.5*MA(ma->invepsep,r,z)*Czep);
            double derp = c*MA(ma->invepser,r,z)*(IT(hz,r,z)*morph);
            double erz = CM(er,r,z) - PMLZ(z_erp,r);
            PMLZ(z_erp,r) += derp;
            CM(er,r,z)+= derp + inveps_plus_Czer*(-c*(CM(hp,r,z)-CM(hp,r,z-1)) - Czer*erz);
            double depz = inveps_plus_Czep*(c*(CM(hr,r,z)-CM(hr,r,z-1)) - Czep*PMLZ(z_epz,r));
            double epr = CM(ep,r,z) - PMLZ(z_epz,r);
            PMLZ(z_epz,r) += depz;
            CM(ep,r,z)+= depz + MA(ma->invepsep,r,z)*
              (-c*(CM(hz,r,z)-CM(hz,r-1,z)));
            CM(ez,r,z)+= c*MA(ma->invepsez,r,z)*
              ((CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor - IT(hr,r,z)*mor);
          }
        }
      }
      if (npmlz) {
        const int z=0; // False boundary layer!
        CM(ez,r,z)+= c*MA(ma->invepsez,r,z)*
          ((CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor - IT(hr,r,z)*mor);
      }
    }
    if (npmlz) { // Added Monday Feb 10 -- looks good.
      int z0 = npmlz;
      for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
        int z = z0;
        for (int iz=0;iz<npmlz;iz++,z+=lr) {
          const int r=rmin_bulk(m)-1;
          double oorph = 1/(r+0.5);
          double morph = m*oorph; /* false boundary layer */
          double Czer = ma->Czer[iz];
          double inveps_plus_Czer = MA(ma->invepser,r,z)/(1+.5*MA(ma->invepser,r,z)*Czer);
          double derp = c*MA(ma->invepser,r,z)*(IT(hz,r,z)*morph);
          double erz = CM(er,r,z) - PMLZ(z_erp,r);
          PMLZ(z_erp,r) += derp;
          CM(er,r,z)+= derp + inveps_plus_Czer*(-c*(CM(hp,r,z)-CM(hp,r,z-1)) - Czer*erz);
        }
      }
    }
    if (m==1 && npmlz) {
      int z0 = npmlz;
      for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
        int z = z0;
        for (int iz=0;iz<npmlz;iz++,z+=lr) {
          const int r=0;
          double Czep = ma->Czep[iz];
          double inveps_plus_Czep = MA(ma->invepsep,r,z)/(1+.5*MA(ma->invepsep,r,z)*Czep);
          double depz = inveps_plus_Czep*(c*(CM(hr,r,z)-CM(hr,r,z-1)) - Czep*PMLZ(z_epz,r));
          PMLZ(z_epz,r) += depz;
          CM(ep,r,z)+= depz + c*MA(ma->invepsep,r,z)*(-(CM(hz,r,z)*2.0));
        }
      }
    }
    // update large r pml for all z (except actual boundary)...
    for (int r=nr-npmlr;r<nr;r++) {
      double oorph = 1/(r+0.5);
      double oor = 1.0/(double)r;
      double morph = m*oorph;
      double mor = m*oor;

      double Cper = ma->Cper[r-nr+npmlr];
      double Crep = ma->Crep[r-nr+npmlr];
      double Crez = ma->Crez[r-nr+npmlr];
      double Cpez = ma->Cpez[r-nr+npmlr];
      for (int z=1;z<nz;z++) {
        double Czep, Czer;
        if (z <= npmlz) {
          Czer = ma->Czer[npmlz - z];
          Czep = ma->Czep[npmlz - z];
        } else if (z >= nz - npmlz) {
          Czer = ma->Czer[z+npmlz-nz];
          Czep = ma->Czep[z+npmlz-nz];
        } else {
          Czer = 0;
          Czep = 0;
        }
        double inveps_plus_Cper = MA(ma->invepser,r,z)/(1+.5*MA(ma->invepser,r,z)*Cper);
        double inveps_plus_Czer = MA(ma->invepser,r,z)/(1+.5*MA(ma->invepser,r,z)*Czer);
        double derp = inveps_plus_Cper*(c*IT(hz,r,z)*morph - Cper*PMLR(erp,r,z));
        double erz = CM(er,r,z) - PMLR(erp,r,z);
        PMLR(erp,r,z) += derp;
        CM(er,r,z)+= derp + inveps_plus_Czer*(-c*(CM(hp,r,z)-CM(hp,r,z-1)) - Czer*erz);

        double inveps_plus_Czep = MA(ma->invepsep,r,z)/(1+.5*MA(ma->invepsep,r,z)*Czep);
        double inveps_plus_Crep = MA(ma->invepsep,r,z)/(1+.5*MA(ma->invepsep,r,z)*Crep);
        double depz = inveps_plus_Czep*(c*(CM(hr,r,z)-CM(hr,r,z-1)) - Czep*PMLR(epz,r,z));
        double epr = CM(ep,r,z) - PMLR(epz,r,z);
        PMLR(epz,r,z) += depz;
        CM(ep,r,z)+= depz + inveps_plus_Crep*(-c*(CM(hz,r,z)-CM(hz,r-1,z)) - Crep*epr);

        double inveps_plus_Crez = MA(ma->invepsez,r,z)/(1+.5*MA(ma->invepsez,r,z)*Crez);
        double inveps_plus_Cpez = MA(ma->invepsez,r,z)/(1+.5*MA(ma->invepsez,r,z)*Cpez);
        double dezr = inveps_plus_Crez*(c*(CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor
                                        - Crez*PMLR(ezr,r,z));
        double ezp = CM(ez,r,z)- PMLR(ezr,r,z);
        PMLR(ezr,r,z) += dezr;
        CM(ez,r,z)+= dezr + inveps_plus_Cpez*(-c*IT(hr,r,z)*mor - Cpez*ezp);
      }
    }
    for (int r=nr-npmlr;r<nr;r++) {
      int z = 0; // False boundary layer.
      double oorph = 1/(r+0.5);
      double oor = 1.0/(double)r;
      double morph = m*oorph;
      double mor = m*oor;
      double Crez = ma->Crez[r-nr+npmlr];
      double Cpez = ma->Cpez[r-nr+npmlr];

      double inveps_plus_Crez = MA(ma->invepsez,r,z)/(1+.5*MA(ma->invepsez,r,z)*Crez);
      double inveps_plus_Cpez = MA(ma->invepsez,r,z)/(1+.5*MA(ma->invepsez,r,z)*Cpez);
      double dezr = inveps_plus_Crez*(c*(CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor
                                      - Crez*PMLR(ezr,r,z));
      double ezp = CM(ez,r,z)- PMLR(ezr,r,z);
      PMLR(ezr,r,z) += dezr;
      CM(ez,r,z)+= dezr + inveps_plus_Cpez*(-c*IT(hr,r,z)*mor - Cpez*ezp);
    }
    if (k >= 0.0) { // k < 0 indicates metallic axial boundaries.
      for (int r=nr-npmlr;r<nr;r++) {
        double oor = 1.0/(double)r;
        double mor = m*oor;
        double oorph = 1/(r+0.5);
        double morph = m*oorph;
        double Cper = ma->Cper[r-nr+npmlr];
        double Czer = 0;
        double Czep = 0;
        double Crep = ma->Crep[r-nr+npmlr];
        double inveps_plus_Cper = MA(ma->invepser,r,0)/(1+.5*MA(ma->invepser,r,0)*Cper);
        double inveps_plus_Czer = MA(ma->invepser,r,0)/(1+.5*MA(ma->invepser,r,0)*Czer);
        double inveps_plus_Czep = MA(ma->invepsep,r,0)/(1+.5*MA(ma->invepsep,r,0)*Czep);
        double inveps_plus_Crep = MA(ma->invepsep,r,0)/(1+.5*MA(ma->invepsep,r,0)*Crep);
        {
          int z=0;
          double derp = inveps_plus_Cper*(c*IT(hz,r,z)*morph - Cper*PMLR(erp,r,z));
          double erz = CM(er,r,z) - PMLR(erp,r,z);
          PMLR(erp,r,z) += derp;
          CM(er,r,z)+= derp + inveps_plus_Czer*(-c*(CM(hp,r,z)-EMIKZ(hp,r,nz-1)) - Czer*erz);

          double depz = inveps_plus_Czep*(c*(CM(hr,r,z)-EMIKZ(hr,r,nz-1)) - Czep*PMLR(epz,r,z));
          double epr = CM(ep,r,z) - PMLR(epz,r,z);
          PMLR(epz,r,z) += depz;
          CM(ep,r,z)+= depz + inveps_plus_Crep*(-c*(CM(hz,r,z)-CM(hz,r-1,z)) - Crep*epr);
        }
        { // z=nz
          int z=nz;
          double derp = inveps_plus_Cper*(c*IT(hz,r,z)*morph - Cper*PMLR(erp,r,z));
          double erz = CM(er,r,z) - PMLR(erp,r,z);
          PMLR(erp,r,z) += derp;
          CM(er,r,z)+= derp + inveps_plus_Czep*(-c*(EIKZ(hp,r,0)-CM(hp,r,z-1)) - Czer*erz);
          double depz = inveps_plus_Czep*(c*(EIKZ(hr,r,0)-CM(hr,r,z-1)) - Czep*PMLR(epz,r,z));
          double epr = CM(ep,r,z) - PMLR(epz,r,z);
          PMLR(epz,r,z) += depz;
          CM(ep,r,z)+= depz + inveps_plus_Crep*(-c*(CM(hz,r,z)-CM(hz,r-1,z))-Crep*epr);
        }
      }
    }
  }
}

void fields::step_e_boundaries() {
  DOCMP {
    if (k >= 0.0 && npmlz == 0) { // k < 0 indicates metallic axial boundaries.
      for (int r=rmin_bulk(m);r<nr-npmlr;r++) {
        const int z=0;
        double oor = 1.0/(double)r;
        double mor = m*oor;
        double oorph = 1/(r+0.5);
        double morph = m*oorph;
        CM(er,r,z)+= c*MA(ma->invepser,r,z)*
          (IT(hz,r,z)*morph - (CM(hp,r,z)-EMIKZ(hp,r,nz-1)));
        CM(er,r,nz)+= c*MA(ma->invepser,r,z)*
          (IT(hz,r,nz)*morph - (EIKZ(hp,r,0)-CM(hp,r,nz-1)));

        CM(ep,r, 0)+= c*MA(ma->invepsep,r,z)*
          ((CM(hr,r,0)-EMIKZ(hr,r,nz-1)) - (CM(hz,r, 0)-CM(hz,r-1, 0)));
        CM(ep,r,nz)+= c*MA(ma->invepsep,r,z)*
          (( EIKZ(hr,r,0)-CM(hr,r,nz-1)) - (CM(hz,r,nz)-CM(hz,r-1,nz)));
      }
      {
        const int r=rmin_bulk(m)-1, z=0;
        double oorph = 1/(r+0.5);
        double morph = m*oorph; /* false boundary layer */
        CM(er,r,z)+= c*MA(ma->invepser,r,z)*
          (IT(hz,r,z)*morph - (CM(hp,r,z)-EMIKZ(hp,r,nz-1)));
        CM(er,r,nz)+= c*MA(ma->invepser,r,z)*
          (IT(hz,r,nz)*morph - (EIKZ(hp,r,0)-CM(hp,r,nz-1)));
      }
      if (m==1) { // corner case...
        const int r=0;
        CM(ep,r, 0)+= c*MA(ma->invepsep,r,0)*
          ((CM(hr,r,0)-EMIKZ(hr,r,nz-1)) - (CM(hz,r, 0)*2.0));
        CM(ep,r,nz)+= c*MA(ma->invepsep,r,0)*
          (( EIKZ(hr,r,0)-CM(hr,r,nz-1)) - (CM(hz,r,nz)*2.0));
      }
    }
    if (m==0) {
      for (int z=0;z<nz;z++) { // Works also in PML!  :)
        const int r=0;
        CM(ez,r,z)+= c*MA(ma->invepsez,r,z)*
          (CM(hp,r,z) - IT(hr,r+1,z)*m)/* /1.0 */;
      }
    } else if (m==1) {
      for (int z=1+npmlz;z<nz-npmlz;z++) {
        const int r=0;
        CM(ep,r,z)+= c*MA(ma->invepsep,r,z)*
          ((CM(hr,r,z)-CM(hr,r,z-1)) - (CM(hz,r,z)*2.0));
      }
    }
  }
}

void fields::step_h_source(const src *s) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(t);
  if (A == 0.0) {
    step_h_source(s->next);
    return;
  }
  IM(hr,s->r,s->z) += imag(A*s->er);
  IM(hp,s->r,s->z) += imag(A*s->ep);
  IM(hz,s->r,s->z) += imag(A*s->ez);
  RE(hr,s->r,s->z) += real(A*s->er);
  RE(hp,s->r,s->z) += real(A*s->ep);
  RE(hz,s->r,s->z) += real(A*s->ez);
  if (s->z == 0) {
    IM(hz,s->r,nz) += imag(eiknz*A*s->ez);
    RE(hz,s->r,nz) += real(eiknz*A*s->ez);
  }
  step_h_source(s->next);
}

void fields::step_e_source(const src *s) {
  if (s == NULL) return;
  complex<double> A = s->get_amplitude_at_time(t);
  if (A == 0.0) {
    step_h_source(s->next);
    return;
  }
  IM(er,s->r,s->z) += imag(A*s->er);
  IM(ep,s->r,s->z) += imag(A*s->ep);
  IM(ez,s->r,s->z) += imag(A*s->ez);
  RE(er,s->r,s->z) += real(A*s->er);
  RE(ep,s->r,s->z) += real(A*s->ep);
  RE(ez,s->r,s->z) += real(A*s->ez);
  if (s->z == 0) {
    IM(er,s->r,nz) += imag(eiknz*A*s->er);
    RE(er,s->r,nz) += real(eiknz*A*s->er);
    IM(ep,s->r,nz) += imag(eiknz*A*s->ep);
    RE(ep,s->r,nz) += real(eiknz*A*s->ep);
  }
  step_e_source(s->next);
}

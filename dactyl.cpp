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

inline int small_r_metal(int m) {
  return m-1;
}

inline int rmin_bulk(int m) {
  int r = 1 + small_r_metal(m);
  if (r < 1) r = 1;
  return r;
}

mat::~mat() {
  delete[] invepser;
  delete[] invepsep;
  delete[] invepsez;
  delete[] eps;

  delete[] Crez;
  delete[] Crep;
  delete[] Crhz;
  delete[] Crhp;
  delete[] Cper;
  delete[] Cpez;
  delete[] Cphr;
  delete[] Cphz;

  delete[] Czer;
  delete[] Czep;
  delete[] Czhr;
  delete[] Czhp;
}

static double sig(double r);

void mat::use_pml(int numpmlr, int numpmlz) {
  npmlz = numpmlz;
  npmlr = numpmlr;
  // Delete any previously allocated conductivity arrays...
  delete[] Crez;
  delete[] Crep;
  delete[] Crhz;
  delete[] Crhp;
  delete[] Cper;
  delete[] Cpez;
  delete[] Cphr;
  delete[] Cphz;

  delete[] Czer;
  delete[] Czep;
  delete[] Czhr;
  delete[] Czhp;
  // Allocate the conductivity arrays:
  Crez = new double[npmlr];
  Crep = new double[npmlr];
  Crhz = new double[npmlr];
  Crhp = new double[npmlr];
  Cper = new double[npmlr];
  Cpez = new double[npmlr];
  Cphr = new double[npmlr];
  Cphz = new double[npmlr];

  Czer = new double[npmlz];
  Czep = new double[npmlz];
  Czhr = new double[npmlz];
  Czhp = new double[npmlz];
  // Initialize them... (for now I'm setting them to zero...)
  const double Cmax = 0.5;
  for (int r=0;r<npmlr;r++) {
    double rr = (r)/(double)npmlr;
    double rp = (1+r)/(double)npmlr;
    Crez[r] = Crep[r] = Cmax*0.5*(sig(rp)+sig(rr));
    Crhz[r] = Crhp[r] = Cmax*sig(rp);
  }
  double sigintegrated = 0.0;
  for (int r=0;r<npmlr;r++) {
    double rr = (r)/(double)npmlr;
    double rp = (1+r)/(double)npmlr;
    sigintegrated += Cmax*0.5*(sig(rp)+sig(rr));
    Cper[r] = Cphz[r] = sigintegrated/(nr-npmlr+r+0.5);
    if (r==0) Cpez[r] = Cphr[r] = 0.5*Cper[r];
    else Cpez[r] = Cphr[r] = 0.5*(Cper[r]+Cper[r-1]);

    Cphz[r] = 0.0;
    Cper[r] = 0.0;
    Cphr[r] = 0.0;
    Cpez[r] = 0.0;
    printf("Big sig(Ez)[%d  ] is %10lg, little sig is %10lg\n", r, Crez[r], Cpez[r]);
    printf("Big sig(Hz)[%d  ] is %10lg, little sig is %10lg\n", r, Crhz[r], Cphz[r]);
  }
  for (int z=0;z<npmlz;z++) {
    double rr = (z)/(double)npmlz;
    double rp = (1+z)/(double)npmlz;
    Czer[z] = Czep[z] = Cmax*0.5*(sig(rp)+sig(rr));
    Czhr[z] = Czhp[z] = Cmax*sig(rp);
  }
}

static double sig(double r) {
  return pow(r, 2);
}

mat::mat(double feps(double r, double z),
         double rmax, double zmax, double ta) {
  int r,z;
  outdir = ".";
  numpols = 0;
  a = ta;
  nr = (int) (rmax*a);
  nz = (int) (zmax*a);
  npmlz = 0;
  npmlr = 0;
  if (nz == 0) nz = 1;
  eps = new double[nr*(nz+1)];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (r=0;r<nr;r++) {
    for (z=0;z<nz+1;z++) {
      MA(eps,r,z) = feps((r+0.5)/a,(z+0.5));
    }
  }
  invepser = new double[nr*(nz+1)];
  invepsep = new double[nr*(nz+1)];
  invepsez = new double[nr*(nz+1)];

  // Initialize eps to NaNs so we'll notice if we don't set it properly.
  for (int i=0;i<nr*(nz+1);i++)
    invepser[i] = invepsep[i] = invepsez[i] = sqrt(-1);
  for (r=1;r<nr;r++) {
    for (z=1;z<=nz;z++) {
      MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,z-1));
      MA(invepsep,r,z) = 4./(MA(eps,r,z)+MA(eps,r-1,z)+
                            MA(eps,r,z-1)+MA(eps,r-1,z-1));
      MA(invepsez,r,z) = 2./(MA(eps,r-1,z)+MA(eps,r,z));
    }
  }
  for (r=1;r<nr;r++) {
    {
      const int z=0;
      MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,nz-1));
      MA(invepsep,r,z) = 4./(MA(eps,r,z)+MA(eps,r-1,z)+
                          MA(eps,r,nz-1)+MA(eps,r-1,nz-1));
      MA(invepsez,r,z) = 2./(MA(eps,r-1,z)+MA(eps,r,z));
    }
    {
      const int z=nz;
      MA(invepser,r,z) = 2./(MA(eps,r,0)+MA(eps,r,z-1));
      MA(invepsep,r,z) = 4./(MA(eps,r,0)+MA(eps,r-1,0)+
                          MA(eps,r,z-1)+MA(eps,r-1,z-1));
    }
  }
  for (z=1;z<=nz;z++) {
    {
      const int r=0;
      MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,z-1));
      MA(invepsep,r,z) = 2./(MA(eps,r,z)+ MA(eps,r,z-1));
      MA(invepsez,r,z) = 1./MA(eps,r,z);
    }
  }
  {
    const int r=0,z=0;
    MA(invepser,r,z) = 2./(MA(eps,r,z)+MA(eps,r,nz-1));
    MA(invepsep,r,z) = 2./(MA(eps,r,z)+MA(eps,r,nz-1));
    MA(invepsez,r,z) = 1./MA(eps,r,z);
  }
  {
    const int r=0,z=nz;
    MA(invepser,r,z) = 2./(MA(eps,r,0)+MA(eps,r,nz-1));
    MA(invepsep,r,z) = 2./(MA(eps,r,0)+MA(eps,r,z-1));
    MA(invepsez,r,z) = 1./MA(eps,r,0);
  }
  // Allocate the conductivity arrays:
  Crez = Crep = Crhz = Crhp = Cper = Cpez = Cphr = Cphz = NULL;
  Czer = Czep = Czhr = Czhp = NULL;
}

fields::~fields() {
  DOCMP {
    delete[] hr[cmp];
    delete[] hp[cmp];
    delete[] hz[cmp];
    delete[] er[cmp];
    delete[] ep[cmp];
    delete[] ez[cmp];
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
  for (int i=0;i<numpols;i++) delete[] pol[i];
  if (nfreq) delete[] freqs;
}

void fields::use_bloch(double tk) {
  if (npmlz) npmlz = 0;
  k = tk;
  cosknz = cos(k*inva*nz);
  sinknz = sin(k*inva*nz);
}

fields::fields(const mat *the_ma, int tm) {
  int r, z;
  ma = the_ma;
  nr = ma->nr;
  nz = ma->nz;
  outdir = ma->outdir;
  m = tm;
  k = -1;
  a = ma->a;
  inva = 1.0/a;
  t = 0;
  numpols = ma->numpols;
  h_sources = e_sources = NULL;
  hr[0] = new double[(nr+1)*(nz+1)];
  hp[0] = new double[nr*(nz+1)];
  hz[0] = new double[nr*(nz+1)];
  er[0] = new double[nr*(nz+1)];
  ep[0] = new double[(nr+1)*(nz+1)];
  ez[0] = new double[(nr+1)*(nz+1)];
  hr[1] = new double[(nr+1)*(nz+1)];
  hp[1] = new double[nr*(nz+1)];
  hz[1] = new double[nr*(nz+1)];
  er[1] = new double[nr*(nz+1)];
  ep[1] = new double[(nr+1)*(nz+1)];
  ez[1] = new double[(nr+1)*(nz+1)];
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
    for (r=0;r<nr  ;r++) for (z=0;z<nz+1;z++)
      CM(hp,r,z) = 0.0;
    for (r=0;r<nr  ;r++) for (z=0;z<nz+1;z++)
      CM(hz,r,z) = 0.0;
    for (r=0;r<nr  ;r++) for (z=0;z<nz+1;z++)
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

void fields::add_src_pt(int r, int z,
                        double Pr, double Pp, double Pz,
                        double freq, double width, double peaktime,
                        double cutoff, int is_h) {
  const double pi=3.14159265;
  if (m!=0 && r <= rmin_bulk(m)-1) return;
  if (r >= nr - npmlr) return;
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
}

void fields::add_hr_source(double freq, double width, double peaktime,
                   double cutoff, int z, double amp(double r)) {
  int r;
  for (r=0;r<nr;r++)
    if (amp(r*inva) != 0.0)
      add_src_pt(r, z, amp((r+0.5)*inva), 0.0, 0.0, freq, width, peaktime, cutoff, 1);
}

void fields::add_hp_source(double freq, double width, double peaktime,
                   double cutoff, int z, double amp(double r)) {
  int r;
  for (r=0;r<nr;r++)
    if (amp((r+0.5)*inva) != 0.0)
      add_src_pt(r, z, 0.0 , amp(r*inva), 0.0, freq, width, peaktime, cutoff, 1);
}

void fields::add_hz_source(double freq, double width, double peaktime,
                   double cutoff, int z, double amp(double r)) {
  int r;
  for (r=0;r<nr;r++)
    if (amp((r+0.5)*inva) != 0.0)
      add_src_pt(r, z, 0.0, 0.0, amp(r*inva), freq, width, peaktime, cutoff, 1);
}

void fields::add_er_source(double freq, double width, double peaktime,
                   double cutoff, int z, double amp(double r)) {
  int r;
  for (r=0;r<nr;r++)
    if (amp((r+0.5)*inva) != 0.0)
      add_src_pt(r, z, amp((r+0.5)*inva), 0.0, 0.0, freq, width, peaktime, cutoff);
}

void fields::add_ep_source(double freq, double width, double peaktime,
                   double cutoff, int z, double amp(double r)) {
  int r;
  for (r=0;r<nr;r++)
    if (amp(r*inva) != 0.0)
      add_src_pt(r, z, 0.0 , amp(r*inva), 0.0, freq, width, peaktime, cutoff);
}

void fields::add_ez_source(double freq, double width, double peaktime,
                   double cutoff, int z, double amp(double r)) {
  int r;
  for (r=0;r<nr;r++)
    if (amp(r*inva) != 0.0)
      add_src_pt(r, z, 0.0, 0.0, amp(r*inva), freq, width, peaktime, cutoff);
}

void fields::step() {
  t += 1;

  step_h_bulk();
  step_h_pml();
  step_h_boundaries();
  step_h_source(h_sources);

  step_e_bulk();
  step_e_pml();
  step_e_boundaries();
  step_e_source(e_sources);
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
              double dhrp = c*(-IT(ez,r,z)*mor);
              double hrz = CM(hr,r,z) - PMLZ(z_hrp,r);
              PMLZ(z_hrp,r) += dhrp;
              CM(hr,r,z)+= dhrp + (c*(CM(ep,r,z+1)-CM(ep,r,z))*(1-0.5*Czhr) - Czhr*hrz);
              
              double dhpz = (-c*(CM(er,r,z+1)-CM(er,r,z))*(1-0.5*Czhp)-Czhp*PMLZ(z_hpz,r));
              PMLZ(z_hpz,r) += dhpz;
              CM(hp,r,z)+= dhpz + c*(CM(ez,r+1,z)-CM(ez,r,z));
            }
          }
        }
      }
      /* if (npmlz != 0) { // False boundary layer!
        const int z=0, iz=npmlz, lr=-1;
        double Czhr = ma->Czhr[iz];
        double dhrp = c*(-IT(ez,r,z)*mor);
        double hrz = CM(hr,r,z) - PMLZ(z_hrp,r);
        PMLZ(z_hrp,r) += dhrp;
        CM(hr,r,z)+= dhrp + c*((CM(ep,r,z+1)-CM(ep,r,z))*(1-0.5*Czhr) - Czhr*hrz);

        double Czhp = ma->Czhp[iz];
        double dhpz = c*(-(CM(er,r,z+1)-CM(er,r,z))*(1-0.5*Czhp)-Czhp*PMLZ(z_hpz,r));
        PMLZ(z_hpz,r) += dhpz;
        CM(hp,r,z)+= dhpz + c*(CM(ez,r+1,z)-CM(ez,r,z));
        }*/
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
          double dhpz = (-c*(CM(er,r,z+1)-CM(er,r,z))*(1-0.5*Czhp)-Czhp*PMLZ(z_hpz,r));
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
          double dhrp = c*(-IT(ez,r+1,z)*m/* /1.0 */);
          double hrz = CM(hr,r,z) - PMLZ(z_hrp,r);
          PMLZ(z_hrp,r) += dhrp;
          CM(hr,r,z)+= dhrp +
            (c*(CM(ep,r,z+1)-CM(ep,r,z))*(1-0.5*Czhr) - Czhr*hrz);
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
          double dhrp = (- c*IT(ez,r,z)*mor*(1-0.5*Cphr)-Cphr*PMLR(hrp,r,z) );
          double hrz = CM(hr,r,z) - PMLR(hrp,r,z);
          PMLR(hrp,r,z) += dhrp;
          CM(hr,r,z)+= dhrp + (c*(CM(ep,r,z+1)-CM(ep,r,z))*(1-0.5*Czhr) - Czhr*hrz);
          double dhpz = (-c*(CM(er,r,z+1)-CM(er,r,z))*(1-0.5*Czhp)-Czhp*PMLR(hpz,r,z));
          double hpr = CM(hp,r,z)-PMLR(hpz,r,z);
          PMLR(hpz,r,z) += dhpz;
          CM(hp,r,z)+= dhpz + 
            (c*(CM(ez,r+1,z)-CM(ez,r,z))*(1-0.5*Crhp) - Crhp*hpr);
          double dhzr = (-c*(CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph*(1-0.5*Crhz)
                           - Crhz*PMLR(hzr,r,z));
          double hzp = CM(hz,r,z) - PMLR(hzr,r,z);
          PMLR(hzr,r,z) += dhzr;
          CM(hz,r,z)+= dhzr +
            (c*IT(er,r,z)*morph*(1-0.5*Cphz) - Cphz*hzp);
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
        double dhrp = (- c*IT(ez,r,z)*mor*(1-0.5*Cphr)-Cphr*PMLR(hrp,r,z) );
        double hrz = CM(hr,r,z) - PMLR(hrp,r,z);
        PMLR(hrp,r,z) += dhrp;
        CM(hr,r,z)+= dhrp + (c*(CM(ep,r,z+1)-CM(ep,r,z))*(1-0.5*Czhr) - Czhr*hrz);
        double dhpz = (-c*(CM(er,r,z+1)-CM(er,r,z))*(1-0.5*Czhp)-Czhp*PMLR(hpz,r,z));
        double hpr = CM(hp,r,z)-PMLR(hpz,r,z);
        PMLR(hpz,r,z) += dhpz;
        CM(hp,r,z)+= dhpz + 
          (c*(CM(ez,r+1,z)-CM(ez,r,z))*(1-0.5*Crhp) - Crhp*hpr);
      }
    }
    if (k >= 0.0 && npmlr) { // k < 0 indicates metallic axial boundaries.
      for (int z=0;z<=nz;z+=nz) { // z=0,nz
        for (int r=nr-npmlr;r<nr;r++) {
          double oorph = 1.0/(r+0.5);
          double morph = m*oorph;
          double Crhz = ma->Crhz[r-nr+npmlr];
          double Cphz = ma->Cphz[r-nr+npmlr];
          double dhzr = (-c*(CM(ep,r+1,z)*(r+1.)-CM(ep,r,z)*r)*oorph*(1-0.5*Crhz)
                           - Crhz*PMLR(hzr,r,z));
          double hzp = CM(hz,r,z) - PMLR(hzr,r,z);
          PMLR(hzr,r,z) += dhzr;
          CM(hz,r,z)+= dhzr +
            (c*IT(er,r,z)*morph*(1-0.5*Cphz) - Cphz*hzp);
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
            double derp = c*MA(ma->invepser,r,z)*(IT(hz,r,z)*morph);
            double erz = CM(er,r,z) - PMLZ(z_erp,r);
            PMLZ(z_erp,r) += derp;
            CM(er,r,z)+= derp + MA(ma->invepser,r,z)*
              (-c*(CM(hp,r,z)-CM(hp,r,z-1))*(1-0.5*MA(ma->invepser,r,z)*Czer) - Czer*erz);
            double depz = MA(ma->invepsep,r,z)*(c*(CM(hr,r,z)-CM(hr,r,z-1))*(1-.5*MA(ma->invepsep,r,z)*Czep)
                                               - Czep*PMLZ(z_epz,r));
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
          double Czep = ma->Czep[iz];
          double derp = c*MA(ma->invepser,r,z)*(IT(hz,r,z)*morph);
          double erz = CM(er,r,z) - PMLZ(z_erp,r);
          PMLZ(z_erp,r) += derp;
          CM(er,r,z)+= derp + MA(ma->invepser,r,z)*
            (-c*(CM(hp,r,z)-CM(hp,r,z-1))*(1-0.5*MA(ma->invepser,r,z)*Czer) - Czer*erz);
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
          double depz = MA(ma->invepsep,r,z)*(c*(CM(hr,r,z)-CM(hr,r,z-1))*(1-.5*MA(ma->invepsep,r,z)*Czep)
                                             - Czep*PMLZ(z_epz,r));
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
        double derp = MA(ma->invepser,r,z)*(c*IT(hz,r,z)*morph*(1-.5*MA(ma->invepser,r,z)*Cper)
                                           - Cper*PMLR(erp,r,z));
        double erz = CM(er,r,z) - PMLR(erp,r,z);
        PMLR(erp,r,z) += derp;
        CM(er,r,z)+= derp + MA(ma->invepser,r,z)*
          (-c*(CM(hp,r,z)-CM(hp,r,z-1))*(1-0.5*MA(ma->invepser,r,z)*Czer) - Czer*erz);
        double depz = MA(ma->invepsep,r,z)*(c*(CM(hr,r,z)-CM(hr,r,z-1))*(1-.5*MA(ma->invepsep,r,z)*Czep)
                                           - Czep*PMLR(epz,r,z));
        double epr = CM(ep,r,z) - PMLR(epz,r,z);
        PMLR(epz,r,z) += depz;
        CM(ep,r,z)+= depz + MA(ma->invepsep,r,z)*
          (-c*(CM(hz,r,z)-CM(hz,r-1,z))*(1-.5*MA(ma->invepsep,r,z)*Crep)-Crep*epr);
        double dezr = MA(ma->invepsez,r,z)*
          (c*(CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor*(1-.5*MA(ma->invepsez,r,z)*Crez)
           - Crez*PMLR(ezr,r,z));
        double ezp = CM(ez,r,z)- PMLR(ezr,r,z);
        PMLR(ezr,r,z) += dezr;
        CM(ez,r,z)+= dezr + MA(ma->invepsez,r,z)*
          (-c*IT(hr,r,z)*mor*(1-.5*MA(ma->invepsez,r,z)*Cpez)-Cpez*ezp);
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
      double dezr = MA(ma->invepsez,r,z)*
        (c*(CM(hp,r,z)*(r+0.5)-CM(hp,r-1,z)*(r-0.5))*oor*(1-.5*MA(ma->invepsez,r,z)*Crez)
         - Crez*PMLR(ezr,r,z));
      double ezp = CM(ez,r,z)- PMLR(ezr,r,z);
      PMLR(ezr,r,z) += dezr;
      CM(ez,r,z)+= dezr + MA(ma->invepsez,r,z)*
        (-c*IT(hr,r,z)*mor*(1-.5*MA(ma->invepsez,r,z)*Cpez)-Cpez*ezp);
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
        {
          int z=0;
          double derp = MA(ma->invepser,r,z)*(c*IT(hz,r,z)*morph*(1-.5*MA(ma->invepser,r,z)*Cper)
                                             - Cper*PMLR(erp,r,z));
          double erz = CM(er,r,z) - PMLR(erp,r,z);
          PMLR(erp,r,z) += derp;
          CM(er,r,z)+= derp + MA(ma->invepser,r,z)*
            (-c*(CM(hp,r,z)-EMIKZ(hp,r,nz-1))*(1-0.5*MA(ma->invepser,r,z)*Czer) - Czer*erz);
          double depz = MA(ma->invepsep,r,z)*(c*(CM(hr,r,z)-EMIKZ(hr,r,nz-1))*(1-.5*MA(ma->invepsep,r,z)*Czep)
                                             - Czep*PMLR(epz,r,z));
          double epr = CM(ep,r,z) - PMLR(epz,r,z);
          PMLR(epz,r,z) += depz;
          CM(ep,r,z)+= depz + MA(ma->invepsep,r,z)*
            (-c*(CM(hz,r,z)-CM(hz,r-1,z))*(1-.5*MA(ma->invepsep,r,z)*Crep)-Crep*epr);
        }
        { // z=nz
          int z=nz;
          double derp = MA(ma->invepser,r,z)*(c*IT(hz,r,z)*morph*(1-.5*MA(ma->invepser,r,z)*Cper)
                                             - Cper*PMLR(erp,r,z));
          double erz = CM(er,r,z) - PMLR(erp,r,z);
          PMLR(erp,r,z) += derp;
          CM(er,r,z)+= derp + MA(ma->invepser,r,z)*
            (-c*(EIKZ(hp,r,0)-CM(hp,r,z-1))*(1-0.5*MA(ma->invepser,r,z)*Czer) - Czer*erz);
          double depz = MA(ma->invepsep,r,z)*(c*(EIKZ(hr,r,0)-CM(hr,r,z-1))*(1-.5*MA(ma->invepsep,r,z)*Czep)
                                             - Czep*PMLR(epz,r,z));
          double epr = CM(ep,r,z) - PMLR(epz,r,z);
          PMLR(epz,r,z) += depz;
          CM(ep,r,z)+= depz + MA(ma->invepsep,r,z)*
            (-c*(CM(hz,r,z)-CM(hz,r-1,z))*(1-.5*MA(ma->invepsep,r,z)*Crep)-Crep*epr);
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

        CM(ep,r,0)+= c*MA(ma->invepsep,r,z)*
          ((CM(hr,r,0)-EMIKZ(hr,r,nz-1)) - (CM(hz,r,0)-CM(hz,r-1,0)));
        CM(ep,r,nz)+= c*MA(ma->invepsep,r,z)*
          ((EIKZ(hr,r,0)-CM(hr,r,nz-1)) - (CM(hz,r,nz)-CM(hz,r-1,nz)));
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
        CM(ep,r,0)+= c*MA(ma->invepsep,r,0)*
          ((CM(hr,r,0)-EMIKZ(hr,r,nz-1)) - (CM(hz,r,0)*2.0));
        CM(ep,r,nz)+= c*MA(ma->invepsep,r,0)*
          ((EIKZ(hr,r,0)-CM(hr,r,nz-1)) - (CM(hz,r,nz)*2.0));
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
  double Ar, Ai;
  double tt = t - s->peaktime;
  if (fabs(tt) > s->cutoff) {
    step_h_source(s->next);
    return;
  }
  const double pi = 3.14159265358979323846;
  double envelope = exp(-tt*tt/(2*s->width*s->width));
  Ar = cos(2*pi*s->freq*tt)*envelope;
  if (s->is_real) Ai = 0;
  else Ai = -sin(2*pi*s->freq*tt)*envelope;
  IM(hr,s->r,s->z) += Ai*s->er;
  IM(hp,s->r,s->z) += Ai*s->ep;
  IM(hz,s->r,s->z) += Ai*s->ez;
  RE(hr,s->r,s->z) += Ar*s->er;
  RE(hp,s->r,s->z) += Ar*s->ep;
  RE(hz,s->r,s->z) += Ar*s->ez;
  if (s->z == 0) {
    IM(hz,s->r,nz) += s->ez*(Ai*cosknz + Ar*sinknz);
    RE(hz,s->r,nz) += s->ez*(Ar*cosknz - Ai*sinknz);
  }
  step_h_source(s->next);
}

void fields::step_e_source(const src *s) {
  if (s == NULL) return;
  double Ar, Ai;
  double tt = t - s->peaktime;
  if (fabs(tt) > s->cutoff) {
    Ai = Ar = 0.0;
    step_e_source(s->next);
    return;
  }
  const double pi = 3.14159265358979323846;
  double envelope = exp(-tt*tt/(2*s->width*s->width));
  Ar = cos(2*pi*s->freq*tt)*envelope;
  Ai = -sin(2*pi*s->freq*tt)*envelope;
  IM(er,s->r,s->z) += Ai*s->er;
  IM(ep,s->r,s->z) += Ai*s->ep;
  IM(ez,s->r,s->z) += Ai*s->ez;
  RE(er,s->r,s->z) += Ar*s->er;
  RE(ep,s->r,s->z) += Ar*s->ep;
  RE(ez,s->r,s->z) += Ar*s->ez;
  if (s->z == 0) {
    IM(er,s->r,nz) += s->er*(Ai*cosknz + Ar*sinknz);
    IM(ep,s->r,nz) += s->ep*(Ai*cosknz + Ar*sinknz);
    RE(er,s->r,nz) += s->er*(Ar*cosknz - Ai*sinknz);
    RE(ep,s->r,nz) += s->ep*(Ar*cosknz - Ai*sinknz);
  }
  step_e_source(s->next);
}

/* Below are the output routines. */

double fields::total_energy() {
  double energy = 0, norm = 0, pi=3.14159265;
  DOCMP {
    for (int r=0;r<nr;r++) {
      double rph = r+0.5;
      for (int z=0;z<nz;z++) {
        energy += rph*(1./MA(ma->invepser,r,z))*CM(er,r,z)*CM(er,r,z);
        energy += r*(1./MA(ma->invepsep,r,z))*CM(ep,r,z)*CM(ep,r,z);
        energy += r*(1./MA(ma->invepsez,r,z))*CM(ez,r,z)*CM(ez,r,z);
        energy += r*CM(hr,r,z)*CM(hr,r,z);
        energy += rph*CM(hp,r,z)*CM(hp,r,z);
        energy += rph*CM(hz,r,z)*CM(hz,r,z);
        norm += (r + rph)/2;
      }
    }
  }
  return energy/norm/(8*pi);
}

double fields::zflux(int ri, int ro, int z) {
  double flux = 0;
  double rph, rtemp;
 
  if (ri > ro) SWAP(ri,ro)
  DOCMP
    for (int r=ri;r<=ro;r++) {
      rph   = ((double)r) + 0.5;
      flux += rph*CM(er,r,z)*CM(hp,r,z)-((double)r)*CM(ep,r,z)*CM(hr,r,z);
    }
  return flux;
}

double fields::rflux(int zl, int zu, int r) {
  double flux = 0;
  double rph, rtemp;

  if (zl > zu) SWAP(zl,zu)
  DOCMP
    for (int z=zl;z<=zu;z++)
      flux += CM(ep,r,z)*CM(hz,r,z) - CM(hp,r,z) * CM(ez,r,z);
  return sqrt(((double)r)*(((double)r)+0.5))*flux;
}

void fields::output_point(FILE *o, double r, double z, const char *name) {
  fprintf(o, "%s\t%8lg", name, t*inva);
  fprintf(o, "\t%8lg\t%8lg",
          RE(er,(int)(r*a),(int)(z*a)), IM(er,(int)(r*a),(int)(z*a)));
  fprintf(o, "\t%8lg\t%8lg",
          RE(ep,(int)(r*a),(int)(z*a)), IM(ep,(int)(r*a),(int)(z*a)));
  fprintf(o, "\t%8lg\t%8lg",
          RE(ez,(int)(r*a),(int)(z*a)), IM(ez,(int)(r*a),(int)(z*a)));
  fprintf(o, "\t%8lg\t%8lg",
          RE(hr,(int)(r*a),(int)(z*a)), IM(hr,(int)(r*a),(int)(z*a)));
  fprintf(o, "\t%8lg\t%8lg",
          RE(hp,(int)(r*a),(int)(z*a)), IM(hp,(int)(r*a),(int)(z*a)));
  fprintf(o, "\t%8lg\t%8lg",
          RE(hz,(int)(r*a),(int)(z*a)), IM(hz,(int)(r*a),(int)(z*a)));
  fprintf(o, "\n");
}

static double get_phase(double *f[2], int nr, int nz) {
  complex<double> mean=0;
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz;z++) {
      mean += complex<double>(RE(f,r,z), IM(f,r,z));
    }
  }
  complex<double> meanr=0;
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz;z++) {
      if (RE(f,r,z) > 0) meanr += complex<double>(RE(f,r,z), IM(f,r,z));
      else meanr += complex<double>(-RE(f,r,z), -IM(f,r,z));
    }
  }
  complex<double> meani=0;
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz;z++) {
      if (IM(f,r,z) > 0) meani += complex<double>(RE(f,r,z), IM(f,r,z));
      else meani += complex<double>(-RE(f,r,z), -IM(f,r,z));
    }
  }
  if (abs(mean) > abs(meanr) && abs(mean) > abs(meani)) return arg(mean);
  if (abs(meanr) > abs(meani)) return arg(meanr);
  if (abs(meani) > 0.0) return arg(meani);
  return 0;
}

static void output_complex_slice(double *f[2], int nr, int nz, 
                                 const char *name) {
  int r, z;
  FILE *out = fopen(name, "w");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  double phase = get_phase(f, nr, nz);
  double c = cos(phase), s = sin(phase);
  for (r=0;r<nr;r++) {
    for (z=0;z<nz;z++) {
      fprintf(out, "%d\t%d\t%lg\n", z, r, c*RE(f,r,z)+s*IM(f,r,z));
    }
  }
  fclose(out);
}

static void output_slice(const double *f, int nr, int nz, const char *name) {
  int r, z;
  FILE *out = fopen(name, "w");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  for (r=0;r<nr;r++) {
    for (z=0;z<nz;z++) {
      fprintf(out, "%d\t%d\t%lg\n", z, r, MA(f,r,z));
    }
  }
  fclose(out);
}

void mat::output_sigma_slice(const char *filename) {
  double *sigma = new double[nr*(nz+1)];
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz;z++) {
      MA(sigma,r,z) = 0.0;
    }
  }
  if (npmlz) {
    for (int r=0;r<nr;r++) {
      int z0 = npmlz-1;
      for (int lr=-1;lr<2;lr+=2,z0=nz-npmlz) {
        int z = z0;
        for (int iz=0;iz<npmlz;iz++,z+=lr) {
          MA(sigma,r,z) += Czhp[iz]*Czhp[iz];
        }
      }
    }
  }
  if (npmlr) {
    for (int r=nr-npmlr;r<nr;r++) {
      for (int z=0;z<nz;z++) {
        MA(sigma,r,z) += Cphz[r-nr+npmlr]*Cphz[r-nr+npmlr];
        MA(sigma,r,z) += Crhz[r-nr+npmlr]*Crhz[r-nr+npmlr];
      }
    }
  }
  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz;z++) {
      MA(sigma,r,z) = sqrt(MA(sigma,r,z));
    }
  }
  output_slice(sigma, nr, nz, filename);
  delete[] sigma;
}

void mat::output_slices(const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  if (!n) {
    printf("Allocation failure!\n");
    exit(1);
  }
  snprintf(n, buflen, "%s/%sepsilon.sli", outdir, nname);
  output_slice(eps, nr, nz, n);
  for (int i=0;i<numpols;i++) {
    snprintf(n, buflen, "%s/%sfreq-%d.sli", outdir, nname, i);
    output_slice(freq[i], nr, nz, n);
  }
  snprintf(n, buflen, "%s/%ssigma.sli", outdir, nname);
  output_sigma_slice(n);
  free(n);
}

void fields::output_real_imaginary_slices(const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) {
    printf("Allocation failure!\n");
    exit(1);
  }
  char *r_or_i = "-re";
  for (int cmp=0;cmp<2;cmp++) {
    if (a == 1) snprintf(n, buflen, "%s/%shr%s-%06.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%shr%s-%07.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(hr[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%shp%s-%06.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%shp%s-%07.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(hp[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%shz%s-%06.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%shz%s-%07.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(hz[cmp], nr, nz, n);
    
    if (a == 1) snprintf(n, buflen, "%s/%ser%s-%06.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%ser%s-%07.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(er[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%sep%s-%06.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%sep%s-%07.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(ep[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%sez%s-%06.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%sez%s-%07.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(ez[cmp], nr, nz, n);
    r_or_i = "-im";
  }

  free(n);
}

void fields::output_slices(const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) {
    printf("Allocation failure!\n");
    exit(1);
  }
  for (i=0;i<numpols;i++) {
    if (a == 1) snprintf(n, buflen, "%s/%spol-%d-%06.0f.sli", outdir, nname, i, t*inva);
    else snprintf(n, buflen, "%s/%spol-%d-%07.2f.sli", outdir, nname, i, t*inva);
    output_slice(pol[i], nr, nz, n);
  }

  if (a == 1) snprintf(n, buflen, "%s/%shr-%06.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%shr-%07.2f.sli", outdir, nname, t*inva);
  output_complex_slice(hr, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%shp-%06.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%shp-%07.2f.sli", outdir, nname, t*inva);
  output_complex_slice(hp, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%shz-%06.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%shz-%07.2f.sli", outdir, nname, t*inva);
  output_complex_slice(hz, nr, nz, n);
  
  if (a == 1) snprintf(n, buflen, "%s/%ser-%06.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%ser-%07.2f.sli", outdir, nname, t*inva);
  output_complex_slice(er, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%sep-%06.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%sep-%07.2f.sli", outdir, nname, t*inva);
  output_complex_slice(ep, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%sez-%06.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%sez-%07.2f.sli", outdir, nname, t*inva);
  output_complex_slice(ez, nr, nz, n);

  free(n);
}

int fields::setifreqmax_and_iposmax(int ifreq, int ipos) {
  int j, arraysize;

  if (t==0) {
    if (ifreq >= 0)
      ifreqmax = ifreq;
    else if (ifreqmax < 0) {
      printf("Illegal value of ifreqmax (ifreqmax = %d).\n", ifreqmax);
      return 1;
    }
    if (ipos >= 0)
      iposmax  = ipos;
    else if (iposmax < 0) {
      printf("Illegal value of iposmax (iposmax = %d).\n", iposmax);
      return 1;
    }
    arraysize = ifreqmax*iposmax;
    DOCMP {
      erw[cmp] = new double[arraysize];
      for (j=0; j<arraysize;j++) erw[cmp][j] = 0.0;
      epw[cmp] = new double[arraysize];
      for (j=0; j<arraysize;j++) epw[cmp][j] = 0.0;
      ezw[cmp] = new double[arraysize];
      for (j=0; j<arraysize;j++) ezw[cmp][j] = 0.0;
      hrw[cmp] = new double[arraysize];
      for (j=0; j<arraysize;j++) hrw[cmp][j] = 0.0;
      hpw[cmp] = new double[arraysize];
      for (j=0; j<arraysize;j++) hpw[cmp][j] = 0.0;
      hzw[cmp] = new double[arraysize];
      for (j=0; j<arraysize;j++) hzw[cmp][j] = 0.0;
    }
    return 0;
  } else {
    printf("Can't reset ifreqmax now (on timestep t=%d)!\n", t);
    return 1;
  }
}

int fields::set_frequency_range(double wl, double wu, double deltaw) {
  double wrange, wtemp;

  wrange = wu - wl;
  if (wrange < 0) {
    wtemp = wl;
    wl = wu;
    wu = wtemp;
    wrange = -wrange;
  }
  /*
    if ((int)(wrange / deltaw + 1) > ifreqmax) {
    printf("Warning, number of desired frequencies exceeds ifreqmax.\n");
    deltaw = wrange / ((double) (ifreqmax-1));
    printf("Resetting delta omega to %g.\n", deltaw);
    }
  */
  nfreq = (int) (wrange / deltaw + 1.000001);
  freqs = new double[nfreq];
  for (int i=0;i<nfreq;i++)
    freqs[i] = wl + ((double) i) * deltaw;
  setifreqmax_and_iposmax(nfreq, -1);
  return 0;
}

int fields::add_zfluxplane(int ri, int ro, int z) {
  nzfluxplane[nzflux] = new int[3];
  nzfluxplane[nzflux][0] = ri;
  nzfluxplane[nzflux][1] = ro;
  nzfluxplane[nzflux][2] = z;
  nzflux++;
  setifreqmax_and_iposmax(-1,iposmax+(ro-ri+1));
  return 0;
}

int fields::add_rfluxplane(int zl, int zu, int r) {
  nrfluxplane[nrflux] = new int[3];
  if (zl > zu) SWAP(zl,zu)
  nrfluxplane[nrflux][0] = zl;
  nrfluxplane[nrflux][1] = zu;
  nrfluxplane[nrflux][2] = r;
  nrflux++;
  setifreqmax_and_iposmax(-1,iposmax+(zu-zl+1));
  return 0;
}



void fields::dft_flux() {
  int n, r, ri, ro, z, zl, zu;
  int ipos;
  complex<double> cer, cep, cez, chr, chp, chz;

  ipos = 0;
  for (n=0;n<nzflux;n++) {
    ri = nzfluxplane[n][0];
    ro = nzfluxplane[n][1];
    z  = nzfluxplane[n][2];
    for (r = ri; r <= ro; r++) {
      //may want to average as appropriate over Yee-type lattice
	cer = complex<double>(RE(er,r,z),IM(er,r,z));
	cep = complex<double>(RE(ep,r,z),IM(ep,r,z));
	cez = complex<double>(RE(ez,r,z),IM(ez,r,z));
	chr = complex<double>(RE(hr,r,z),IM(hr,r,z));
	chp = complex<double>(RE(hp,r,z),IM(hp,r,z));
	chz = complex<double>(RE(hz,r,z),IM(hz,r,z));
	ttow(cer, FW(erw,ipos), ((double) t));
	ttow(cep, FW(epw,ipos), ((double) t));
	ttow(cez, FW(ezw,ipos), ((double) t));
	ttow(chr, FW(hrw,ipos), ((double)t)+0.5);
	ttow(chp, FW(hpw,ipos), ((double)t)+0.5);
	ttow(chz, FW(hzw,ipos), ((double)t)+0.5);
	ipos++;
    }
  }
  for (n=0;n<nrflux;n++) {
    zl = nrfluxplane[n][0];
    zu = nrfluxplane[n][1];
    r  = nrfluxplane[n][2];
    for (z = zl; z <= zu; z++) {
      //may want to average as appropriate over Yee-type lattice
      cer = complex<double>(RE(er,r,z),IM(er,r,z));
      cep = complex<double>(RE(ep,r,z),IM(ep,r,z));
      cez = complex<double>(RE(ez,r,z),IM(ez,r,z));
      chr = complex<double>(RE(hr,r,z),IM(hr,r,z));
      chp = complex<double>(RE(hp,r,z),IM(hp,r,z));
      chz = complex<double>(RE(hz,r,z),IM(hz,r,z));
      ttow(cer, FW(erw,ipos), ((double) t));
      ttow(cep, FW(epw,ipos), ((double) t));
      ttow(cez, FW(ezw,ipos), ((double) t));
      ttow(chr, FW(hrw,ipos), ((double)t)+0.5);
      ttow(chp, FW(hpw,ipos), ((double)t)+0.5);
      ttow(chz, FW(hzw,ipos), ((double)t)+0.5);
      ipos++;
    }
  }
  return;
}


void fields::ttow(complex<double> field, double *retarget, double *imtarget, double time) {
  const double twopi=6.283185307;  // = 2 * 3.141592654
  double twopi_c_over_a, phase;
  double *targetptr[2];

  twopi_c_over_a = twopi * c * inva;
  targetptr[0] = retarget;
  targetptr[1] = imtarget;
  for (int i=0;i<nfreq;i++) {
    phase = twopi_c_over_a * freqs[i] * time;
    DOCMP {
      *(targetptr[cmp]) += FIPHI(field,phase);
      targetptr[cmp]++;
    }
  }
}

void fields::fluxw_output(FILE *outpf, char *header) {
  int r, ri, ro, z, zl, zu, ipos, freq;
  double rph, *(sr[nzflux]), *(s[nrflux]);
  int nz, nr;
  int i, j;

  ipos = 0;
  for (i=0;i<nzflux;i++) {
    sr[i] = new double[nfreq];
    for (j=0;j<nfreq;j++) sr[i][j] = 0.0;
  }
  for (i=0;i<nrflux;i++) {
    s[i] = new double[nfreq];
    for (j=0;j<nfreq;j++) s[i][j] = 0.0;
  }
  for (nz=0;nz<nzflux;nz++) {
    ri = nzfluxplane[nz][0];
    ro = nzfluxplane[nz][1];
    z  = nzfluxplane[nz][2];
    for (r=ri;r<=ro;r++) {
      rph = ((double) r) + 0.5;
      for (freq=0;freq<nfreq;freq++) {
        if (abs(FPW(erw,ipos,freq))>1000) printf("large erw at %d.\n",freq);
        if (abs(FPW(hpw,ipos,freq))>1000) printf("large hpw at %d.\n",freq);
        if (abs(FPW(epw,ipos,freq))>1000) printf("large epw at %d.\n",freq);
        if (abs(FPW(hrw,ipos,freq))>1000) printf("large hrw at %d.\n",freq);
	sr[nz][freq] +=rph * real(conj(FPW(erw,ipos,freq))*FPW(hpw,ipos,freq)) 
	  - r * real(conj(FPW(epw,ipos,freq))*FPW(hrw,ipos,freq));
      }
      ipos++;
    }
  }
  for (nr=0;nr<nrflux;nr++) {
    zl = nrfluxplane[nr][0];
    zu = nrfluxplane[nr][1];
    r  = nrfluxplane[nr][2];
    for (z=zl;z<=zu;z++) {
      for (freq=0;freq<nfreq;freq++) {
        s[nr][freq] += real(conj(FPW(epw,ipos,freq))*FPW(hzw,ipos,freq))
	  - real(conj(FPW(hpw,ipos,freq))*FPW(ezw,ipos,freq));
      }
      ipos++;
    }
    for (freq=0;freq<nfreq;freq++)
      s[nr][freq] *= sqrt( ((double) r) * ((double) (r+0.5)) );
  }
  for (freq=0;freq<nfreq;freq++) {
    fprintf(outpf, "%s: %10.6g, ", header, freqs[freq]);
    for (nz=0;nz<nzflux;nz++)
      fprintf(outpf, "%10.6g ", sr[nz][freq]);
    for (nr=0;nr<nrflux;nr++)
      fprintf(outpf, "%10.6g ", s[nr][freq]);
    fprintf(outpf, "\n");
  }
  for (nr=0;nr<nrflux;nr++)
    delete[] s[nr];
  for (nz=0;nz<nzflux;nz++)
    delete[] sr[nz];
  return;
}

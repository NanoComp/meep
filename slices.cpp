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
    if (a == 1) snprintf(n, buflen, "%s/%shr%s-%08.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%shr%s-%09.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(hr[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%shp%s-%08.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%shp%s-%09.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(hp[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%shz%s-%08.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%shz%s-%09.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(hz[cmp], nr, nz, n);
    
    if (a == 1) snprintf(n, buflen, "%s/%ser%s-%08.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%ser%s-%09.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(er[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%sep%s-%08.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%sep%s-%09.2f.sli", outdir, nname, r_or_i, t*inva);
    output_slice(ep[cmp], nr, nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%sez%s-%08.0f.sli", outdir, nname, r_or_i, t*inva);
    else snprintf(n, buflen, "%s/%sez%s-%09.2f.sli", outdir, nname, r_or_i, t*inva);
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

  if (a == 1) snprintf(n, buflen, "%s/%shr-%08.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%shr-%09.2f.sli", outdir, nname, t*inva);
  output_complex_slice(hr, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%shp-%08.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%shp-%09.2f.sli", outdir, nname, t*inva);
  output_complex_slice(hp, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%shz-%08.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%shz-%09.2f.sli", outdir, nname, t*inva);
  output_complex_slice(hz, nr, nz, n);
  
  if (a == 1) snprintf(n, buflen, "%s/%ser-%08.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%ser-%09.2f.sli", outdir, nname, t*inva);
  output_complex_slice(er, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%sep-%08.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%sep-%09.2f.sli", outdir, nname, t*inva);
  output_complex_slice(ep, nr, nz, n);
  if (a == 1) snprintf(n, buflen, "%s/%sez-%08.0f.sli", outdir, nname, t*inva);
  else snprintf(n, buflen, "%s/%sez-%09.2f.sli", outdir, nname, t*inva);
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

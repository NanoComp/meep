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

#include "tidod.h"
#include "tidod_internals.h"
#include "harminv.h"

#define BAND(b,z,t) ((b)[(z)+(t)*nz])

bandsdata_1d::bandsdata_1d() {
  verbosity = 0;
  maxbands = -1;
  tstart = nz = 0;
  tend = -1;
  hy = ex = NULL;
}

bandsdata_1d::~bandsdata_1d() {
  delete[] hy;
  delete[] ex;
}

int src_1d::find_last_source(int sofar) {
  if (peaktime + cutoff > sofar) sofar = (int)peaktime + cutoff;
  if (next == NULL) return sofar;
  return next->find_last_source(sofar);
}

int fields_1d::find_last_source() {
  int last_source = 0;
  if (e_sources != NULL)
    last_source = e_sources->find_last_source();
  if (h_sources != NULL)
    last_source = max(last_source, h_sources->find_last_source());
  return last_source;  
}

void fields_1d::prepare_for_bands(double endtime, double fmax,
                                  double qmin, double frac_pow_min) {
  int last_source = find_last_source();
  last_source = max(last_source, t + phasein_time);
  if (fmax == 0) fmax = preferred_fmax;
  else preferred_fmax = fmax;
  if (!bands) bands = new bandsdata_1d;
  bands->tstart = last_source+1;
  if (bands->tstart < t) bands->tstart = t;
  bands->tend = t + (int)(endtime*a/c) - 1;
  bands->fpmin = frac_pow_min;
  bands->nz = nz;

  // Set fmin properly...
  double epsmax = 1;
  for (int z=0;z<nz;z++) {
    if (MA(ma->eps,z) > epsmax) epsmax = MA(ma->eps,z);
  }
  bands->fmin = k*c*sqrt(1.0/epsmax);
  bands->qmin = qmin;
  // Set fmax and determine how many timesteps to skip over...
  bands->fmax = fmax;
  {
    // for when there are too many data points...
    double decayconst = bands->fmax*(c*inva)/qmin*8.0;
    double smalltime = 1./(decayconst + bands->fmax*(c*inva));
    bands->scale_factor = (int)(0.06*smalltime);
    if (bands->scale_factor < 1) bands->scale_factor = 1;
    if (verbosity) printf("scale_factor is %d (%lg,%lg)\n",
                          bands->scale_factor, bands->fmax*(c*inva), decayconst);
  }

  if (bands->tend <= bands->tstart) {
    printf("Oi, we don't have any time to take a fourier transform!\n");
    printf("FT start is %d and end is %d\n", bands->tstart, bands->tend);
    exit(1);
  }
  bands->ntime = (1+(bands->tend-bands->tstart)/bands->scale_factor);
  bands->a = a;
  bands->inva = inva;
  bands->hy = new cmplx[nz*bands->ntime];
  bands->ex = new cmplx[nz*bands->ntime];
  if (bands->ex == NULL) {
    printf("Unable to allocate bandstructure array!\n");
    exit(1);
  }
  bands->verbosity = verbosity;
}

void fields_1d::record_bands() {
  if (t > bands->tend || t < bands->tstart) return;
  if (t % bands->scale_factor != 0) return;
  int thet = (t-bands->tstart)/bands->scale_factor;
  if (thet >= bands->ntime) return;
  for (int z=0;z<nz;z++) {
    BAND(bands->hy,z,thet) = cmplx(RE(hy,z), IM(hy,z));
    BAND(bands->ex,z,thet) = cmplx(RE(ex,z), IM(ex,z));
  }
}

complex<double> fields_1d::get_band(int nn, int maxbands) {
  complex<double> *fad = clever_cluster_bands(maxbands);
  complex<double> thef = fad[nn-1];
  delete[] fad;
  return thef;
}

void fields_1d::grace_bands(grace *g, int maxbands) {
  double *approx_power = new double[maxbands];
  //complex<double> *fad = get_the_bands(maxbands, approx_power);
  complex<double> *fad = clever_cluster_bands(maxbands, approx_power);

  int num_found = 0;
  for (int i=0;i<maxbands;i++) if (fad[i] != 0) num_found = i+1;

  for (int i = 0; i < num_found; ++i) {
    g->output_out_of_order(i, k, fabs(real(fad[i])), fabs(imag(fad[i])),
                           approx_power[i]);
  }
  delete[] fad;
  delete[] approx_power;
}

static inline int get_closest(double f[], int fmax) {
  double deltamin = 1e300;
  for (int i=0;i<fmax-1;i++) {
    if (f[i+1]-f[i] < deltamin) deltamin = f[i+1]-f[i];
  }
  for (int i=0;i<fmax-1;i++) {
    if (f[i+1]-f[i] == deltamin) return i;
  }
  return 0;
}

static inline int go_higher(double f[], int fmax, int lo, int hi) {
  if (lo == 0) return 1;
  if (hi == fmax-1) return 0;
  if (f[lo]-f[lo-1] > f[hi+1]-f[hi]) return 1;
  else return 0;
}

static inline int am_done(double f[], int fmax, int lo, int hi) {
  double wid = f[hi]-f[lo] + 0.001;
  return f[lo]-f[lo-1] > wid && f[hi+1]-f[hi] > wid;
}

static void get_cluster(double f[], int fmax, int maxsize, double maxwid,
                        int *out_lo, int *out_hi) {
  int lo = get_closest(f,fmax);
  int hi = lo+1;
  int minsize = maxsize/2+1;
  if (minsize < 3) minsize = 3;
  for (int i=0;i<minsize-2;i++) {
    if (go_higher(f,fmax,lo,hi)) hi++;
    else lo--;
  }
  while (hi + 1 - lo < maxsize) {
    if (am_done(f,fmax,lo,hi)) break;
    if (go_higher(f,fmax,lo,hi)) {
      if (f[hi+1]-f[lo] > maxwid) break;
      hi++;
    } else {
      if (f[hi]-f[lo-1] > maxwid) break;
      lo--;
    }
  }
  *out_lo = lo;
  *out_hi = hi;
}

int fields_1d::cluster_some_bands_cleverly(double *tf, double *td, complex<double> *ta,
                                        int num_freqs, int fields_considered, int maxbands,
                                        complex<double> *fad, double *approx_power) {
  const double total_time = (bands->tend-bands->tstart)*c/a;
  const double deltaf = 1.0/total_time;
  int freqs_so_far = num_freqs;
  printf("About to sort by frequency... (%d frequencies)\n", freqs_so_far);
  // Sort by frequency...
  for (int i = 1; i < freqs_so_far; i++) {
    for (int j=i; j>0;j--) {
      if (tf[j]<tf[j-1]) {
        double t1 = tf[j], t2 = td[j];
        tf[j] = tf[j-1];
        td[j] = td[j-1];
        tf[j-1] = t1;
        td[j-1] = t2;
        complex<double> temp = ta[j];
        ta[j] = ta[j-1];
        ta[j-1] = temp;
      }
    }
  }
  printf("Looking for clusters...\n");
  int num_found = 0;
  double totwid = 0.001;
  while (freqs_so_far >= fields_considered/2 + 1) {
    int hi, lo;
    get_cluster(tf,freqs_so_far,fields_considered,deltaf,&lo,&hi);
    int mid = lo + (hi-lo)/2;
    if (tf[hi]-tf[lo] < deltaf) {
      printf("Got a cluster from %lg to %lg (%d freqs)\n", tf[lo], tf[hi], 1+hi-lo);
      fad[num_found] = complex<double>(tf[mid],td[mid]);
      if (approx_power) {
        approx_power[num_found] = 0;
        for (int i=lo;i<=hi;i++) {
          if (abs(ta[i])*abs(ta[i]) > approx_power[num_found]) {
            approx_power[num_found] = abs(ta[i])*abs(ta[i]);
          }
        }
      }
      totwid += tf[hi]-tf[lo];
      num_found++;
      if (num_found >= maxbands) num_found--;
    } else {
      printf("Rejected a cluster from %lg to %lg (%d freqs)\n", tf[lo], tf[hi], 1+hi-lo);
      if (verbosity > 1) printf("width is %g vs %g\n", tf[hi] - tf[lo], deltaf);
      lo = get_closest(tf,freqs_so_far);
      hi = lo+1;
      if (verbosity > 1) printf("dropping %g and %g\n", tf[hi], tf[lo]);
    }
    freqs_so_far -= 1 + hi - lo;
    for (int i=lo;i<freqs_so_far;i++) {
      tf[i] = tf[i + 1 + hi - lo];
      td[i] = td[i + 1 + hi - lo];
      ta[i] = ta[i + 1 + hi - lo];
    }
  }
  for (int i=0;i<freqs_so_far;i++) {
    if (verbosity > 1) printf("Have a leftover freq: %g\n", tf[i]);
  }
  return num_found;
}

complex<double> *fields_1d::clever_cluster_bands(int maxbands, double *approx_power) {
  const double total_time = (bands->tend-bands->tstart)*c/a;
  const double deltaf = 1.0/total_time;
  bands->maxbands = maxbands;
  const int max_harminvs = 120;
  const int max_freqs = max_harminvs*maxbands;
  double *tf = new double[max_freqs];
  double *td = new double[max_freqs];
  cmplx *ta = new cmplx[max_freqs];
  const int ntime = bands->ntime;
  cmplx *simple_data = new cmplx[ntime];
  if (!simple_data || !ta) {
    printf("Error allocating...\n");
    exit(1);
  }
  /*
  freqs_here = bands->get_freqs(bands->P, ntime,
                               ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
  if (freqs_here) {
    fields_considered++;
    freqs_so_far += freqs_here;
  }
  */
  complex<double> *fad = new complex<double>[maxbands];
  for (int i=0;i<maxbands;i++) fad[i] = 0.0;
  int freqs_so_far = 0;
  int fields_considered = 0;

  for (int z=1;z<nz;z+=1+(int)(bands->scale_factor/c*3.99)) {
    cmplx *bdata;
    for (int whichf = 0; whichf < 1; whichf++) {
      switch (whichf) {
      case 0: bdata = bands->ex; break;
      case 1: bdata = bands->hy; break;
      }
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bdata,z,t);
      int freqs_here = bands->get_freqs(simple_data, ntime,
                                        ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
      if (freqs_so_far + maxbands > max_freqs) break;
    }
  }
  int num_found = cluster_some_bands_cleverly(tf, td, ta, freqs_so_far, fields_considered,
                                              maxbands, fad, approx_power);
    
  delete[] ta;
  delete[] tf;
  delete[] td;
  // Sorting by frequency again...
  for (int i=0;i<num_found;i++) {
    for (int j=i; j>0;j--) {
      if (real(fad[j])<real(fad[j-1])) {
        complex<double> t1 = fad[j];
        fad[j] = fad[j-1];
        fad[j-1] = t1;
        double temp = approx_power[j];
        approx_power[j] = approx_power[j-1];
        approx_power[j-1] = temp;
      }
    }
  }
  return fad;
}

int bandsdata_1d::get_freqs(cmplx *data, int n,
                         cmplx *amps, double *freq_re, double *freq_im) {
  
  int num = do_harminv(data, n, scale_factor, a, fmin, fmax, maxbands, amps, 
		       freq_re, freq_im);
  // First deal with any negative frequency solutions.
  const double total_time = n*scale_factor*c/a;
  for (int i=0;i<num-1;i++) {
    if (freq_re[i] < 0) {
      for (int j=i+1;j<num;j++) {
        if (abs(freq_re[j]+freq_re[i]) < 2.0/total_time) {
          if (verbosity > 2 && freq_re[i] != 0.0) {
            printf("Got a plus/minus freq match at %lg\n",freq_re[j]); 
            printf("Total time: %lg and delta freq limit %lg\n",
                   total_time, 2.0/total_time);
          }
          freq_re[i] = -0.0; // It will get cleaned up later...
        }
      }
      freq_re[i] = -freq_re[i];
      if (verbosity > 2 && freq_re[i] != 0.0)
        printf("Flipping sign of a negative freq:  %lg %lg\n", freq_re[i], freq_im[i]);
    }
  }
  // Now sort the silly solutions again...
  for (int i=0;i<num-1;i++) { // This is a really bad sort algorithm...
    for (int j=i+1;j<num;j++) {
      if (freq_re[i] > freq_re[j]) {
        double t = freq_re[i];
        freq_re[i] = freq_re[j];
        freq_re[j] = t;
        t = freq_im[i];
        freq_im[i] = freq_im[j];
        freq_im[j] = t;
        complex<double> tc = amps[i];
        amps[i] = amps[j];
        amps[j] = tc;
      }
    }
  }
  // Now get rid of any spurious low frequency solutions...
  int orignum = num;
  for (int i=0;i<orignum;i++) {
    if (freq_re[0] < fmin*.9) {
      if (verbosity > 2 && freq_re[0] != 0.0) {
        printf("Trashing a spurious low frequency solution with freq %lg %lg\n",
               freq_re[0], freq_im[0]);
        //printf("For your info, fmin is %lg\n", fmin);
      }
      for (int j=0;j<num-1;j++) {
        freq_re[j]=freq_re[j+1];
        freq_im[j]=freq_im[j+1];
        amps[j]=amps[j+1];
      }
      num--;
    }
  }
  // Now get rid of any spurious transient solutions...
  for (int i=num-1;i>=0;i--) {
    double qminhere = 1.0/(1.0/qmin + 0.25/(freq_re[i]*total_time));
    double qhere = 0.5*fabs(fabs(freq_re[i])/freq_im[i]);
    if (qhere < qminhere) {
      num--;
      if (verbosity > 2) {
        printf("Trashing a spurious low Q solution with freq %lg %lg (%lg vs %lg)\n",
               freq_re[i], freq_im[i], qhere, qminhere);
      }
      for (int j=i;j<num;j++) {
        freq_re[j] = freq_re[j+1];
        freq_im[j] = freq_im[j+1];
        amps[j] = amps[j+1];
      }
    }
  }
  return num;
}

int do_harminv(cmplx *data, int n, int sampling_rate, double a, 
	       double fmin, double fmax, int maxbands,
	       cmplx *amps, double *freq_re, double *freq_im, double *errors) {
  // data is a size n array.

  // check for all zeros in input
  {
    int all_zeros = 1;
    for (int i=0; i<n; i++)
      if (data[i] != 0) all_zeros = 0;
    if (all_zeros)
      return 0;
  }

  harminv_data hd = 
    harminv_data_create(n, data, fmin*sampling_rate*c/a, fmax*sampling_rate*c/a, maxbands);

  int prev_nf, cur_nf;
  harminv_solve(hd);
  prev_nf = cur_nf = harminv_get_num_freqs(hd);

  /* keep re-solving as long as spurious solutions are eliminated */
  do {
    prev_nf = cur_nf;
    harminv_solve_again(hd);
    cur_nf = harminv_get_num_freqs(hd);
  } while (cur_nf < prev_nf);
  if (cur_nf > prev_nf)
    fprintf(stderr,
            "harminv: warning, number of solutions increased from %d to %d!\n",
            prev_nf, cur_nf);
  
  cmplx *tmpamps = harminv_compute_amplitudes(hd);
  double *tmperrors = harminv_compute_frequency_errors(hd);

  freq_re[0] = a*harminv_get_freq(hd, 0)/c/sampling_rate;
  freq_im[0] = -1/(2*pi)*a*harminv_get_decay(hd, 0)/c/sampling_rate;
  for (int i = 1; i < harminv_get_num_freqs(hd); ++i) {
    freq_re[i] = a*harminv_get_freq(hd, i)/c/sampling_rate;
    freq_im[i] = -1/(2*pi)*a*harminv_get_decay(hd, i)/c/sampling_rate;
    for (int j=i; j>0;j--) {
      if (freq_re[j]<freq_re[j-1]) {
        double t1 = freq_re[j], t2 = freq_im[j], e = tmperrors[j];
        cmplx a = tmpamps[j];
        tmpamps[j] = tmpamps[j-1];
        tmperrors[j] = tmperrors[j-1];
        freq_re[j] = freq_re[j-1];
        freq_im[j] = freq_im[j-1];
        freq_re[j-1] = t1;
        freq_im[j-1] = t2;
        tmpamps[j-1] = a;
        tmperrors[j-1] = e;
      }
    }
  }
  int num = harminv_get_num_freqs(hd);
  for (int i = 0; i < num; ++i) {
    amps[i] = tmpamps[i];
    if (errors)
      errors[i] = tmperrors[i];
  }
  free(tmpamps);
  free(tmperrors);
  harminv_data_destroy(hd);
  return num;
}



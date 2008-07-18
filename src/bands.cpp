/* Copyright (C) 2005-2008 Massachusetts Institute of Technology  
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

#include "config.h"
#ifdef HAVE_HARMINV
#  include <harminv.h>
#endif

namespace meep {

#define BAND(b,r,t) ((b)[(r)+(t)*nr])

bandsdata::bandsdata() {
  verbosity = 0;
  maxbands = -1;
  tstart = 0;
  for (int i=0;i<num_bandpts;i++) index[i] = -1;
  tend = -1;
  for (int i=0;i<num_bandpts;i++) for (int c=0;c<10;c++) f[i][c] = NULL;
  P = NULL;
}

bandsdata::~bandsdata() {
  for (int i=0;i<num_bandpts;i++) for (int c=0;c<10;c++) delete[] f[i][c];
  delete[] P;
}

double fields::last_source_time() {
  double last_time = 0;
    if (sources != NULL)
      last_time = max(last_time, sources->last_time_max());
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      last_time = max(last_time,chunks[i]->last_source_time());
  return max_to_all(last_time);
}

double fields_chunk::last_source_time() {
  return 0;
}

void fields::prepare_for_bands(const vec &p, double endtime, double fmax,
                               double qmin, double frac_pow_min) {
  int last_source = (int)(last_source_time()/dt+0.5);
  last_source = max(last_source, t + phasein_time);
  if (!bands) bands = new bandsdata;
  bands->tstart = last_source+1;
  if (bands->tstart < t) bands->tstart = t;
  bands->tend = t + (int)(endtime/dt) - 1;

  {
    int ind[8];
    double w[8];
    int indind = 0;
    while (bands->index[indind] != -1 && indind < num_bandpts) indind++;
    for (int h=0;h<num_chunks;h++)
      if (chunks[h]->v.contains(p)) {
        chunks[h]->v.interpolate(chunks[h]->v.eps_component(), p, ind, w);
        for (int i=0;i<8&&w[i]&&indind<num_bandpts;i++) {
          bands->chunk[indind] = chunks[h];
          bands->index[indind++] = ind[i];
        }
        break;
      }
  }
  bands->fpmin = frac_pow_min;

  // Set fmin properly...
  const double epsmax = max_eps();

  double cutoff_freq = 0.0;
  if (v.dim == Dcyl) {
    cutoff_freq = 1.84*a*dt/(2*pi)/v.nr()/sqrt(epsmax);
    if (m == 0) cutoff_freq *= 0.5;
  }
  bands->fmin = sqrt(cutoff_freq*cutoff_freq + abs(k[Z])*abs(k[Z])*(a*dt)*(a*dt)/epsmax); // FIXME
  bands->fmin = cutoff_freq/dt;
  bands->qmin = qmin;
  // Set fmax and determine how many timesteps to skip over...
  bands->fmax = fmax;
  {
    // for when there are too many data points...
    double decayconst = bands->fmax*dt/qmin*8.0;
    double smalltime = 1./(decayconst + bands->fmax*dt);
    bands->scale_factor = (int)(0.06*smalltime);
    if (bands->scale_factor < 1) bands->scale_factor = 1;
    if (verbosity) master_printf("scale_factor is %d (%g,%g)\n",
                                 bands->scale_factor, bands->fmax*dt,
                                 decayconst);
  }

  if (bands->tend <= bands->tstart) {
    printf("Oi, we don't have any time to take a fourier transform!\n");
    abort("FT start is %d and end is %d\n", bands->tstart, bands->tend);
  }
  bands->ntime = (1+(bands->tend-bands->tstart)/bands->scale_factor);
  bands->dt = dt * bands->scale_factor;
  for (int c=0;c<10;c++)
    for (int i=0;i<num_bandpts;i++)
      if (v.has_field((component)c)) {
        delete[] bands->f[i][c];
        bands->f[i][c] = new complex<double>[bands->ntime];
        if (bands->f[i][c] == NULL)
          abort("Unable to allocate bandstructure array!\n");
        for (int j=0;j<bands->ntime;j++) bands->f[i][c][j] = 0.0;
      }
  bands->P = new complex<double>[bands->ntime];
  for (int i=0;i<bands->ntime;i++) bands->P[i] = 0.0;
  bands->verbosity = verbosity;
  for (int h=0;h<num_chunks;h++)
    chunks[h]->bands = bands;
}

void fields::record_bands() {
  if (!bands) return;
  if (t > bands->tend || t < bands->tstart) return;
  if (t % bands->scale_factor != 0) return;
  for (int i=0;i<num_chunks;i++)
    chunks[i]->record_bands(t);
}

void fields_chunk::record_bands(int tcount) {
  int thet = (tcount-bands->tstart)/bands->scale_factor;
  if (thet >= bands->ntime) return;
  for (int p=0; p<num_bandpts && bands->index[p]!=-1; p++)
    if (this == bands->chunk[p])
      for (int c=0;c<10;c++)
        if (v.has_field((component)c)) {
          complex<double> tmp;
          if (f[c][0] && f[c][1]) tmp = complex<double>(f[c][0][bands->index[p]],
                                                        f[c][1][bands->index[p]]);
          bands->f[p][c][thet] = broadcast(n_proc(), tmp);
        }
}

#define HARMOUT(o,n,f) ((o)[(n)+(f)*maxbands])

complex<double> fields::get_band(int nn, int maxbands) {
  //complex<double> *fad = get_the_bands(maxbands, approx_power);
  complex<double> *fad = clever_cluster_bands(maxbands);
  complex<double> thef = fad[nn-1];
  delete[] fad;
  return thef;
}

void fields::grace_bands(grace *g, int maxbands) {
  double *approx_power = new double[maxbands];
  //complex<double> *fad = get_the_bands(maxbands, approx_power);
  complex<double> *fad = clever_cluster_bands(maxbands, approx_power);

  int num_found = 0;
  for (int i=0;i<maxbands;i++) if (fad[i] != 0.0) num_found = i+1;
  for (int i = 0; i < num_found; ++i) {
    g->output_out_of_order(i, abs(k[Z]), fabs(real(fad[i])), fabs(imag(fad[i])),
                           approx_power[i]); // FIXME
  }
  delete[] fad;
  delete[] approx_power;
}

void fields::output_bands(FILE *o, const char *name, int maxbands) {
  out_bands(o, name, maxbands);
}

void fields::out_bands(FILE *o, const char *name, int maxbands) {
  double *approx_power = new double[maxbands];
  //complex<double> *fad = get_the_bands(maxbands, approx_power);
  complex<double> *fad = clever_cluster_bands(maxbands, approx_power);

  complex<double> *eigen = new complex<double>[maxbands*6];
  if (!eigen) abort("Error allocating...\n");

  for (int whichf = 0; whichf < 6; whichf++) {
    for (int n=0;n<maxbands;n++) {
      HARMOUT(eigen,n,whichf) = 0;
    }
  }
  int num_found = 0;
  for (int i=0;i<maxbands;i++) if (fad[i] != 0.0) num_found = i+1;

  if (am_master()) {
    for (int i = 0; i < num_found; ++i) {
      // k k k m index freq decay Q approx_power
      fprintf(o, "%s\t%g\t%g\t%g\t%g\t%d\t%g \t%g \t%g \t%g\n", 
	      name, 
	      real(k[0]), real(k[1]), real(k[2]),   
	      m, i, fabs(real(fad[i])), imag(fad[i]),
	      -fabs(real(fad[i])) / (2 * imag(fad[i])),
	      approx_power[i]);
    }
    fflush(o);
  }
  delete[] approx_power;
  delete[] fad;
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
  int lodone = lo == 0 || f[lo]-f[lo-1] > wid;
  int hidone = hi == fmax-1 || f[hi+1]-f[hi] > wid;
  return lodone && hidone;
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

int fields::cluster_some_bands_cleverly(double *tf, double *td, complex<double> *ta,
                                        int num_freqs, int fields_considered,
                                        int maxbands,
                                        complex<double> *fad, double *approx_power) {
  const double total_time = (bands->tend-bands->tstart)*dt;
  const double deltaf = 1.0/total_time;
  int freqs_so_far = num_freqs;
  if (!quiet) master_printf("About to sort by frequency... (%d frequencies)\n", freqs_so_far);
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
  if (!quiet) master_printf("Looking for clusters...\n");
  int num_found = 0;
  double totwid = 0.001;
  while (freqs_so_far >= fields_considered/2 + 1) {
    int hi, lo;
    get_cluster(tf,freqs_so_far,fields_considered,deltaf,&lo,&hi);
    int mid = lo + (hi-lo)/2;
    if (tf[hi]-tf[lo] < deltaf) {
      if (!quiet) master_printf("Got a cluster from %g to %g (%d freqs)\n",
				tf[lo], tf[hi], 1+hi-lo);
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
      if (!quiet)
	master_printf("Rejected a cluster from %g to %g (%d/%d freqs)\n",
		      tf[lo], tf[hi], 1+hi-lo, fields_considered);
      if (verbosity > 1) master_printf("width is %g vs %g\n",
                                       tf[hi] - tf[lo], deltaf);
      lo = get_closest(tf,freqs_so_far);
      hi = lo+1;
      if (verbosity > 1) master_printf("dropping %g and %g\n", tf[hi], tf[lo]);
    }
    freqs_so_far -= 1 + hi - lo;
    for (int i=lo;i<freqs_so_far;i++) {
      tf[i] = tf[i + 1 + hi - lo];
      td[i] = td[i + 1 + hi - lo];
      ta[i] = ta[i + 1 + hi - lo];
    }
  }
  for (int i=0;i<freqs_so_far;i++) {
    if (verbosity > 1) master_printf("Have a leftover freq: %g\n", tf[i]);
  }
  return num_found;
}

complex<double> *fields::clever_cluster_bands(int maxbands, double *approx_power) {
  bands->maxbands = maxbands;
  const int max_harminvs = 120;
  const int max_freqs = max_harminvs*maxbands;
  double *tf = new double[max_freqs];
  double *td = new double[max_freqs];
  complex<double> *ta = new complex<double>[max_freqs];
  const int ntime = bands->ntime;
  if (!ta) abort("Error allocating...\n");
  int num_found = 0;
  complex<double> *fad = new complex<double>[maxbands];
  for (int i=0;i<maxbands;i++) fad[i] = 0.0;
  int freqs_so_far = 0;
  int fields_considered = 0;

  for (int p=0; p<num_bandpts && bands->index[p]!=-1; p++)
    for (int whichf = 0; whichf < 10; whichf++)
      if (v.has_field((component)whichf) && maxbands < max_freqs - freqs_so_far) {
        if (verbosity>1) master_printf("Looking at field %d\n", whichf);
        int freqs_here = bands->get_freqs(bands->f[p][whichf], ntime,
                                          ta+freqs_so_far,
                                          tf+freqs_so_far,
                                          td+freqs_so_far);
        if (freqs_here) {
          fields_considered++;
          freqs_so_far += freqs_here;
        }
        if (freqs_so_far + maxbands > max_freqs) break;
      }
  if (k == 0 && v.dim == Dcyl && m != 0) fields_considered /= 2;
  num_found = cluster_some_bands_cleverly(tf, td, ta, freqs_so_far, fields_considered,
                                          maxbands, fad, approx_power);
  delete[] ta;
  delete[] tf;
  delete[] td;
  // Get rid of bands with too little power in them...
  {
    double maxp = 0.0;
    for (int i=0;i<num_found;i++)
      maxp = max(maxp, approx_power[i]);
    double minp = maxp*bands->fpmin;
    for (int i=0;i<num_found;i++)
      if (approx_power[i] < minp) {
        for (int j=i; j<num_found-1;j++) {
          fad[j] = fad[j+1];
          approx_power[j] = approx_power[j+1];
        }
        num_found--;
        i--;
        fad[num_found] = 0.0;
      }
  }
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

int bandsdata::get_freqs(complex<double> *data, int n, complex<double> *amps,
                         double *freq_re, double *freq_im) {
  
  const double total_time = n*dt;
  const double qminhere = 1.0/(1.0/qmin + 0.25/(fmin*total_time));
  return do_harminv(data, n, dt, fmin, fmax, maxbands,
		    amps,  freq_re, freq_im, NULL, 1.1, qminhere);
}

int do_harminv(complex<double> *data, int n, double dt, 
	       double fmin, double fmax, int maxbands,
	       complex<double> *amps, double *freq_re, double *freq_im, double *errors,
	       double spectral_density, double Q_thresh, double rel_err_thresh, double err_thresh, double rel_amp_thresh, double amp_thresh) {
#ifndef HAVE_HARMINV
  abort("compiled without Harminv library, required for do_harminv");
  return 0;
#else
  int numfreqs = int(fabs(fmax-fmin)*dt*n*spectral_density); // c.f. harminv
  if (numfreqs > 150) numfreqs = 150; // prevent matrices from getting too big
  if (numfreqs < 2) numfreqs = 2;

  if (maxbands > numfreqs)
       numfreqs = maxbands;

  // check for all zeros in input
  // data is a size n array.
  {
    int i;
    for (i = 0; i < n && data[i] == 0.0; i++)
      ;
    if (i == n)
      return 0;
  }

#if 0
  // debugging: save data file and arguments for standalone harminv program
  {
    FILE *f = fopen("harminv.dat", "w");
    fprintf(f, "# -f %d -t %g %g-%g -Q %e -e %e -E %e -a %e -A %e -F\n",
	    numfreqs, dt, fmin, fmax,
	    Q_thresh, rel_err_thresh, err_thresh, rel_amp_thresh, amp_thresh);
    for (int i = 0; i < n; ++i)
      fprintf(f, "%g%+gi\n", real(data[i]), imag(data[i]));
    fclose(f);
  }
#endif

  harminv_data hd = harminv_data_create(n, data, fmin*dt, fmax*dt, numfreqs);
  harminv_solve(hd);

  int nf = harminv_get_num_freqs(hd);
  if (nf == 0) return 0;
  int *fsort = new int[nf]; // indices of frequencies, sorted as needed
  
  for (int i = 0; i < nf; ++i)
    fsort[i] = i;
  for (int i = 0; i < nf; ++i) // sort in increasing order of error
    for (int j = i + 1; j < nf; ++j) 
      if (harminv_get_freq_error(hd, fsort[i]) >
	  harminv_get_freq_error(hd, fsort[j])) {
	int k = fsort[i];
	fsort[i] = fsort[j];
	fsort[j] = k;
      }
  
  double min_err = harminv_get_freq_error(hd, fsort[0]);
  double max_amp = abs(harminv_get_amplitude(hd, 0));
  for (int i = 1; i < nf; ++i) {
    double amp = abs(harminv_get_amplitude(hd, i));
    if (max_amp < amp)
      max_amp = amp;
  }
  { // eliminate modes that fall outside the various thresholds:
    int j = 0;
    for (int i = 0; i < nf; ++i) {
      double f = abs(harminv_get_freq(hd, fsort[i]) / dt);
      double err = harminv_get_freq_error(hd, fsort[i]);
      double amp = abs(harminv_get_amplitude(hd, fsort[i]));
      if (f >= fmin && f <= fmax
	  && abs(harminv_get_Q(hd, fsort[i])) > Q_thresh
	  && err < err_thresh
	  && err < rel_err_thresh * min_err
	  && amp > amp_thresh
	  && amp > rel_amp_thresh * max_amp) {
	fsort[j++] = fsort[i];
      }
    }
    nf = j;
  }
  { // eliminate positive/negative frequency pairs
    // set indices to -1 for frequencies to be eliminated
    for (int i = 0; i < nf; ++i) 
      if (fsort[i] != -1) { // i hasn't been eliminated yet
	double f = harminv_get_freq(hd, fsort[i]);
	if (f < 0.0) {
	  double kdiff = -2 * f;
	  int kpos = i;
	  for (int k = 0; k < nf; ++k) // search for closest positive freq.
	    if (fsort[k] != -1) { // k hasn't been eliminated yet
	      double fdiff = abs(harminv_get_freq(hd, fsort[k]) + f);
	      if (fdiff < kdiff) {
		kpos = k;
		kdiff = fdiff;
	      }
	    }
	  if (kpos != i && kdiff < 2.0 / n) { // consider them the same
	    // pick the one with the smaller error
	    if (harminv_get_freq_error(hd, fsort[i]) <
		harminv_get_freq_error(hd, fsort[kpos]))
	      fsort[kpos] = -1;
	    else
	      fsort[i] = -1;
	  }
	}
      }
    int j = 0;
    for (int i = 0; i < nf; ++i) // remove the eliminated indices
      if (fsort[i] != -1)
	fsort[j++] = fsort[i];
    nf = j;
  }
  
  if (nf > maxbands)
    nf = maxbands;
  
  // sort again, this time in increasing order of freq:
  for (int i = 0; i < nf; ++i) // simple O(nf^2) sort
    for (int j = i + 1; j < nf; ++j) 
      if (abs(harminv_get_freq(hd, fsort[i])) >
	  abs(harminv_get_freq(hd, fsort[j]))) {
	int k = fsort[i];
	fsort[i] = fsort[j];
	fsort[j] = k;
      }
  
  for (int i = 0; i < nf; ++i) {
    complex<double> freq = harminv_get_omega(hd, fsort[i]) / (2*pi*dt);
    freq_re[i] = abs(real(freq));
    freq_im[i] = imag(freq);
    amps[i] = harminv_get_amplitude(hd, fsort[i]);
    if (errors)
      errors[i] = harminv_get_freq_error(hd, fsort[i]);
  }

  delete[] fsort;
  harminv_data_destroy(hd);
  return nf;
#endif
}

} // namespace meep

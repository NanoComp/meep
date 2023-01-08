/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "meep.hpp"
#include "meep_internals.hpp"

#include "config.h"
#ifdef HAVE_HARMINV
#include <harminv.h>
#endif

using namespace std;

namespace meep {

double fields::last_source_time() {
  double last_time = 0;
  if (sources != NULL) last_time = std::max(last_time, sources->last_time_max());
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) last_time = std::max(last_time, chunks[i]->last_source_time());
  return max_to_all(last_time);
}

double fields_chunk::last_source_time() { return 0; }

/* backwards compatibility with harminv < 1.4 */
#if HARMINV_VERSION_MAJOR < 1 || (HARMINV_VERSION_MAJOR == 1 && HARMINV_VERSION_MINOR < 4)
#define harminv_get_amplitude(pa, d, k) *(pa) = harminv_get_amplitude(d, k)
#define harminv_get_omega(pw, d, k) *(pw) = harminv_get_omega(d, k)
#endif

int do_harminv(complex<double> *data, int n, double dt, double fmin, double fmax, int maxbands,
               complex<double> *amps, double *freq_re, double *freq_im, double *errors,
               double spectral_density, double Q_thresh, double rel_err_thresh, double err_thresh,
               double rel_amp_thresh, double amp_thresh) {
#ifndef HAVE_HARMINV
  (void)data;
  (void)n;
  (void)dt;
  (void)fmin;
  (void)fmax;
  (void)maxbands;
  (void)amps;
  (void)freq_re;
  (void)freq_im;
  (void)errors;
  (void)spectral_density;
  (void)Q_thresh;
  (void)rel_err_thresh;
  (void)err_thresh;
  (void)rel_amp_thresh;
  (void)amp_thresh;
  meep::abort("compiled without Harminv library, required for do_harminv");
  return 0;
#else
  int numfreqs = int(fabs(fmax - fmin) * dt * n * spectral_density); // c.f. harminv
  if (numfreqs > 150) numfreqs = 150; // prevent matrices from getting too big
  if (numfreqs < 2) numfreqs = 2;

  if (maxbands > numfreqs) numfreqs = maxbands;

  // check for all zeros in input
  // data is a size n array.
  {
    int i;
    for (i = 0; i < n && data[i] == 0.0; i++)
      ;
    if (i == n) return 0;
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

  harminv_data hd = harminv_data_create(n, data, fmin * dt, fmax * dt, numfreqs);
  harminv_solve(hd);

  int nf = harminv_get_num_freqs(hd);
  if (nf == 0) return 0;
  int *fsort = new int[nf]; // indices of frequencies, sorted as needed

  for (int i = 0; i < nf; ++i)
    fsort[i] = i;
  for (int i = 0; i < nf; ++i) // sort in increasing order of error
    for (int j = i + 1; j < nf; ++j)
      if (harminv_get_freq_error(hd, fsort[i]) > harminv_get_freq_error(hd, fsort[j])) {
        int k = fsort[i];
        fsort[i] = fsort[j];
        fsort[j] = k;
      }

  double min_err = harminv_get_freq_error(hd, fsort[0]);
  complex<double> aa;
  harminv_get_amplitude(&aa, hd, 0);
  double max_amp = abs(aa);
  for (int i = 1; i < nf; ++i) {
    harminv_get_amplitude(&aa, hd, i);
    double amp = abs(aa);
    if (max_amp < amp) max_amp = amp;
  }
  { // eliminate modes that fall outside the various thresholds:
    int j = 0;
    for (int i = 0; i < nf; ++i) {
      double f = abs(harminv_get_freq(hd, fsort[i]) / dt);
      double err = harminv_get_freq_error(hd, fsort[i]);
      harminv_get_amplitude(&aa, hd, fsort[i]);
      double amp = abs(aa);
      if (f >= fmin && f <= fmax && abs(harminv_get_Q(hd, fsort[i])) > Q_thresh &&
          err < err_thresh && err < rel_err_thresh * min_err && amp > amp_thresh &&
          amp > rel_amp_thresh * max_amp) {
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
            if (fsort[k] != -1) {      // k hasn't been eliminated yet
              double fdiff = abs(harminv_get_freq(hd, fsort[k]) + f);
              if (fdiff < kdiff) {
                kpos = k;
                kdiff = fdiff;
              }
            }
          if (kpos != i && kdiff < 2.0 / n) { // consider them the same
            // pick the one with the smaller error
            if (harminv_get_freq_error(hd, fsort[i]) < harminv_get_freq_error(hd, fsort[kpos]))
              fsort[kpos] = -1;
            else
              fsort[i] = -1;
          }
        }
      }
    int j = 0;
    for (int i = 0; i < nf; ++i) // remove the eliminated indices
      if (fsort[i] != -1) fsort[j++] = fsort[i];
    nf = j;
  }

  if (nf > maxbands) nf = maxbands;

  // sort again, this time in increasing order of freq:
  for (int i = 0; i < nf; ++i) // simple O(nf^2) sort
    for (int j = i + 1; j < nf; ++j)
      if (abs(harminv_get_freq(hd, fsort[i])) > abs(harminv_get_freq(hd, fsort[j]))) {
        int k = fsort[i];
        fsort[i] = fsort[j];
        fsort[j] = k;
      }

  for (int i = 0; i < nf; ++i) {
    complex<double> freq;
    harminv_get_omega(&freq, hd, fsort[i]);
    freq /= (2 * pi * dt);
    freq_re[i] = abs(real(freq));
    freq_im[i] = imag(freq);
    harminv_get_amplitude(&(amps[i]), hd, fsort[i]);
    if (errors) errors[i] = harminv_get_freq_error(hd, fsort[i]);
  }

  delete[] fsort;
  harminv_data_destroy(hd);
  return nf;
#endif
}

} // namespace meep

/* Copyright (C) 2005-2012 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

dft_ldos::dft_ldos(double freq_min, double freq_max, int Nfreq)
{
  if (Nfreq <= 1) {
    omega_min = (freq_min + freq_max) * pi;
    domega = 0;
    Nomega = 1;
  }
  else {
    omega_min = freq_min * 2*pi;
    domega = (freq_max - freq_min) * 2*pi / Nfreq;
    Nomega = Nfreq;
  }
  Fdft = new complex<realnum>[Nomega];
  Jdft = new complex<realnum>[Nomega];
}

double *dft_ldos::ldos() const {
  double *sum = new double[Nomega];
  for (int i = 0; i < Nomega; ++i)
    sum[i] = -real(Fdft[i] * conj(Jdft[i]));
  double *out = new double[Nomega];
  sum_to_all(sum, out, Nomega);
  delete[] sum;
  return out;
}

complex<double> *dft_ldos::F() const {
  complex<double> *out = new complex<double>[Nomega];
  sum_to_all(Fdft, out, Nomega);
  return out;
}

complex<double> *dft_ldos::J() const {
  complex<double> *out = new complex<double>[Nomega];
  sum_to_all(Jdft, out, Nomega);
  return out;
}

void dft_ldos::update(fields &f)
{
  complex<realnum> EJ = 0.0; // integral E * J*
  complex<realnum> HJ = 0.0; // integral H * J* for magnetic currents

  double scale = (f.dt/sqrt(2*pi));

  for (int ic=0;ic<f.num_chunks;ic++) if (f.chunks[ic]->is_mine()) {
      for (src_vol *sv = f.chunks[ic]->sources[D_stuff]; sv; sv = sv->next) {
	component c = direction_component(Ex, component_direction(sv->c));
	realnum *fr = f.chunks[ic]->f[c][0];
	realnum *fi = f.chunks[ic]->f[c][1];
	if (fr && fi) // complex E
	  for (int j=0; j<sv->npts; j++) {
	    const int idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    EJ += complex<realnum>(fr[idx],fi[idx]) * conj(A);
	  }
	else if (fr) { // E is purely real
	  for (int j=0; j<sv->npts; j++) {
	    const int idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    EJ += fr[idx] * conj(A);
	  }
	}
      }
      for (src_vol *sv = f.chunks[ic]->sources[B_stuff]; sv; sv = sv->next) {
	component c = direction_component(Hx, component_direction(sv->c));
	realnum *fr = f.chunks[ic]->f[c][0];
	realnum *fi = f.chunks[ic]->f[c][1];
	if (fr && fi) // complex H
	  for (int j=0; j<sv->npts; j++) {
	    const int idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    HJ += complex<realnum>(fr[idx],fi[idx]) * conj(A);
	  }
	else if (fr) { // H is purely real
	  for (int j=0; j<sv->npts; j++) {
	    const int idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    HJ += fr[idx] * conj(A);
	  }
	}
      }
    }
  for (int i = 0; i < Nomega; ++i) {
    complex<realnum> Ephase = polar(1.0, (omega_min+i*domega)*f.time())*scale;
    complex<realnum> Hphase = polar(1.0, (omega_min+i*domega)*(f.time()-f.dt/2))*scale;
    Fdft[i] += Ephase * EJ + Hphase * HJ;

    // NOTE: take only 1st time dependence: assumes all sources have same J(t)
    if (f.sources)
      Jdft[i] += Ephase * f.sources->current();
  }
}

}

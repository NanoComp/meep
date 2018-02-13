/* Copyright (C) 2005-2015 Massachusetts Institute of Technology.
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

using namespace std;

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
    domega = (freq_max - freq_min) * 2*pi / (Nfreq - 1);
    Nomega = Nfreq;
  }
  Fdft = new complex<realnum>[Nomega];
  Jdft = new complex<realnum>[Nomega];
  for (int i = 0; i < Nomega; ++i) Fdft[i] = Jdft[i] = 0.0;
  Jsum = 1.0;
}

// |c|^2
static double abs2(complex<double> c) {return real(c)*real(c)+imag(c)*imag(c);}

double *dft_ldos::ldos() const {
  // we try to get the overall scale factor right (at least for a point source)
  // so that we can compare against the analytical formula for testing
  // ... in most practical cases, the scale factor won't matter because
  //     the user will compute the relative LDOS of 2 cases (e.g. LDOS/vacuum)

  // overall scale factor
  double Jsum_all = sum_to_all(Jsum);
  double scale = 4.0/pi // from definition of LDOS comparison to power
    * -0.5 // power = -1/2 Re[E* J]
    / (Jsum_all * Jsum_all); // normalize to unit-integral current

  double *sum = new double[Nomega];
  for (int i = 0; i < Nomega; ++i) /* 4/pi * work done by unit dipole */
    sum[i] = scale * real(Fdft[i] * conj(Jdft[i])) / abs2(Jdft[i]);
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
  complex<double> EJ = 0.0; // integral E * J*
  complex<double> HJ = 0.0; // integral H * J* for magnetic currents

  double scale = (f.dt/sqrt(2*pi));

  // compute Jsum for LDOS normalization purposes
  // ...don't worry about the tiny inefficiency of recomputing this repeatedly
  Jsum = 0.0;

  for (int ic=0;ic<f.num_chunks;ic++) if (f.chunks[ic]->is_mine()) {
      for (src_vol *sv = f.chunks[ic]->sources[D_stuff]; sv; sv = sv->next) {
	component c = direction_component(Ex, component_direction(sv->c));
	realnum *fr = f.chunks[ic]->f[c][0];
	realnum *fi = f.chunks[ic]->f[c][1];
	if (fr && fi) // complex E
	  for (size_t j=0; j<sv->npts; j++) {
	    const ptrdiff_t idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    EJ += complex<double>(fr[idx],fi[idx]) * conj(A);
	    Jsum += abs(A);
	  }
	else if (fr) { // E is purely real
	  for (size_t j=0; j<sv->npts; j++) {
	    const ptrdiff_t idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    EJ += double(fr[idx]) * conj(A);
	    Jsum += abs(A);
	  }
	}
      }
      for (src_vol *sv = f.chunks[ic]->sources[B_stuff]; sv; sv = sv->next) {
	component c = direction_component(Hx, component_direction(sv->c));
	realnum *fr = f.chunks[ic]->f[c][0];
	realnum *fi = f.chunks[ic]->f[c][1];
	if (fr && fi) // complex H
	  for (size_t j=0; j<sv->npts; j++) {
	    const ptrdiff_t idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    HJ += complex<double>(fr[idx],fi[idx]) * conj(A);
	    Jsum += abs(A);
	  }
	else if (fr) { // H is purely real
	  for (size_t j=0; j<sv->npts; j++) {
	    const ptrdiff_t idx = sv->index[j];
	    const complex<double> A = sv->A[j];
	    HJ += double(fr[idx]) * conj(A);
	    Jsum += abs(A);
	  }
	}
      }
    }
  for (int i = 0; i < Nomega; ++i) {
    complex<double> Ephase = polar(1.0, (omega_min+i*domega)*f.time())*scale;
    complex<double> Hphase = polar(1.0, (omega_min+i*domega)*(f.time()-f.dt/2))*scale;
    Fdft[i] += Ephase * EJ + Hphase * HJ;

    // NOTE: take only 1st time dependence: assumes all sources have same J(t)
    if (f.sources) {
      if (f.is_real) // todo: not quite right if A is complex
	Jdft[i] += Ephase * real(f.sources->current());
      else
	Jdft[i] += Ephase * f.sources->current();
    }
  }

  // correct for dV factors
  Jsum *= sqrt(f.gv.dV(f.gv.icenter(),1).computational_volume());

}

}

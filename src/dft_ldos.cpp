/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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

dft_ldos::dft_ldos(double freq_min, double freq_max, int Nfreq) {
  freq = meep::linspace(freq_min, freq_max, Nfreq);
  Fdft = new complex<double>[Nfreq];
  Jdft = new complex<double>[Nfreq];
  for (int i = 0; i < Nfreq; ++i)
    Fdft[i] = Jdft[i] = 0.0;
  Jsum = 1.0;
  saved_overall_scale = 1.0;
}

dft_ldos::dft_ldos(const std::vector<double> freq_) {
  const size_t Nfreq = freq_.size();
  freq = freq_;
  Fdft = new complex<double>[Nfreq];
  Jdft = new complex<double>[Nfreq];
  for (size_t i = 0; i < Nfreq; ++i)
    Fdft[i] = Jdft[i] = 0.0;
  Jsum = 1.0;
  saved_overall_scale = 1.0;
}

dft_ldos::dft_ldos(const double *freq_, size_t Nfreq) : freq(Nfreq) {
  for (size_t i = 0; i < Nfreq; ++i)
    freq[i] = freq_[i];
  Fdft = new complex<double>[Nfreq];
  Jdft = new complex<double>[Nfreq];
  for (size_t i = 0; i < Nfreq; ++i)
    Fdft[i] = Jdft[i] = 0.0;
  Jsum = 1.0;
  saved_overall_scale = 1.0;
}

// |c|^2
static double abs2(complex<double> c) { return real(c) * real(c) + imag(c) * imag(c); }

double *dft_ldos::ldos() {
  // we try to get the overall scale factor right (at least for a point source)
  // so that we can compare against the analytical formula for testing
  // ... in most practical cases, the scale factor won't matter because
  //     the user will compute the relative LDOS of 2 cases (e.g. LDOS/vacuum)

  // overall scale factor
  double Jsum_all = sum_to_all(Jsum);
  saved_overall_scale = 4.0 / pi                 // from definition of LDOS comparison to power
                        * -0.5                   // power = -1/2 Re[E* J]
                        / (Jsum_all * Jsum_all); // normalize to unit-integral current

  const size_t Nfreq = freq.size();
  double *sum = new double[Nfreq];
  for (size_t i = 0; i < Nfreq; ++i) /* 4/pi * work done by unit dipole */
    sum[i] = saved_overall_scale * real(Fdft[i] * conj(Jdft[i])) / abs2(Jdft[i]);
  double *out = new double[Nfreq];
  sum_to_all(sum, out, Nfreq);
  delete[] sum;
  return out;
}

complex<double> *dft_ldos::F() const {
  const size_t Nfreq = freq.size();
  complex<double> *out = new complex<double>[Nfreq];
  sum_to_all(Fdft, out, Nfreq);
  return out;
}

complex<double> *dft_ldos::J() const {
  const size_t Nfreq = freq.size();
  complex<double> *out = new complex<double>[Nfreq];
  // note: Jdft is the same on all processes, so no sum_to_all
  memcpy(out, Jdft, Nfreq * sizeof(complex<double>));
  return out;
}

void dft_ldos::update(fields &f) {
  complex<double> EJ = 0.0; // integral E * J*
  complex<double> HJ = 0.0; // integral H * J* for magnetic currents

  double scale = (f.dt / sqrt(2 * pi));

  // compute Jsum for LDOS normalization purposes
  // ...don't worry about the tiny inefficiency of recomputing this repeatedly
  Jsum = 0.0;

  for (int ic = 0; ic < f.num_chunks; ic++)
    if (f.chunks[ic]->is_mine()) {
      for (const src_vol &sv : f.chunks[ic]->get_sources(D_stuff)) {
        component c = direction_component(Ex, component_direction(sv.c));
        realnum *fr = f.chunks[ic]->f[c][0];
        realnum *fi = f.chunks[ic]->f[c][1];
        if (fr && fi) // complex E
          for (size_t j = 0; j < sv.num_points(); j++) {
            const ptrdiff_t idx = sv.index_at(j);
            const complex<double> &A = sv.amplitude_at(j);
            EJ += complex<double>(fr[idx], fi[idx]) * conj(A);
            Jsum += abs(A);
          }
        else if (fr) { // E is purely real
          for (size_t j = 0; j < sv.num_points(); j++) {
            const ptrdiff_t idx = sv.index_at(j);
            const complex<double> &A = sv.amplitude_at(j);
            EJ += double(fr[idx]) * conj(A);
            Jsum += abs(A);
          }
        }
      }
      for (const src_vol &sv : f.chunks[ic]->get_sources(B_stuff)) {
        component c = direction_component(Hx, component_direction(sv.c));
        realnum *fr = f.chunks[ic]->f[c][0];
        realnum *fi = f.chunks[ic]->f[c][1];
        if (fr && fi) // complex H
          for (size_t j = 0; j < sv.num_points(); j++) {
            const ptrdiff_t idx = sv.index_at(j);
            const complex<double> &A = sv.amplitude_at(j);
            HJ += complex<double>(fr[idx], fi[idx]) * conj(A);
            Jsum += abs(A);
          }
        else if (fr) { // H is purely real
          for (size_t j = 0; j < sv.num_points(); j++) {
            const ptrdiff_t idx = sv.index_at(j);
            const complex<double> &A = sv.amplitude_at(j);
            HJ += double(fr[idx]) * conj(A);
            Jsum += abs(A);
          }
        }
      }
    }
  for (size_t i = 0; i < freq.size(); ++i) {
    complex<double> Ephase = polar(1.0, 2 * pi * freq[i] * f.time()) * scale;
    complex<double> Hphase = polar(1.0, 2 * pi * freq[i] * (f.time() - f.dt / 2)) * scale;
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
  Jsum *= sqrt(f.gv.dV(f.gv.icenter(), 1).computational_volume());
}

} // namespace meep

/* Copyright (C) 2005-2019 Massachusetts Institute of Technology
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
#if defined(HAVE_LIBFFTW3)
#include <fftw3.h>
#elif defined(HAVE_LIBDFFTW)
#include <dfftw.h>
#elif defined(HAVE_LIBFFTW)
#include <fftw.h>
#endif
#if defined(HAVE_LIBFFTW3) || defined(HAVE_LIBFFTW) || defined(HAVE_LIBDFFTW)
#define HAVE_SOME_FFTW 1
#else
#define HAVE_SOME_FFTW 0
#endif

/* Below are the monitor point routines. */

using namespace std;

namespace meep {

monitor_point::monitor_point() { next = NULL; }

monitor_point::~monitor_point() {
  if (next) delete next;
}

inline complex<double> getcm(const realnum *const f[2], size_t i) {
  return complex<double>(f[0][i], f[1][i]);
}

void fields::get_point(monitor_point *pt, const vec &loc) const {
  if (pt == NULL) abort("Error:  get_point passed a null pointer!\n");
  for (int i = 0; i < 10; i++)
    pt->f[i] = 0.0;
  pt->loc = loc;
  pt->t = time();
  FOR_COMPONENTS(c) {
    if (gv.has_field(c)) pt->f[c] = get_field(c, loc);
  }
}

complex<double> fields::get_field(int c, const vec &loc) const {
  return (is_derived(c) ? get_field(derived_component(c), loc) : get_field(component(c), loc));
}

double fields::get_field(derived_component c, const vec &loc) const {
  component c1 = Ex, c2 = Ex;
  double sum = 0;
  switch (c) {
    case Sx:
    case Sy:
    case Sz:
    case Sr:
    case Sp:
      switch (c) {
        case Sx:
          c1 = Ey;
          c2 = Hz;
          break;
        case Sy:
          c1 = Ez;
          c2 = Hx;
          break;
        case Sz:
          c1 = Ex;
          c2 = Hy;
          break;
        case Sr:
          c1 = Ep;
          c2 = Hz;
          break;
        case Sp:
          c1 = Ez;
          c2 = Hr;
          break;
        default: break; // never
      }
      sum += real(conj(get_field(c1, loc)) * get_field(c2, loc));
      sum -= real(conj(get_field(direction_component(Ex, component_direction(c2)), loc)) *
                  get_field(direction_component(Hx, component_direction(c1)), loc));
      return sum;
    case EnergyDensity:
    case D_EnergyDensity:
    case H_EnergyDensity:
      if (c != H_EnergyDensity) FOR_ELECTRIC_COMPONENTS(c1) {
          if (gv.has_field(c1)) {
            c2 = direction_component(Dx, component_direction(c1));
            sum += real(conj(get_field(c1, loc)) * get_field(c2, loc));
          }
        }
      if (c != D_EnergyDensity) FOR_MAGNETIC_COMPONENTS(c1) {
          if (gv.has_field(c1)) {
            complex<double> f = get_field(c1, loc);
            sum += real(conj(f) * f);
          }
        }
      return sum * 0.5;
    default: abort("unknown derived_component in get_field");
  }
}

complex<double> fields::get_field(component c, const vec &loc, bool parallel) const {
  switch (c) {
    case Dielectric: return get_eps(loc);
    case Permeability: return get_mu(loc);
    default:
      ivec ilocs[8];
      double w[8];
      complex<double> val[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      gv.interpolate(c, loc, ilocs, w);
      for (int argh = 0; argh < 8 && w[argh]; argh++)
        val[argh] = w[argh] * get_field(c, ilocs[argh], false);
      complex<double> res = 0.0;
      for (int i = 0; i < 8; i++)
        res += val[i];
      return parallel ? sum_to_all(res) : res;
  }
}

complex<double> fields::get_field(component c, const ivec &origloc, bool parallel) const {
  ivec iloc = origloc;
  complex<double> kphase = 1.0;
  locate_point_in_user_volume(&iloc, &kphase);
  for (int sn = 0; sn < S.multiplicity(); sn++)
    for (int i = 0; i < num_chunks; i++)
      if (chunks[i]->gv.owns(S.transform(iloc, sn))) {
        complex<double> val = S.phase_shift(c, sn) * kphase *
               chunks[i]->get_field(S.transform(c, sn), S.transform(iloc, sn));
        return parallel ? sum_to_all(val) : val;
      }
  return 0.0;
}

complex<double> fields_chunk::get_field(component c, const ivec &iloc) const {
  if (is_mine())
    return f[c][0] ? (f[c][1] ? getcm(f[c], gv.index(c, iloc)) : f[c][0][gv.index(c, iloc)]) : 0.0;
  else
    return 0.0;
}

double fields::get_chi1inv(component c, direction d, const ivec &origloc, bool parallel) const {
  ivec iloc = origloc;
  complex<double> aaack = 1.0;
  locate_point_in_user_volume(&iloc, &aaack);
  for (int sn = 0; sn < S.multiplicity(); sn++)
    for (int i = 0; i < num_chunks; i++)
      if (chunks[i]->gv.owns(S.transform(iloc, sn))) {
        signed_direction ds = S.transform(d, sn);
        double val = chunks[i]->get_chi1inv(S.transform(c, sn), ds.d, S.transform(iloc, sn)) *
                     (ds.flipped ^ S.transform(component_direction(c), sn).flipped ? -1 : 1);
        return parallel ? sum_to_all(val) : val;
      }
  return 0.0;
}

double fields_chunk::get_chi1inv(component c, direction d, const ivec &iloc) const {
  if (is_mine())
    return s->chi1inv[c][d] ? s->chi1inv[c][d][gv.index(c, iloc)]
                           : (d == component_direction(c) ? 1.0 : 0);
  return 0.0;
}

double fields::get_chi1inv(component c, direction d, const vec &loc, bool parallel) const {
  ivec ilocs[8];
  double w[8];
  double val[8] = {0,0,0,0,0,0,0,0};
  gv.interpolate(c, loc, ilocs, w);
  for (int argh = 0; argh < 8 && w[argh] != 0; argh++)
    val[argh] = w[argh] * get_chi1inv(c, d, ilocs[argh], false);
  double res = 0.0;
  for (int i = 0; i < 8; i++)
    res += val[i];
  return parallel ? sum_to_all(res) : res;
}

double fields::get_eps(const vec &loc) const {
  double tr = 0;
  int nc = 0;
  FOR_ELECTRIC_COMPONENTS(c) {
    if (gv.has_field(c)) {
      tr += get_chi1inv(c, component_direction(c), loc, false);
      ++nc;
    }
  }
  return nc / sum_to_all(tr);
}

double fields::get_mu(const vec &loc) const {
  double tr = 0;
  int nc = 0;
  FOR_MAGNETIC_COMPONENTS(c) {
    if (gv.has_field(c)) {
      tr += get_chi1inv(c, component_direction(c), loc, false);
      ++nc;
    }
  }
  return nc / sum_to_all(tr);
}

double structure::get_chi1inv(component c, direction d, const ivec &origloc, bool parallel) const {
  ivec iloc = origloc;
  for (int sn = 0; sn < S.multiplicity(); sn++)
    for (int i = 0; i < num_chunks; i++)
      if (chunks[i]->gv.owns(S.transform(iloc, sn))) {
        signed_direction ds = S.transform(d, sn);
        double val = chunks[i]->get_chi1inv(S.transform(c, sn), ds.d, S.transform(iloc, sn)) *
               (ds.flipped ^ S.transform(component_direction(c), sn).flipped ? -1 : 1);
        return parallel ? sum_to_all(val) : val;
      }
  return 0.0;
}

double structure_chunk::get_chi1inv(component c, direction d, const ivec &iloc) const {
  if (is_mine())
    return chi1inv[c][d] ? chi1inv[c][d][gv.index(c, iloc)] : (d == component_direction(c) ? 1.0 : 0);
  return 0.0;
}

double structure::get_chi1inv(component c, direction d, const vec &loc, bool parallel) const {
  ivec ilocs[8];
  double w[8];
  double val[8] = {0,0,0,0,0,0,0,0};
  gv.interpolate(c, loc, ilocs, w);
  for (int argh = 0; argh < 8 && w[argh]; argh++)
    val[argh] = w[argh] * get_chi1inv(c, d, ilocs[argh], false);
  double res = 0.0;
  for (int i = 0; i < 8; i++)
    res += val[i];
  return parallel ? sum_to_all(res) : res;
}

double structure::get_eps(const vec &loc) const {
  double tr = 0;
  int nc = 0;
  FOR_ELECTRIC_COMPONENTS(c) {
    if (gv.has_field(c)) {
      tr += get_chi1inv(c, component_direction(c), loc, false);
      ++nc;
    }
  }
  return nc / sum_to_all(tr);
}

double structure::get_mu(const vec &loc) const {
  double tr = 0;
  int nc = 0;
  FOR_MAGNETIC_COMPONENTS(c) {
    if (gv.has_field(c)) {
      tr += get_chi1inv(c, component_direction(c), loc, false);
      ++nc;
    }
  }
  return nc / sum_to_all(tr);
}

monitor_point *fields::get_new_point(const vec &loc, monitor_point *the_list) const {
  monitor_point *p = new monitor_point();
  get_point(p, loc);
  p->next = the_list;
  return p;
}

complex<double> monitor_point::get_component(component w) { return f[w]; }

double monitor_point::poynting_in_direction(direction d) {
  direction d1 = cycle_direction(loc.dim, d, 1);
  direction d2 = cycle_direction(loc.dim, d, 2);

  // below Ex and Hx are used just to say that we want electric or magnetic component
  complex<double> E1 = get_component(direction_component(Ex, d1));
  complex<double> E2 = get_component(direction_component(Ex, d2));
  complex<double> H1 = get_component(direction_component(Hx, d1));
  complex<double> H2 = get_component(direction_component(Hx, d2));

  return (real(E1) * real(H2) - real(E2) * real(H1)) + (imag(E1) * imag(H2) - imag(E2) * imag(H1));
}

double monitor_point::poynting_in_direction(vec dir) {
  if (dir.dim != loc.dim) abort("poynting_in_direction: dir.dim != loc.dim\n");
  dir = dir / abs(dir);
  double result = 0.0;
  LOOP_OVER_DIRECTIONS(dir.dim, d) { result += dir.in_direction(d) * poynting_in_direction(d); }
  return result;
}

void monitor_point::fourier_transform(component w, complex<double> **a, complex<double> **f,
                                      int *numout, double fmin, double fmax, int maxbands) {
  int n = 1;
  monitor_point *p = next;
  double tmax = t, tmin = t;
  while (p) {
    n++;
    if (p->t > tmax) tmax = p->t;
    if (p->t < tmin) tmin = p->t;
    p = p->next;
  }
  p = this;
  complex<double> *d = new complex<double>[n];
  for (int i = 0; i < n; i++, p = p->next) {
    d[i] = p->get_component(w);
  }
  if (fmin > 0.0) { // Get rid of any static fields_chunk!
    complex<double> mean = 0.0;
    for (int i = 0; i < n; i++)
      mean += d[i];
    mean /= n;
    for (int i = 0; i < n; i++)
      d[i] -= mean;
  }
#if HAVE_SOME_FFTW
  if ((fmin > 0.0 || fmax > 0.0) && maxbands > 0) {
#else
  if ((fmin <= 0.0 && fmax <= 0.0) || maxbands <= 0) {
    maxbands = n;
    fmin = 0;
    fmax = (n - 1) * (1.0 / (tmax - tmin));
  }
#endif
    *a = new complex<double>[maxbands];
    *f = new complex<double>[maxbands];
    *numout = maxbands;
    delete[] d;
    for (int i = 0; i < maxbands; i++) {
      double df = (maxbands == 1) ? 0.0 : (fmax - fmin) / (maxbands - 1);
      (*f)[i] = fmin + i * df;
      (*a)[i] = 0.0;
      p = this;
      while (p) {
        double inside = 2 * pi * real((*f)[i]) * p->t;
        (*a)[i] += p->get_component(w) * complex<double>(cos(inside), sin(inside));
        p = p->next;
      }
      (*a)[i] /= (tmax - tmin);
    }
#if HAVE_SOME_FFTW
  } else {
    *numout = n;
    *a = new complex<double>[n];
    *f = d;
    fftw_complex *in = (fftw_complex *)d, *out = (fftw_complex *)*a;
    fftw_plan p;
#ifdef HAVE_LIBFFTW3
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
#else
    p = fftw_create_plan(n, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_one(p, in, out);
    fftw_destroy_plan(p);
#endif
    for (int i = 0; i < n; i++) {
      (*f)[i] = i * (1.0 / (tmax - tmin));
      if (real((*f)[i]) > 0.5 * n / (tmax - tmin)) (*f)[i] -= n / (tmax - tmin);
      (*a)[i] *= (tmax - tmin) / n;
    }
  }
#endif
}

void monitor_point::harminv(component w, complex<double> **a, complex<double> **f, int *numout,
                            double fmin, double fmax, int maxbands) {
  int n = 1;
  monitor_point *p = next;
  double tmax = t, tmin = t;
  while (p) {
    n++;
    if (p->t > tmax) tmax = p->t;
    if (p->t < tmin) tmin = p->t;
    p = p->next;
  }
  p = this;
  complex<double> *d = new complex<double>[n];
  for (int i = 0; i < n; i++, p = p->next) {
    d[i] = p->get_component(w);
  }
  *a = new complex<double>[n];
  double *f_re = new double[n];
  double *f_im = new double[n];
  *numout = do_harminv(d, n, (tmax - tmin) / (n - 1), fmin, fmax, maxbands, *a, f_re, f_im, NULL);
  *f = new complex<double>[*numout];
  for (int i = 0; i < *numout; i++)
    (*f)[i] = complex<double>(f_re[i], f_im[i]);
  delete[] f_re;
  delete[] f_im;
  delete[] d;
}

} // namespace meep

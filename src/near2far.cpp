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

/* Near-to-far field transformation: compute DFT of tangential fields on
   a "near" surface, and use these (via the equivalence principle) to
   compute the fields on a "far" surface via the homogeneous-medium Green's
   function in 2d or 3d. */

#include "meep_internals.hpp"
#include <assert.h>
#include "config.h"
#include <math.h>

using namespace std;

namespace meep {

dft_near2far::dft_near2far(dft_chunk *F_, double fmin, double fmax, int Nf, double eps_, double mu_,
                           const volume &where_, const direction periodic_d_[2],
                           const int periodic_n_[2], const double periodic_k_[2],
                           const double period_[2])
    : F(F_), eps(eps_), mu(mu_), where(where_) {
  freq = meep::linspace(fmin, fmax, Nf);
  for (int i = 0; i < 2; ++i) {
    periodic_d[i] = periodic_d_[i];
    periodic_n[i] = periodic_n_[i];
    periodic_k[i] = periodic_k_[i];
    period[i] = period_[i];
  }
}

dft_near2far::dft_near2far(dft_chunk *F_, const std::vector<double> &freq_, double eps_, double mu_,
                           const volume &where_, const direction periodic_d_[2],
                           const int periodic_n_[2], const double periodic_k_[2],
                           const double period_[2])
    : F(F_), eps(eps_), mu(mu_), where(where_) {
  freq = freq_;
  for (int i = 0; i < 2; ++i) {
    periodic_d[i] = periodic_d_[i];
    periodic_n[i] = periodic_n_[i];
    periodic_k[i] = periodic_k_[i];
    period[i] = period_[i];
  }
}

dft_near2far::dft_near2far(dft_chunk *F_, const double *freq_, size_t Nfreq, double eps_,
                           double mu_, const volume &where_, const direction periodic_d_[2],
                           const int periodic_n_[2], const double periodic_k_[2],
                           const double period_[2])
    : F(F_), eps(eps_), mu(mu_), where(where_) {
  freq.resize(Nfreq);
  for (size_t i = 0; i < Nfreq; ++i)
    freq[i] = freq_[i];
  for (int i = 0; i < 2; ++i) {
    periodic_d[i] = periodic_d_[i];
    periodic_n[i] = periodic_n_[i];
    periodic_k[i] = periodic_k_[i];
    period[i] = period_[i];
  }
}

dft_near2far::dft_near2far(const dft_near2far &f) : F(f.F), eps(f.eps), mu(f.mu), where(f.where) {
  freq = f.freq;
  for (int i = 0; i < 2; ++i) {
    periodic_d[i] = f.periodic_d[i];
    periodic_n[i] = f.periodic_n[i];
    periodic_k[i] = f.periodic_k[i];
    period[i] = f.period[i];
  }
}

void dft_near2far::remove() {
  while (F) {
    dft_chunk *nxt = F->next_in_dft;
    delete F;
    F = nxt;
  }
}

void dft_near2far::operator-=(const dft_near2far &st) {
  if (F && st.F) *F -= *st.F;
}

void dft_near2far::save_hdf5(h5file *file, const char *dprefix) {
  save_dft_hdf5(F, "F", file, dprefix);
}

void dft_near2far::load_hdf5(h5file *file, const char *dprefix) {
  load_dft_hdf5(F, "F", file, dprefix);
}

void dft_near2far::save_hdf5(fields &f, const char *fname, const char *dprefix,
                             const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::WRITE, prefix);
  save_hdf5(ff, dprefix);
  delete ff;
}

void dft_near2far::load_hdf5(fields &f, const char *fname, const char *dprefix,
                             const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::READONLY, prefix);
  load_hdf5(ff, dprefix);
  delete ff;
}

void dft_near2far::scale_dfts(complex<double> scale) {
  if (F) F->scale_dft(scale);
}

typedef void (*greenfunc)(std::complex<double> *EH, const vec &x, double freq, double eps,
                          double mu, const vec &x0, component c0, std::complex<double> f0);

/* Given the field f0 correponding to current-source component c0 at
   x0, compute the E/H fields EH[6] (6 components) at x for a frequency
   freq in the homogeneous 3d medium eps and mu.

   Adapted from code by M. T. Homer Reid in his SCUFF-EM package
   (file scuff-em/src/libs/libIncField/PointSource.cc), which is GPL v2+. */
void green3d(std::complex<double> *EH, const vec &x, double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0) {
  vec rhat = x - x0;
  double r = abs(rhat);
  rhat = rhat / r;

  if (rhat.dim != D3) meep::abort("wrong dimensionality in green3d");

  double n = sqrt(eps * mu);
  double k = 2 * pi * freq * n;
  std::complex<double> ikr = std::complex<double>(0.0, k * r);
  double ikr2 = -(k * r) * (k * r);
  /* note that SCUFF-EM computes the fields from the dipole moment p,
     whereas we need it from the current J = -i*omega*p, so our result
     is divided by -i*omega compared to SCUFF */
  std::complex<double> expfac = f0 * polar(k * n / (4 * pi * r), k * r + pi * 0.5);
  double Z = sqrt(mu / eps);

  vec p = zero_vec(rhat.dim);
  p.set_direction(component_direction(c0), 1);
  double pdotrhat = p & rhat;
  vec rhatcrossp = vec(rhat.y() * p.z() - rhat.z() * p.y(), rhat.z() * p.x() - rhat.x() * p.z(),
                       rhat.x() * p.y() - rhat.y() * p.x());

  /* compute the various scalar quantities in the point source formulae */
  std::complex<double> term1 = 1.0 - 1.0 / ikr + 1.0 / ikr2;
  std::complex<double> term2 = (-1.0 + 3.0 / ikr - 3.0 / ikr2) * pdotrhat;
  std::complex<double> term3 = (1.0 - 1.0 / ikr);

  /* now assemble everything based on source type */
  if (is_electric(c0)) {
    expfac /= eps;

    EH[0] = expfac * (term1 * p.x() + term2 * rhat.x());
    EH[1] = expfac * (term1 * p.y() + term2 * rhat.y());
    EH[2] = expfac * (term1 * p.z() + term2 * rhat.z());

    EH[3] = expfac * term3 * rhatcrossp.x() / Z;
    EH[4] = expfac * term3 * rhatcrossp.y() / Z;
    EH[5] = expfac * term3 * rhatcrossp.z() / Z;
  }
  else if (is_magnetic(c0)) {
    expfac /= mu;

    EH[0] = -expfac * term3 * rhatcrossp.x() * Z;
    EH[1] = -expfac * term3 * rhatcrossp.y() * Z;
    EH[2] = -expfac * term3 * rhatcrossp.z() * Z;

    EH[3] = expfac * (term1 * p.x() + term2 * rhat.x());
    EH[4] = expfac * (term1 * p.y() + term2 * rhat.y());
    EH[5] = expfac * (term1 * p.z() + term2 * rhat.z());
  }
  else
    meep::abort("unrecognized source type");
}

// hankel function J + iY
#if defined(HAVE_JN)
static std::complex<double> hankel(int n, double x) {
  return std::complex<double>(jn(n, x), yn(n, x));
}
#elif defined(HAVE_LIBGSL)
#include <gsl/gsl_sf_bessel.h>
static std::complex<double> hankel(int n, double x) {
  return std::complex<double>(gsl_sf_bessel_Jn(n, x), gsl_sf_bessel_Yn(n, x));
}
#else  /* !HAVE_LIBGSL */
static std::complex<double> hankel(int n, double x) {
  (void)n;
  (void)x; // unused
  meep::abort("GNU GSL library is required for Hankel functions");
}
#endif /* !HAVE_LIBGSL */

/* like green3d, but 2d Green's functions */
void green2d(std::complex<double> *EH, const vec &x, double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0) {
  vec rhat = x - x0;
  double r = abs(rhat);
  rhat = rhat / r;

  if (rhat.dim != D2) meep::abort("wrong dimensionality in green2d");

  double omega = 2 * pi * freq;
  double k = omega * sqrt(eps * mu);
  std::complex<double> ik = std::complex<double>(0.0, k);
  double kr = k * r;
  double Z = sqrt(mu / eps);
  std::complex<double> H0 = hankel(0, kr) * f0;
  std::complex<double> H1 = hankel(1, kr) * f0;
  std::complex<double> ikH1 = 0.25 * ik * H1;

  if (component_direction(c0) == meep::Z) {
    if (is_electric(c0)) { // Ez source
      EH[0] = EH[1] = 0.0;
      EH[2] = (-0.25 * omega * mu) * H0;

      EH[3] = -rhat.y() * ikH1;
      EH[4] = rhat.x() * ikH1;
      EH[5] = 0.0;
    }
    else /* (is_magnetic(c0)) */ { // Hz source
      EH[0] = rhat.y() * ikH1;
      EH[1] = -rhat.x() * ikH1;
      EH[2] = 0.0;

      EH[3] = EH[4] = 0.0;
      EH[5] = (-0.25 * omega * eps) * H0;
    }
  }
  else { /* in-plane source */
    std::complex<double> H2 = hankel(2, kr) * f0;

    vec p = zero_vec(rhat.dim);
    p.set_direction(component_direction(c0), 1);

    double pdotrhat = p & rhat;
    double rhatcrossp = rhat.x() * p.y() - rhat.y() * p.x();

    if (is_electric(c0)) { // Exy source
      EH[0] = -(rhat.x() * (pdotrhat / r * 0.25 * Z)) * H1 +
              (rhat.y() * (rhatcrossp * omega * mu * 0.125)) * (H0 - H2);
      EH[1] = -(rhat.y() * (pdotrhat / r * 0.25 * Z)) * H1 -
              (rhat.x() * (rhatcrossp * omega * mu * 0.125)) * (H0 - H2);
      EH[2] = 0.0;

      EH[3] = EH[4] = 0.0;
      EH[5] = -rhatcrossp * ikH1;
    }
    else /* (is_magnetic(c0)) */ { // Hxy source
      EH[0] = EH[1] = 0.0;
      EH[2] = rhatcrossp * ikH1;

      EH[3] = -(rhat.x() * (pdotrhat / r * 0.25 / Z)) * H1 +
              (rhat.y() * (rhatcrossp * omega * eps * 0.125)) * (H0 - H2);
      EH[4] = -(rhat.y() * (pdotrhat / r * 0.25 / Z)) * H1 -
              (rhat.x() * (rhatcrossp * omega * eps * 0.125)) * (H0 - H2);
      EH[5] = 0.0;
    }
  }
}

// cylindrical Green's function constructed by integrating green3d as the source
// term rotates around the z axis with exp(im*phi) dependence, integrated to a tolerance tol.
// (note: this is the Green's function divided by 2pi*x0.r(), to compensate for a 2piR factor
//  in the near2far add_dft weight.)
void greencyl(std::complex<double> *EH, const vec &x, double freq, double eps, double mu,
              const vec &x0, component c0, std::complex<double> f0, double m, double tol) {
  if (x0.dim != Dcyl) meep::abort("wrong dimensionality in greencyl");
  vec x_3d(x.dim == Dcyl ? x.r() : x.x(), x.y(), x.z());
  direction d = component_direction(c0);
  component cx = direction_component(c0, X), cy = direction_component(c0, Y);
  for (int j = 0; j < 6; ++j)
    EH[j] = 0;

  /* Perform phi integral.  Since phi integrand is smooth, quadrature with equally spaced points
     should converge exponentially fast with the number N of quadrature points.  We
     repeatedly double N until convergence to tol is achieved, re-using previous points. */
  const int N0 = 4;
  double dphi = 2.0 / N0; // factor of 2*pi*r is already included in add_dft weight
  for (int N = N0; N <= 65536; N *= 2) {
    std::complex<double> EH_sum[6];
    dphi *= 0.5; // delta phi is halved because N doubles
    double dphi2pi = dphi * 2 * pi;
    for (int j = 0; j < 6; ++j)
      EH_sum[j] = EH[j] * 0.5; // re-use previous quadrature points (with halved dphi)
    /* N-point quadrature points i = 0..N-1.  After the first iteration (N==N0), we
       only need to sum over odd i, since the even i were summed for the previous N. */
    for (int i = (N > N0); i < N; i += 1 + (N > N0)) {
      double phi = i * dphi2pi, c = cos(phi), s = sin(phi);
      vec x0_phi(x0.r() * c, x0.r() * s, x0.z()); // source point rotated by phi
      std::complex<double> EH_phi[6], f0_exp_imphi = f0 * std::polar(1.0, m * phi) * dphi;
      /* if the source direction is in the r or phi directions, then we must rotate
        the direction of the source current in the xy plane as well */
      if (d == Z) { // source currents in z direction don't rotate
        green3d(EH_phi, x_3d, freq, eps, mu, x0_phi, c0, f0_exp_imphi);
        for (int j = 0; j < 6; ++j)
          EH_sum[j] += EH_phi[j];
      }
      else if (d == R) { // r_hat = c x_hat + s y_hat
        green3d(EH_phi, x_3d, freq, eps, mu, x0_phi, cx, f0_exp_imphi * c);
        for (int j = 0; j < 6; ++j)
          EH_sum[j] += EH_phi[j];
        green3d(EH_phi, x_3d, freq, eps, mu, x0_phi, cy, f0_exp_imphi * s);
        for (int j = 0; j < 6; ++j)
          EH_sum[j] += EH_phi[j];
      }
      else { // (d == P):  phi_hat = c y_hat - s x_hat
        green3d(EH_phi, x_3d, freq, eps, mu, x0_phi, cx, f0_exp_imphi * (-s));
        for (int j = 0; j < 6; ++j)
          EH_sum[j] += EH_phi[j];
        green3d(EH_phi, x_3d, freq, eps, mu, x0_phi, cy, f0_exp_imphi * c);
        for (int j = 0; j < 6; ++j)
          EH_sum[j] += EH_phi[j];
      }
    }
    // accumulate the new and old sums and check how much the integral has changed in L1 norm
    double sumdiff = 0, sumabs = 0;
    for (int j = 0; j < 6; ++j) {
      sumdiff += abs(EH[j] - EH_sum[j]);
      sumabs += abs(EH_sum[j]);
      EH[j] = EH_sum[j];
    }
    if (sumdiff <= sumabs * tol) break; // doubling N changed sum by less than tol
  }
}

void dft_near2far::farfield_lowlevel(std::complex<double> *EH, const vec &x) {
  if (x.dim != D3 && x.dim != D2 && x.dim != Dcyl)
    meep::abort("only 2d or 3d or cylindrical far-field computation is supported");
  greenfunc green = x.dim == D2 ? green2d : green3d;

  const size_t Nfreq = freq.size();
  for (size_t i = 0; i < 6 * Nfreq; ++i)
    EH[i] = 0.0;

  for (dft_chunk *f = F; f; f = f->next_in_dft) {
    assert(Nfreq == f->omega.size());

    component c0 = component(f->vc); /* equivalent source component */

    vec rshift(f->shift * (0.5 * f->fc->gv.inva));
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0; i < Nfreq; ++i) {
      std::complex<double> EH6[6];
      size_t idx_dft = 0;
      LOOP_OVER_IVECS(f->fc->gv, f->is, f->ie, idx) {
        IVEC_LOOP_LOC(f->fc->gv, x0);
        x0 = f->S.transform(x0, f->sn) + rshift;
        vec xs(x0);
        for (int i0 = -periodic_n[0]; i0 <= periodic_n[0]; ++i0) {
          if (periodic_d[0] != NO_DIRECTION)
            xs.set_direction(periodic_d[0], x0.in_direction(periodic_d[0]) + i0 * period[0]);
          double phase0 = i0 * periodic_k[0];
          for (int i1 = -periodic_n[1]; i1 <= periodic_n[1]; ++i1) {
            if (periodic_d[1] != NO_DIRECTION)
              xs.set_direction(periodic_d[1], x0.in_direction(periodic_d[1]) + i1 * period[1]);
            double phase = phase0 + i1 * periodic_k[1];
            std::complex<double> cphase = std::polar(1.0, phase);
            if (x.dim == Dcyl)
              greencyl(EH6, x, freq[i], eps, mu, xs, c0, f->dft[Nfreq * idx_dft + i], f->fc->m,
                       1e-3);
            else
              green(EH6, x, freq[i], eps, mu, xs, c0, f->dft[Nfreq * idx_dft + i]);
            for (int j = 0; j < 6; ++j)
              EH[i * 6 + j] += EH6[j] * cphase;
          }
        }
        idx_dft++;
      }
    }
  }
}

std::complex<double> *dft_near2far::farfield(const vec &x) {
  std::complex<double> *EH, *EH_local;
  const size_t Nfreq = freq.size();
  EH_local = new std::complex<double>[6 * Nfreq];
  farfield_lowlevel(EH_local, x);
  EH = new std::complex<double>[6 * Nfreq];
  sum_to_all(EH_local, EH, 6 * Nfreq);
  delete[] EH_local;
  return EH;
}

double *dft_near2far::get_farfields_array(const volume &where, int &rank, size_t *dims, size_t &N,
                                          double resolution) {
  /* compute output grid size etc. */
  double dx[3] = {0, 0, 0};
  direction dirs[3] = {X, Y, Z};

  rank = 0;
  N = dims[0] = dims[1] = dims[2] = 1;

  LOOP_OVER_DIRECTIONS(where.dim, d) {
    dims[rank] = int(floor(where.in_direction(d) * resolution));
    if (dims[rank] <= 1)
      dims[rank] = 1;
    else
      dx[rank] = where.in_direction(d) / (dims[rank] - 1);
    N *= dims[rank];
    dirs[rank++] = d;
  }
  if (where.dim == Dcyl) dirs[2] = P; // otherwise Z is listed twice

  const size_t Nfreq = freq.size();
  if (N * Nfreq < 1) return NULL; /* nothing to output */

  /* 6 x 2 x N x Nfreq array of fields in row-major order */
  double *EH = new double[6 * 2 * N * Nfreq];
  double *EH_ = new double[6 * 2 * N * Nfreq]; // temp array for sum_to_all

  /* fields for farfield_lowlevel for a single output point x */
  std::complex<double> *EH1 = new std::complex<double>[6 * Nfreq];

  double start = wall_time();
  size_t last_point = 0;

  vec x(where.dim);
  for (size_t i0 = 0; i0 < dims[0]; ++i0) {
    x.set_direction(dirs[0], where.in_direction_min(dirs[0]) + i0 * dx[0]);
    for (size_t i1 = 0; i1 < dims[1]; ++i1) {
      x.set_direction(dirs[1], where.in_direction_min(dirs[1]) + i1 * dx[1]);
      for (size_t i2 = 0; i2 < dims[2]; ++i2) {
        x.set_direction(dirs[2], where.in_direction_min(dirs[2]) + i2 * dx[2]);
        double t;
        if (verbosity > 0 && (t = wall_time()) > start + MEEP_MIN_OUTPUT_TIME) {
          size_t this_point = (dims[1] * i0 + i1) * dims[2] + i2 + 1;
          master_printf("get_farfields_array working on point %zu of %zu (%d%% done), %g s/point\n",
                        this_point, N, (int)((double)this_point / N * 100),
                        (t - start) / (std::max(1, (int)(this_point - last_point))));
          start = t;
          last_point = this_point;
        }
        farfield_lowlevel(EH1, x);
        if (verbosity > 1) all_wait(); // Allow consistent progress updates from master
        ptrdiff_t idx = (i0 * dims[1] + i1) * dims[2] + i2;
        for (size_t i = 0; i < Nfreq; ++i)
          for (int k = 0; k < 6; ++k) {
            EH_[((k * 2 + 0) * N + idx) * Nfreq + i] = real(EH1[i * 6 + k]);
            EH_[((k * 2 + 1) * N + idx) * Nfreq + i] = imag(EH1[i * 6 + k]);
          }
      }
    }
  }
  sum_to_all(EH_, EH, 6 * 2 * N * Nfreq);

  /* collapse singleton dimensions */
  int ireduced = 0;
  for (int i = 0; i < rank; ++i) {
    if (dims[i] > 1) dims[ireduced++] = dims[i];
  }
  rank = ireduced;

  delete[] EH_;
  delete[] EH1;
  return EH;
}

void dft_near2far::save_farfields(const char *fname, const char *prefix, const volume &where,
                                  double resolution) {
  size_t dims[4] = {1, 1, 1, 1};
  int rank = 0;
  size_t N = 1;

  double *EH = get_farfields_array(where, rank, dims, N, resolution);
  if (!EH) return; /* nothing to output */

  const size_t Nfreq = freq.size();
  /* frequencies are the last dimension */
  if (Nfreq > 1) dims[rank++] = Nfreq;

  /* output to a file with one dataset per component & real/imag part */
  if (am_master()) {
    const int buflen = 1024;
    static char filename[buflen];
    snprintf(filename, buflen, "%s%s%s.h5", prefix ? prefix : "", prefix && prefix[0] ? "-" : "",
             fname);
    h5file ff(filename, h5file::WRITE, false);
    component c[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
    char dataname[128];
    for (int k = 0; k < 6; ++k)
      for (int reim = 0; reim < 2; ++reim) {
        snprintf(dataname, 128, "%s.%c", component_name(c[k]), "ri"[reim]);
        ff.write(dataname, rank, dims, EH + (k * 2 + reim) * N * Nfreq);
      }
  }

  delete[] EH;
}

double *dft_near2far::flux(direction df, const volume &where, double resolution) {
  if (coordinate_mismatch(where.dim, df) || where.dim == Dcyl)
    meep::abort("cannot get flux for near2far: co-ordinate mismatch");

  size_t dims[4] = {1, 1, 1, 1};
  int rank = 0;
  size_t N = 1;

  double *EH = get_farfields_array(where, rank, dims, N, resolution);

  const size_t Nfreq = freq.size();
  double *F = new double[Nfreq];
  std::complex<double> ff_EH[6];
  std::complex<double> cE[2], cH[2];

  for (size_t i = 0; i < Nfreq; ++i)
    F[i] = 0;

  for (size_t idx = 0; idx < N; ++idx) {
    for (size_t i = 0; i < Nfreq; ++i) {
      for (int k = 0; k < 6; ++k)
        ff_EH[k] = std::complex<double>(*(EH + ((k * 2 + 0) * N + idx) * Nfreq + i),
                                        *(EH + ((k * 2 + 1) * N + idx) * Nfreq + i));
      switch (df) {
        case X: cE[0] = ff_EH[1], cE[1] = ff_EH[2], cH[0] = ff_EH[5], cH[1] = ff_EH[4]; break;
        case Y: cE[0] = ff_EH[2], cE[1] = ff_EH[0], cH[0] = ff_EH[3], cH[1] = ff_EH[5]; break;
        case Z: cE[0] = ff_EH[0], cE[1] = ff_EH[1], cH[0] = ff_EH[4], cH[1] = ff_EH[3]; break;
        case R:
        case P:
        case NO_DIRECTION: meep::abort("invalid flux direction");
      }
      for (int j = 0; j < 2; ++j)
        F[i] += real(cE[j] * conj(cH[j])) * (1 - 2 * j);
    }
  }

  double dV = 1;
  LOOP_OVER_DIRECTIONS(where.dim, d) {
    int dim = int(floor(where.in_direction(d) * resolution));
    if (dim > 1) dV *= where.in_direction(d) / (dim - 1);
  }

  for (size_t i = 0; i < Nfreq; ++i)
    F[i] *= dV;

  delete[] EH;

  return F;
}

static double approxeq(double a, double b) { return fabs(a - b) < 0.5e-11 * (fabs(a) + fabs(b)); }

dft_near2far fields::add_dft_near2far(const volume_list *where, const double *freq, size_t Nfreq,
                                      int decimation_factor, int Nperiods) {

  dft_chunk *F = 0; /* E and H chunks*/
  double eps = 0, mu = 0;
  volume everywhere = where->v;

  direction periodic_d[2] = {NO_DIRECTION, NO_DIRECTION};
  int periodic_n[2] = {0, 0};
  double periodic_k[2] = {0, 0}, period[2] = {0, 0};

  for (const volume_list *w = where; w; w = w->next) {
    everywhere = everywhere | where->v;
    direction nd = component_direction(w->c);
    if (nd == NO_DIRECTION) nd = normal_direction(w->v);
    if (nd == NO_DIRECTION) meep::abort("unknown dft_near2far normal");
    direction fd[2];

    double weps = real(get_eps(w->v.center()));
    double wmu = real(get_mu(w->v.center()));
    if (w != where && !(approxeq(eps, weps) && approxeq(mu, wmu)))
      meep::abort("dft_near2far requires surfaces in a homogeneous medium");
    eps = weps;
    mu = wmu;

    /* two transverse directions to normal (in cyclic order to get
       correct sign s below) */
    switch (nd) {
      case X:
        fd[0] = Y;
        fd[1] = Z;
        break;
      case Y:
        fd[0] = Z;
        fd[1] = X;
        break;
      case R:
        fd[0] = P;
        fd[1] = Z;
        break;
      case P:
        fd[0] = Z;
        fd[1] = R;
        break;
      case Z:
        if (gv.dim == Dcyl)
          fd[0] = R, fd[1] = P;
        else
          fd[0] = X, fd[1] = Y;
        break;
      default: meep::abort("invalid normal direction in dft_near2far!");
    }

    if (Nperiods > 1) {
      for (int i = 0; i < 2; ++i) {
        double user_width = user_volume.num_direction(fd[i]) / a;
        if (has_direction(v.dim, fd[i]) && boundaries[High][fd[i]] == Periodic &&
            boundaries[Low][fd[i]] == Periodic &&
            float(w->v.in_direction(fd[i])) >= float(user_width)) {
          periodic_d[i] = fd[i];
          periodic_n[i] = Nperiods;
          period[i] = user_width;
          periodic_k[i] = 2 * pi * real(k[fd[i]]) * period[i];
        }
      }
    }

    for (int i = 0; i < 2; ++i) {   /* E or H */
      for (int j = 0; j < 2; ++j) { /* first or second component */
        component c = direction_component(i == 0 ? Ex : Hx, fd[j]);

        /* find equivalent source component c0 and sign s */
        component c0 = direction_component(i == 0 ? Hx : Ex, fd[1 - j]);
        double s = j == 0 ? 1 : -1; /* sign of n x c */
        if (is_electric(c)) s = -s;

        F = add_dft(c, w->v, freq, Nfreq, true, s * w->weight, F, false, 1.0, false, c0,
                    decimation_factor);
      }
    }
  }

  return dft_near2far(F, freq, Nfreq, eps, mu, everywhere, periodic_d, periodic_n, periodic_k,
                      period);
}

// Modified from farfield_lowlevel
std::vector<struct sourcedata> dft_near2far::near_sourcedata(const vec &x_0, double *farpt_list,
                                                             size_t nfar_pts,
                                                             const std::complex<double> *dJ) {
  if (x_0.dim != D3 && x_0.dim != D2 && x_0.dim != Dcyl)
    meep::abort("only 2d or 3d or cylindrical far-field computation is supported");
  greenfunc green = x_0.dim == D2 ? green2d : green3d;

  const size_t Nfreq = freq.size();
  std::vector<struct sourcedata> temp;

  for (dft_chunk *f = F; f; f = f->next_in_dft) {
    assert(Nfreq == f->omega.size());
    std::vector<ptrdiff_t> idx_arr;
    std::vector<std::complex<double> > amp_arr;
    component c0 = component(f->vc); /* equivalent source component */

    vec rshift(f->shift * (0.5 * f->fc->gv.inva));
    std::complex<double> EH6[6];
    size_t idx_dft = 0;
    sourcedata temp_struct = {component(f->c), idx_arr, f->fc->chunk_idx, amp_arr};

    LOOP_OVER_IVECS(f->fc->gv, f->is, f->ie, idx) {
      IVEC_LOOP_ILOC(f->fc->gv, ix0);
      IVEC_LOOP_LOC(f->fc->gv, x0);
      x0 = f->S.transform(x0, f->sn) + rshift;
      vec xs(x0);
      double w;
      w = IVEC_LOOP_WEIGHT(f->s0, f->s1, f->e0, f->e1, f->dV0 + f->dV1 * loop_i2);

      temp_struct.idx_arr.push_back(idx);
      for (size_t i = 0; i < Nfreq; ++i) {
        std::complex<double> EH0 = std::complex<double>(0, 0);
        for (int i0 = -periodic_n[0]; i0 <= periodic_n[0]; ++i0) {
          if (periodic_d[0] != NO_DIRECTION)
            xs.set_direction(periodic_d[0], x0.in_direction(periodic_d[0]) + i0 * period[0]);
          double phase0 = i0 * periodic_k[0];
          for (int i1 = -periodic_n[1]; i1 <= periodic_n[1]; ++i1) {
            if (periodic_d[1] != NO_DIRECTION)
              xs.set_direction(periodic_d[1], x0.in_direction(periodic_d[1]) + i1 * period[1]);
            double phase = phase0 + i1 * periodic_k[1];
            std::complex<double> cphase = std::polar(1.0, phase);
            for (size_t ipt = 0; ipt < nfar_pts; ++ipt) {
              vec x = vec(farpt_list[3 * ipt], farpt_list[3 * ipt + 1], farpt_list[3 * ipt + 2]);
              if (x_0.dim == Dcyl)
                greencyl(EH6, x, freq[i], eps, mu, xs, c0, w, f->fc->m, 1e-3);
              else
                green(EH6, x, freq[i], eps, mu, xs, c0, w);
              for (int j = 0; j < 6; ++j)
                EH0 += EH6[j] * cphase * (f->stored_weight) * dJ[6 * Nfreq * ipt + 6 * i + j];
            }
          }
        }
        idx_dft++;
        if (is_electric(temp_struct.near_fd_comp)) EH0 *= -1;
        EH0 /= f->S.multiplicity(ix0);
        if (x_0.dim == Dcyl) {
          if (xs.r() != 0) EH0 /= xs.r();
          // Somehow, a factor of (2pi)r for r of the near2far region was double counted.
          // 2pi is canceled out by the integral over design region, where an extra factor of 2pi*r'
          // (r' of the design region) is needed. See meepgeom.cpp
        }
        temp_struct.amp_arr.push_back(EH0);
      }
    }
    temp.push_back(temp_struct);
  }
  return temp;
}

} // namespace meep

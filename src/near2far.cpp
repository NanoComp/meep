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

/* Near-to-far field transformation: compute DFT of tangential fields on
   a "near" surface, and use these (via the equivalence principle) to
   compute the fields on a "far" surface via the homogeneous-medium Green's
   function in 2d or 3d. */

#include <meep.hpp>
#include <assert.h>
#include "config.h"
#include <math.h>

using namespace std;

namespace meep {

dft_near2far::dft_near2far(dft_chunk *F_,
                           double fmin, double fmax, int Nf,
                           double eps_, double mu_)
{
  if (Nf <= 1) fmin = fmax = (fmin + fmax) * 0.5;
  freq_min = fmin;
  Nfreq = Nf;
  dfreq = Nf <= 1 ? 0.0 : (fmax - fmin) / (Nf - 1);
  F = F_;
  eps = eps_; mu = mu_;
}

dft_near2far::dft_near2far(const dft_near2far &f) {
  freq_min = f.freq_min; Nfreq = f.Nfreq; dfreq = f.dfreq;
  F = f.F;
  eps = f.eps; mu = f.mu;
}

void dft_near2far::remove()
{
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

typedef void (*greenfunc)(std::complex<double> *EH, const vec &x,
                          double freq, double eps, double mu,
                          const vec &x0, component c0, std::complex<double>);

/* Given the field f0 correponding to current-source component c0 at
   x0, compute the E/H fields EH[6] (6 components) at x for a frequency
   freq in the homogeneous 3d medium eps and mu. 

   Adapted from code by M. T. Homer Reid in his SCUFF-EM package
   (file scuff-em/src/libs/libIncField/PointSource.cc), which is GPL v2+. */   
void green3d(std::complex<double> *EH, const vec &x,
             double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0)
{
    vec rhat = x - x0;
    double r = abs(rhat);
    rhat = rhat / r;

    if (rhat.dim != D3) abort("wrong dimensionality in green3d");

    double n = sqrt(eps*mu);
    double k = 2*pi*freq*n;
    std::complex<double> ikr = std::complex<double>(0.0, k*r);
    double ikr2   = -(k*r)*(k*r);
    /* note that SCUFF-EM computes the fields from the dipole moment p,
       whereas we need it from the current J = -i*omega*p, so our result
       is divided by -i*omega compared to SCUFF */
    std::complex<double> expfac = f0 * polar(k*n/(4*pi*r), k*r + pi*0.5);
    double Z = sqrt(mu/eps);

    vec p = zero_vec(rhat.dim);
    p.set_direction(component_direction(c0), 1);
    double pdotrhat = p & rhat;
    vec rhatcrossp = vec(rhat.y() * p.z() -
                         rhat.z() * p.y(),
                         rhat.z() * p.x() -
                         rhat.x() * p.z(),
                         rhat.x() * p.y() -
                         rhat.y() * p.x());
   
    /* compute the various scalar quantities in the point source formulae */
    std::complex<double> term1 =  1.0 - 1.0/ikr + 1.0/ikr2;
    std::complex<double> term2 = (-1.0 + 3.0/ikr - 3.0/ikr2) * pdotrhat;
    std::complex<double> term3 = (1.0 - 1.0/ikr);
 
    /* now assemble everything based on source type */
    if (is_electric(c0)) {
        expfac /= eps;

        EH[0] = expfac * (term1*p.x() + term2*rhat.x());
        EH[1] = expfac * (term1*p.y() + term2*rhat.y());
        EH[2] = expfac * (term1*p.z() + term2*rhat.z());
        
        EH[3] = expfac*term3*rhatcrossp.x() / Z;
        EH[4] = expfac*term3*rhatcrossp.y() / Z;
        EH[5] = expfac*term3*rhatcrossp.z() / Z;
    }
    else if (is_magnetic(c0)) {
        expfac /= mu;

        EH[0] = -expfac*term3*rhatcrossp.x() * Z;
        EH[1] = -expfac*term3*rhatcrossp.y() * Z;
        EH[2] = -expfac*term3*rhatcrossp.z() * Z;
        
        EH[3] = expfac * (term1*p.x() + term2*rhat.x());
        EH[4] = expfac * (term1*p.y() + term2*rhat.y());
        EH[5] = expfac * (term1*p.z() + term2*rhat.z());
    }
    else
        abort("unrecognized source type");
}

#ifdef HAVE_LIBGSL
#  include <gsl/gsl_sf_bessel.h>
// hankel function J + iY
static std::complex<double> hankel(int n, double x) {
    return std::complex<double>(gsl_sf_bessel_Jn(n, x),
                                gsl_sf_bessel_Yn(n, x));
}
#else /* !HAVE_LIBGSL */
static std::complex<double> hankel(int n, double x) {
    (void) n; (void) x; // unused
    abort("GNU GSL library is required for Hankel functions");
}
#endif /* !HAVE_LIBGSL */

/* like green3d, but 2d Green's functions */
void green2d(std::complex<double> *EH, const vec &x,
             double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0)
{
    vec rhat = x - x0;
    double r = abs(rhat);
    rhat = rhat / r;

    if (rhat.dim != D2) abort("wrong dimensionality in green2d");

    double omega = 2*pi*freq;
    double k = omega*sqrt(eps*mu);
    std::complex<double> ik = std::complex<double>(0.0, k);
    double kr = k*r;
    double Z = sqrt(mu/eps);
    std::complex<double> H0 = hankel(0, kr) * f0;
    std::complex<double> H1 = hankel(1, kr) * f0;
    std::complex<double> ikH1 = 0.25 * ik * H1;

    if (component_direction(c0) == meep::Z) {
        if (is_electric(c0)) { // Ez source
            EH[0] = EH[1] = 0.0;
            EH[2] = (-0.25*omega*mu) * H0;

            EH[3] = -rhat.y() * ikH1;
            EH[4] =  rhat.x() * ikH1;
            EH[5] = 0.0;
        }
        else /* (is_magnetic(c0)) */ { // Hz source
            EH[0] =  rhat.y() * ikH1;
            EH[1] = -rhat.x() * ikH1;
            EH[2] = 0.0;

            EH[3] = EH[4] = 0.0;
            EH[5] = (-0.25*omega*eps) * H0;
        }
    }
    else { /* in-plane source */
        std::complex<double> H2 = hankel(2, kr) * f0;

        vec p = zero_vec(rhat.dim);
        p.set_direction(component_direction(c0), 1);

        double pdotrhat = p & rhat;
        double rhatcrossp = rhat.x() * p.y() - rhat.y() * p.x();

        if (is_electric(c0)) { // Exy source
            EH[0] = -(rhat.x() * (pdotrhat/r * 0.25*Z)) * H1 +
                (rhat.y() * (rhatcrossp * omega*mu * 0.125)) * (H0 - H2);
            EH[1] = -(rhat.y() * (pdotrhat/r * 0.25*Z)) * H1 -
                (rhat.x() * (rhatcrossp * omega*mu * 0.125)) * (H0 - H2);
            EH[2] = 0.0;

            EH[3] = EH[4] = 0.0;
            EH[5] = -rhatcrossp * ikH1;
        }
        else /* (is_magnetic(c0)) */ { // Hxy source
            EH[0] = EH[1] = 0.0;
            EH[2] = rhatcrossp * ikH1;
            
            EH[3] = -(rhat.x() * (pdotrhat/r * 0.25/Z)) * H1 +
                (rhat.y() * (rhatcrossp * omega*eps * 0.125)) * (H0 - H2);
            EH[4] = -(rhat.y() * (pdotrhat/r * 0.25/Z)) * H1 -
                (rhat.x() * (rhatcrossp * omega*eps * 0.125)) * (H0 - H2);
            EH[5] = 0.0;
        }
    }
}

void dft_near2far::farfield_lowlevel(std::complex<double> *EH, const vec &x)
{
    if (x.dim != D3 && x.dim != D2)
        abort("only 2d or 3d far-field computation is supported");
    greenfunc green = x.dim == D2 ? green2d : green3d;

    std::complex<double> EH6[6];
    for (int i = 0; i < 6 * Nfreq; ++i)
        EH[i] = 0.0;
    
    for (dft_chunk *f = F; f; f = f->next_in_dft) {
        assert(Nfreq == f->Nomega);

        component c0 = component(f->vc); /* equivalent source component */

        vec rshift(f->shift * (0.5*f->fc->gv.inva));
        int idx_dft = 0;
        LOOP_OVER_IVECS(f->fc->gv, f->is, f->ie, idx) {
            IVEC_LOOP_LOC(f->fc->gv, x0);
            x0 = f->S.transform(x0, f->sn) + rshift;
            for (int i = 0; i < Nfreq; ++i) {
                double freq = freq_min + i*dfreq;
                green(EH6, x, freq, eps, mu, x0, c0, f->dft[Nfreq*idx_dft+i]);
                for (int j = 0; j < 6; ++j) EH[i*6 + j] += EH6[j];
            }
            idx_dft++;
        }
    }
}

std::complex<double> *dft_near2far::farfield(const vec &x) {
    std::complex<double> *EH, *EH_local;
    EH_local = new std::complex<double>[6 * Nfreq];
    farfield_lowlevel(EH_local, x);
    EH = new std::complex<double>[6 * Nfreq];
    sum_to_all(EH_local, EH, 6 * Nfreq);
    delete[] EH_local;
    return EH;
}

void dft_near2far::save_farfields(const char *fname, const char *prefix,
                                  const volume &where, double resolution) {
    /* compute output grid size etc. */
    int dims[4] = {1,1,1,1};
    double dx[3] = {0,0,0};
    direction dirs[3] = {X,Y,Z};
    int rank = 0, N = 1;
    LOOP_OVER_DIRECTIONS(where.dim, d) {
        dims[rank] = int(floor(where.in_direction(d) * resolution));
        if (dims[rank] <= 1) {
            dims[rank] = 1;
            dx[rank] = 0;
        }
        else
            dx[rank] = where.in_direction(d) / (dims[rank] - 1);
        N *= dims[rank];
        dirs[rank++] = d;
    }

    if (N * Nfreq < 1) return; /* nothing to output */

    /* 6 x 2 x N x Nfreq array of fields in row-major order */
    realnum *EH = new realnum[6*2*N*Nfreq];
    realnum *EH_ = new realnum[6*2*N*Nfreq]; // temp array for sum_to_master

    /* fields for farfield_lowlevel for a single output point x */
    std::complex<double> *EH1 = new std::complex<double>[6*Nfreq];

    vec x(where.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
        x.set_direction(dirs[0], where.in_direction_min(dirs[0]) + i0*dx[0]);
        for (int i1 = 0; i1 < dims[1]; ++i1) {
            x.set_direction(dirs[1], 
                            where.in_direction_min(dirs[1]) + i1*dx[1]);
            for (int i2 = 0; i2 < dims[2]; ++i2) {
                x.set_direction(dirs[2], 
                                where.in_direction_min(dirs[2]) + i2*dx[2]);
                farfield_lowlevel(EH1, x);
                int idx = (i0 * dims[1] + i1) * dims[2] + i2;
                for (int i = 0; i < Nfreq; ++i)
                    for (int k = 0; k < 6; ++k) {
                        EH_[((k * 2 + 0) * N + idx) * Nfreq + i] =
                            real(EH1[i * 6 + k]);
                        EH_[((k * 2 + 1) * N + idx) * Nfreq + i] =
                            imag(EH1[i * 6 + k]);
                    }
            }
        }
    }

    delete[] EH1;
    sum_to_master(EH_, EH, 6*2*N*Nfreq);
    delete[] EH_;

    /* collapse trailing singleton dimensions */
    while (rank > 0 && dims[rank-1] == 1)
        --rank;
    /* frequencies are the last dimension */
    if (Nfreq > 1)
        dims[rank++] = Nfreq;

    /* output to a file with one dataset per component & real/imag part */
    if (am_master()) {
        const int buflen = 1024;
        static char filename[buflen];
        snprintf(filename, buflen, "%s%s%s.h5",
                 prefix ? prefix : "", prefix && prefix[0] ? "-" : "",
                 fname);
        h5file ff(filename, h5file::WRITE, false);
        component c[6] = {Ex,Ey,Ez,Hx,Hy,Hz};
        char dataname[128];
        for (int k = 0; k < 6; ++k)
            for (int reim = 0; reim < 2; ++reim) {
                snprintf(dataname, 128, "%s.%c", 
                         component_name(c[k]), "ri"[reim]);
                ff.write(dataname, rank, dims, EH + (k*2 + reim)*N*Nfreq);
            }
    }

    delete[] EH;
}

double *dft_near2far::flux(direction df, const volume &where, double resolution) {

    /* compute output grid size etc. */
    int dims[4] = {1,1,1,1};
    double dx[3] = {0,0,0};
    direction dirs[3] = {X,Y,Z};
    int rank = 0, N = 1;
    LOOP_OVER_DIRECTIONS(where.dim, d) {
      dims[rank] = int(floor(where.in_direction(d) * resolution));
      if (dims[rank] <= 1) {
	dims[rank] = 1;
	dx[rank] = 0;
      }
      else
	dx[rank] = where.in_direction(d) / (dims[rank] - 1);
      N *= dims[rank];
      dirs[rank++] = d;
    }

    /* 6 x N x Nfreq array of fields in row-major order */
    std::complex<double> *EH = new std::complex<double>[6*N*Nfreq];
    std::complex<double> *EH_ = new std::complex<double>[6*N*Nfreq]; // temp array for sum_to_master
 
    /* fields for farfield_lowlevel for a single output point x */
    std::complex<double> *EH1 = new std::complex<double>[6*Nfreq];

    vec x(where.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      x.set_direction(dirs[0], where.in_direction_min(dirs[0]) + i0*dx[0]);
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	x.set_direction(dirs[1],
			where.in_direction_min(dirs[1]) + i1*dx[1]);
	for (int i2 = 0; i2 < dims[2]; ++i2) {
	  x.set_direction(dirs[2],
			  where.in_direction_min(dirs[2]) + i2*dx[2]);
	  farfield_lowlevel(EH1, x);
	  int idx = (i0 * dims[1] + i1) * dims[2] + i2;
	  for (int i = 0; i < Nfreq; ++i)
	    for (int k = 0; k < 6; ++k)
	      EH_[(k * N + idx) + (6 * N * i)] = EH1[i * 6 + k];
	}
      }
    }

    delete[] EH1;
    sum_to_master(EH_, EH, 6*N*Nfreq);
    delete[] EH_;

    if (coordinate_mismatch(where.dim, df) || where.dim == Dcyl)
      abort("cannot get flux for near2far: co-ordinate mismatch");

    double *F = new double[Nfreq];
    std::complex<double> *ff_EH[6];
    std::complex<double> *cE[2], *cH[2];
    for (int i = 0; i < Nfreq; ++i) {
      for (int k = 0; k < 6; ++k)
	ff_EH[k] = EH + (k * N) + (6 * N * i);
      switch (df) {
      case X: cE[0] = ff_EH[1], cE[1] = ff_EH[2], cH[0] = ff_EH[5], cH[1] = ff_EH[4]; break;
      case Y: cE[0] = ff_EH[2], cE[1] = ff_EH[0], cH[0] = ff_EH[3], cH[1] = ff_EH[5]; break;
      case Z: cE[0] = ff_EH[0], cE[1] = ff_EH[1], cH[0] = ff_EH[4], cH[1] = ff_EH[3]; break;
      case NO_DIRECTION: abort("cannot get flux in NO_DIRECTION");
      }
      F[i] = 0;
      for (int j = 0; j < 2; ++j) {
        double flux_sum = 0;
	for (int p = 0; p < N; ++p)
	  flux_sum += real(cE[j][p]*conj(cH[j][p]));
	F[i] += flux_sum * (1 - 2*j);
      }
      master_printf("near2far-flux:, %0.2g, %0.5g\n",freq_min+i*dfreq,F[i]);
    }
    delete[] EH;
    return F;
}

static double approxeq(double a, double b) { return fabs(a - b) < 0.5e-11 * (fabs(a) + fabs(b)); }

dft_near2far fields::add_dft_near2far(const volume_list *where,
				double freq_min, double freq_max, int Nfreq){
  dft_chunk *F = 0; /* E and H chunks*/
  double eps = 0, mu = 0;
  
  for (const volume_list *w = where; w; w = w->next) {
      direction nd = component_direction(w->c);
      if (nd == NO_DIRECTION) nd = normal_direction(w->v);
      if (nd == NO_DIRECTION) abort("unknown dft_near2far normal");
      direction fd[2];

      double weps = get_eps(w->v.center());
      double wmu = get_mu(w->v.center());
      if (w != where && !(approxeq(eps, weps) && approxeq(mu, wmu)))
          abort("dft_near2far requires surfaces in a homogeneous medium");
      eps = weps;
      mu = wmu;
      
      /* two transverse directions to normal (in cyclic order to get 
         correct sign s below) */
      switch (nd) {
      case X: fd[0] = Y; fd[1] = Z; break;
      case Y: fd[0] = Z; fd[1] = X; break;
      case R: fd[0] = P; fd[1] = Z; break;
      case P: fd[0] = Z; fd[1] = R; break;
      case Z:
          if (gv.dim == Dcyl)
              fd[0] = R, fd[1] = P;
          else
              fd[0] = X, fd[1] = Y;
          break;
      default: abort("invalid normal direction in dft_near2far!");
      }

      for (int i = 0; i < 2; ++i) { /* E or H */
          for (int j = 0; j < 2; ++j) { /* first or second component */
              component c = direction_component(i == 0 ? Ex : Hx, fd[j]);

              /* find equivalent source component c0 and sign s */
              component c0 = direction_component(i == 0 ? Hx : Ex, fd[1-j]);
              double s = j == 0 ? 1 : -1; /* sign of n x c */
              if (is_electric(c)) s = -s;
              
              F = add_dft(c, w->v, freq_min, freq_max, Nfreq,
                          true, s*w->weight, F, false, 1.0, false, c0);
          }
      }
  }

  return dft_near2far(F, freq_min, freq_max, Nfreq, eps, mu);
}

} // namespace meep

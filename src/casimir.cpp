/* Copyright (C) 2005-2009 Massachusetts Institute of Technology.
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

/* Functions to speed up Casimir-force calculations using FDTD.  It is
   possible to use the unmodified Meep, but if you do this from Scheme
   it is rather slow to perform the stress-tensor integration, and
   also the FFT to compute g(t) is moderately painful in Scheme.  Of
   course, you could just use Meep's C++ interface, but to make it
   more user-friendly we implement the following functions to speed up
   a Scheme front-end implementation of the Casimir calculation. */

#include <math.h>
#include <complex>

#include "meep.hpp"

#include "config.h"

#if defined(HAVE_LIBFFTW)
#  include <fftw.h>
#elif defined(HAVE_LIBFFTW3)
#  include <fftw3.h>
#endif

namespace meep {

typedef complex<double> C;

/* Return an array of values of the g(t) function, for times [0,T] with
   steps dt, for a given Casimir conductivity sigma.  If there is
   any additional frequency dependence of the dielectric function,
   eps_func(omega) should equal eps(omega)/eps(infinity); note that
   the omega argument will be complex.  If Tfft is passed, it is a
   time (should be > T) giving extra resolution for the Fourier transform. */
complex<double> *make_casimir_g(double T, double dt, double sigma, 
				complex<double> (*eps_func)(complex<double> omega),
				double Tfft) {
  if (Tfft <= T) Tfft = T * 10;

  int Nfft = ceil(Tfft / dt);
  C *dg = new C[Nfft];
  
  for (int i = 0; i < Nfft; ++i) dg[i] = 0;
  dg[0] = -sigma;
  for (int i = 1; i < Nfft/2; ++i) {
    double xi = 2*pi*i / (Nfft * dt);
    dg[i] = C(0,-xi) * sqrt(C(1.0, sigma/xi)) * C(1.0, 0.5*sigma/xi)
      - C(sigma, -xi) - 0.5 * sqrt(C(0.0, (sigma*sigma*sigma) / xi));
  }
  if (eps_func)
    for (int i = 1; i < Nfft/2 ; ++i) {
      double xi = 2*pi*i / (Nfft * dt);
      dg[i] = dg[i] * eps_func(xi * sqrt(C(1.0, sigma/xi)));
    }

#if defined(HAVE_LIBFFTW)
  fftw_plan p;
  p = fftw_create_plan(Nfft, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fftw_one(reinterpret_cast<fftw_complex*>(dg),  NULL);
  fftw_destroy_plan(p);
#elif defined(HAVE_LIBFFTW3)
  fftw_plan p;
  p = fftw_plan_dft_1d(Nfft, reinterpret_cast<fftw_complex*>(dg),
		       reinterpret_cast<fftw_complex*>(dg),
		       FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
#else
  abort("make_casimir_g requires some version of FFTW");
#endif

  int N = ceil(T / dt);
  C *g = new C[N];
  g[0] = 0;
  double dxi = 1.0 / (Nfft * dt);
  for (int i = 1; i < N; ++i) {
    double t = i * dt;
    g[i] = dg[i] * dxi + C(0.0,1.0) * ((1/(t*t) + sigma/t) / (2*pi)
				      + 0.25 * sqrt(sigma*sigma*sigma/(t*pi)));
  }

  delete[] dg;
  return g;
}

typedef struct {
  double kx, ky, kz;
  double x0, y0, z0;
  direction xd, yd, zd;
} stress_data;

static complex<double> stress(const complex<double> *fields,
			      const vec &loc,
			      void *data_) {
  stress_data *d = (stress_data *) data_;
  return (fields[0] 
	  * cos(d->kx * (loc.in_direction(d->xd) - d->x0))
	  * cos(d->ky * (loc.in_direction(d->yd) - d->y0))
	  * cos(d->kz * (loc.in_direction(d->zd) - d->z0)));
}

complex<double> fields::casimir_stress_dct_integral(direction dforce,
						    direction dsource,
						    int mx, int my, int mz,
						    field_type ft,
						    geometric_volume where) {
  direction dnormal = where.normal_direction();
  direction dcomponent = NO_DIRECTION;
  double coefficient = 1.0;

  if (where.dim != v.dim)
    abort("invalid dimesionality in casimir_stress_dct_integral");
  if (coordinate_mismatch(v.dim,dforce) || coordinate_mismatch(v.dim,dsource))
    abort("invalid directions in casimir_stress_dct_integral");
  if (dnormal == NO_DIRECTION)
    abort("invalid integration surface in casimir_stress_dct_integral");
  if (ft != E_stuff && ft != H_stuff) 
    abort("invalid field type in casimir_stress_dct_integral");

  if (dforce != dnormal && dsource != dnormal)
    return 0.0;
  else if (dforce != dnormal && dsource == dnormal) {
    // force-source offdiagonal term
    dcomponent = dforce;
  }
  else if (dforce == dnormal && dsource == dnormal) {
    // +source-source/2 diagonal terma
    dcomponent = dsource;
    coefficient = +0.5;
  }
  else /* if (dforce == dnormal && dsource != dnormal) */ {
    // -source-source/2 diagonal term
    dcomponent = dsource;
    coefficient = -0.5;
  }

  component c = direction_component(first_field_component(ft), dcomponent);

  stress_data data;

  data.zd = Z;
  if (v.dim == Dcyl) { data.xd = R; data.yd = P; }
  else { data.xd = X; data.yd = Y; }

  if (has_direction(v.dim, data.xd) && where.in_direction(data.xd) > 0) {
    data.x0 = where.in_direction_min(data.xd);
    data.kx = mx * pi / where.in_direction(data.xd);
    coefficient *= sqrt((mx == 0 ? 1.0 : 2.0) / where.in_direction(data.xd));
  }
  else {
    data.xd = start_at_direction(v.dim); // a dir we are guaranteed to have
    data.x0 = data.kx = 0; // innocuous values: ignore this dir
  }
  if (has_direction(v.dim, data.yd) && where.in_direction(data.yd) > 0) {
    data.y0 = where.in_direction_min(data.yd);
    data.ky = my * pi / where.in_direction(data.yd);
    coefficient *= sqrt((my == 0 ? 1.0 : 2.0) / where.in_direction(data.yd));
  }
  else {
    data.yd = start_at_direction(v.dim); // a dir we are guaranteed to have
    data.y0 = data.ky = 0; // innocuous values: ignore this dir
  }
  if (has_direction(v.dim, data.zd) && where.in_direction(data.zd) > 0) {
    data.z0 = where.in_direction_min(data.zd);
    data.kz = mz * pi / where.in_direction(data.zd);
    coefficient *= sqrt((mz == 0 ? 1.0 : 2.0) / where.in_direction(data.zd));
  }
  else {
    data.zd = start_at_direction(v.dim); // a dir we are guaranteed to have
    data.z0 = data.kz = 0; // innocuous values: ignore this dir
  }
  
  coefficient *= (ft==E_stuff 
		  ? get_eps(where.center()) : get_mu(where.center()));
  
  return coefficient * integrate(1, &c, stress, &data, where);
}

} // namespace meep

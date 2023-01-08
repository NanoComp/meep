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

/* Anisotropic dispersion test program */

/* For comparison, we solve for the dispersion relation
   analytically in Matlab using the fsolve command to solve
   the nonlinear eigenproblem.  The Matlab code is as follows.

function result = detMwk(f,k)
  sigE = [ 2.92724 0.45948 0.70117;
          0.45948 2.89689 0.45083;
          0.70117 0.45083 2.17378; ];
  epsinf = [ 2.41104 0.48709 0.41226;
             0.48709 2.43172 1.62060;
             0.41226 1.62060 3.61498; ];
  sigH = 0;
  muinf = [ 1 0 0; 0 1 0; 0 0 1; ];
  f0 = 1.1; g0 = 1e-5;
  epsinv = inv(epsinf + ((f0^2)/(f0^2 - f^2 -i*f*g0))*sigE);
  muinv = inv(muinf + ((f0^2)/(f0^2 - f^2 -i*f*g0))*sigH);
  E = [ 0 0; 1 0; 0 1; ];
  kx = k * [ 0 0 0; 0 0 -1; 0 1 0; ];
  result = det(E' * (f.^2*eye(3) + kx * muinv * kx * epsinv) * E);
endfunction

k = 0.813;
w1 = fsolve(@(w)detMwk(w,k),0.9);
w2 = fsolve(@(w)detMwk(w,k),w1-0.1);

*/

#include <vector>
#include <meep.hpp>
using namespace meep;
using std::complex;

class anisodisp_material : public material_function {
public:
  virtual void eff_chi1inv_row(component c, double chi1inv_row[3], const volume &v,
                               double tol = DEFAULT_SUBPIXEL_TOL,
                               int maxeval = DEFAULT_SUBPIXEL_MAXEVAL) {
    (void)v;
    (void)tol;
    (void)maxeval; // unused
    if (component_direction(c) == X) {
      chi1inv_row[0] = 0.432818;
      chi1inv_row[1] = -0.076724;
      chi1inv_row[2] = -0.014964;
    }
    else if (component_direction(c) == Y) {
      chi1inv_row[0] = -0.076724;
      chi1inv_row[1] = 0.600041;
      chi1inv_row[2] = -0.260249;
    }
    else {
      chi1inv_row[0] = -0.014964;
      chi1inv_row[1] = -0.260249;
      chi1inv_row[2] = 0.395003;
    }
  }
  virtual void sigma_row(component c, double sigrow[3], const vec &r) {
    (void)r; // unused
    if (component_direction(c) == X) {
      sigrow[0] = 2.92724;
      sigrow[1] = 0.45948;
      sigrow[2] = 0.70117;
    }
    else if (component_direction(c) == Y) {
      sigrow[0] = 0.45948;
      sigrow[1] = 2.89689;
      sigrow[2] = 0.45083;
    }
    else {
      sigrow[0] = 0.70117;
      sigrow[1] = 0.45083;
      sigrow[2] = 2.17378;
    }
  }
};

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  bool ok = true;
  // we can only use one process for this 1-pixel simulation
  if (0 == divide_parallel_processes(count_processors())) {
    verbosity = 0;
    const double res = 200;
    grid_volume gv = vol3d(0, 0, 0, res);
    gv.center_origin();
    anisodisp_material anisodispmat;
    structure *s = new structure(gv, anisodispmat);
    s->add_susceptibility(anisodispmat, E_stuff, lorentzian_susceptibility(1.1, 1e-5));
    fields f(s);
    delete s; // should be safe since structure_chunk in f is refcounted
    f.use_bloch(vec(0.813, 0, 0));
    f.add_point_source(Ez, 0.5, 1.0, 0.0, 4.0, vec(0, 0, 0));
    double T = f.last_source_time();
    int iT = T / f.dt;
    while (f.t < iT) {
      if (f.t % (iT / 10) == 0) master_printf("%g%% done with source\n", f.time() / T * 100);
      f.step();
    }
    double T2 = 200;
    int iT2 = T2 / f.dt;
    std::vector<complex<double> > vals(iT2);
    while (f.t - iT < iT2) {
      if ((f.t - iT) % (iT2 / 10) == 0)
        master_printf("%g%% done with harminv\n", (f.t - iT) * 100.0 / iT2);
      vals[f.t - iT] = f.get_field(Ez, vec(0., 0., 0.));
      f.step();
    }
    complex<double> amps[8];
    double freqs_re[8], freqs_im[8];

    master_printf("done with timestepping, running harminv...\n");
    int num = do_harminv(vals.data(), iT2, f.dt, 0.0, 1.0, 8, amps, freqs_re, freqs_im);

    // compute the error compared to analytical solution
    int i0 = 0;
    for (int i = 0; i < num; i++) {
      master_printf("freq %d is %0.6g, %0.6g\n", i, freqs_re[i], freqs_im[i]);
      if (fabs(freqs_re[i] - 0.41562) < fabs(freqs_re[i0] - 0.41562)) i0 = i;
    }
    master_printf("err. real: %g\n", fabs(freqs_re[i0] - 0.41562) / 0.41562);
    master_printf("err. imag: %g\n", fabs(freqs_im[i0] + 4.8297e-07) / 4.8297e-7);

    double tol = sizeof(realnum) == sizeof(float) ? 0.27 : 0.20;
    ok = fabs(freqs_re[i0] - 0.41562) / 0.41562 < 1e-4 &&
         fabs(freqs_im[i0] + 4.8297e-07) / 4.8297e-7 < tol;
  }
  end_divide_parallel();
  return !and_to_all(ok);
}

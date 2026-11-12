/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
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

/* Check of Green's functions (analytical vs. numerical)
   and near-to-far-field transformation. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <meep.hpp>
#include "meep_internals.hpp"
#include "config.h"
using namespace meep;
using std::complex;

static bool check(double x, double xref, const char *what, int m, int n, double tol = 1e-12) {
  if (fabs(x - xref) > tol * fabs(xref)) {
    master_printf("%s(%d, %d) = %.17g, expected %.17g\n", what, m, n, x, xref);
    return false;
  }
  return true;
}

// reference values from scipy.special.jn_zeros and jnp_zeros
static const int ms[] = {0, 1, 2, 5, 10, 30};
static const double jn_zeros[][4] = {
    {2.4048255576957724, 5.520078110286311, 8.653727912911013, 11.791534439014281},
    {3.8317059702075125, 7.015586669815619, 10.173468135062722, 13.323691936314223},
    {5.135622301840683, 8.417244140399866, 11.61984117214906, 14.795951782351262},
    {8.771483815959954, 12.338604197466944, 15.70017407971167, 18.98013387517992},
    {14.47550068655454, 18.43346366696658, 22.0469853646978, 25.509450554182827},
    {36.09833695674772, 41.092778663153425, 45.45266628755661, 49.506177061630915}};
// positive zeros of J_m' (excluding x = 0)
static const double jnp_zeros[][3] = {{3.8317059702075125, 7.015586669815619, 10.173468135062722},
                                      {1.8411837813406595, 5.3314427735250325, 8.536316366346286},
                                      {3.0542369282271404, 6.706133194158459, 9.969467823087596},
                                      {6.415616375700241, 10.519860873772307, 13.9871886301403},
                                      {11.77087667495558, 16.447852748486497, 20.223031412681703},
                                      {32.534223556790145, 38.63609270136744, 43.28070570283707}};

bool test_zeros() {
  bool ok = true;
  for (int i = 0; i < 6; ++i) {
    const int m = ms[i];
    for (int n = 0; n < 4; ++n) {
      ok = check(Jroot(m, n), jn_zeros[i][n], "Jroot", m, n) && ok;
      if (fabs(BesselJ(m, Jroot(m, n))) > 1e-14) {
        master_printf("BesselJ(%d, Jroot(%d, %d)) = %g is not zero\n", m, m, n,
                      BesselJ(m, Jroot(m, n)));
        ok = false;
      }
    }
    // Jmax(m, n) is the zero of J_m' between Jroot(m, n-1) and Jroot(m, n), so
    // for m = 0 the zero of J_0' at x = 0 corresponds to n = 0 and is skipped
    for (int n = 0; n < 3; ++n) {
      const int nmax = m == 0 ? n + 1 : n;
      ok = check(Jmax(m, nmax), jnp_zeros[i][n], "Jmax", m, nmax) && ok;
    }
  }
  return ok;
}

bool test_hankel() {
  // reference values of the Hankel function H^{(1)}_n(x) from scipy.special.hankel1
  struct {
    int n;
    double x, re, im;
  } ref[] = {{0, 0.5, 0.9384698072408126, -0.44451873350670656},
             {0, 3.0, -0.2600519549019336, 0.37685001001279056},
             {0, 25.0, 0.09626678327595814, -0.12724943226800617},
             {1, 0.5, 0.24226845767487393, -1.4714723926702433},
             {1, 3.0, 0.33905895852593654, 0.32467442479180014},
             {1, 25.0, -0.12535024958028995, -0.09882996478323747},
             {2, 0.5, 0.03060402345868268, -5.441370837174266},
             {2, 3.0, 0.48609126058589125, -0.1604003934849238},
             {2, 25.0, -0.10629480324238133, 0.11934303508534717}};
  bool ok = true;
  for (const auto &r : ref) {
    complex<double> h = hankel(r.n, r.x), href(r.re, r.im);
    if (abs(h - href) > 1e-12 * abs(href)) {
      master_printf("hankel(%d, %g) = %.17g%+.17gi, expected %.17g%+.17gi\n", r.n, r.x, real(h),
                    imag(h), r.re, r.im);
      ok = false;
    }
  }
  return ok;
}

static double one(const vec &) { return 1.0; }

// check that initialize_with_nth_tm/te produce J_m(k r) profiles with the
// correct transverse wavevector k for a cylinder of radius R, including when
// the cell is split into several chunks along r
// (uses the reference zeros above, whose first rows are m = 0, 1, 2)
bool test_initialize(int num_chunks) {
  const double R = 1.0, a = 40.0, tol = 1e-2;
  grid_volume gv = volcyl(R, 0.5, a);
  // qualified to avoid ambiguity with C++20's std::identity, since
  // meep_internals.hpp has "using namespace std"
  structure s(gv, one, no_pml(), meep::identity(), num_chunks);
  bool ok = true;
  for (int m = 0; m < 3; ++m) {
    for (int np = 1; np <= 2; ++np) {
      fields ftm(&s, m);
      ftm.use_bloch(0.0);
      ftm.initialize_with_nth_tm(np);
      const double ktm = jn_zeros[m][np - 1] / R;
      fields fte(&s, m);
      fte.use_bloch(0.0);
      fte.initialize_with_nth_te(np);
      const double kte = jnp_zeros[m][np - 1] / R;
      for (double r = 0.1; r < R; r += 0.2) {
        const complex<double> ez = ftm.get_field(Ez, veccyl(r, 0.25));
        const complex<double> hz = fte.get_field(Hz, veccyl(r, 0.25));
        if (abs(ez - BesselJ(m, ktm * r)) > tol) {
          master_printf("TM m=%d n=%d (%d chunks): Ez(r=%g) = %g%+gi, expected %g\n", m, np,
                        num_chunks, r, real(ez), imag(ez), BesselJ(m, ktm * r));
          ok = false;
        }
        if (abs(hz - BesselJ(m, kte * r)) > tol) {
          master_printf("TE m=%d n=%d (%d chunks): Hz(r=%g) = %g%+gi, expected %g\n", m, np,
                        num_chunks, r, real(hz), imag(hz), BesselJ(m, kte * r));
          ok = false;
        }
      }
    }
  }
  return ok;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;

#ifndef HAVE_JN
  master_printf("not compiled with the jn and yn functions; skipping Bessel tests\n");
  return 77; // tells automake that the test was skipped
#endif

  if (!test_zeros()) meep::abort("error in test_zeros\n");
  if (!test_hankel()) meep::abort("error in test_hankel\n");
  if (!test_initialize(1)) meep::abort("error in test_initialize with 1 chunk\n");
  if (!test_initialize(3)) meep::abort("error in test_initialize with 3 chunks\n");

  return 0;
}

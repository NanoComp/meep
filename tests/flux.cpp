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

#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
using namespace meep;

double one(const vec &) { return 1.0; }
static double width = 20.0;
double bump(const vec &pt) { return (fabs(pt.z() - 50.0) > width) ? 1.0 : 12.0; }

double bump2(const vec &pt) { return (fabs(pt.z() - 5.0) > 3.0) ? 1.0 : 12.0; }

double cavity(const vec &pt) {
  const double zz = fabs(pt.z() - 7.5) + 0.3001;
  if (zz > 5.0) return 1.0;
  if (zz < 2.0) return 1.0;
  double norm = zz;
  while (norm > 1.0)
    norm -= 1.0;
  if (norm > 0.3) return 1.0;
  return 12.0;
}

int compare(double a, double b, double eps, double thresh, const char *n) {
  if (fabs(a - b) > fabs(b) * eps && fabs(b) > thresh) {
    master_printf("%s differs by\t%g out of\t%g\n", n, a - b, b);
    master_printf("This gives a fractional error of %g\n", fabs(a - b) / fabs(b));
    return 0;
  }
  else {
    if (fabs(a - b) > fabs(b) * eps * 1.1 && fabs(b) > thresh)
      master_printf("%s fractional error is %g, close to %g threshold.\n", n, fabs(a - b) / fabs(b),
                    eps);
    return 1;
  }
}

static inline double min(double a, double b) { return (a < b) ? a : b; }

int flux_1d(const double zmax, double eps(const vec &)) {
  const double a = 10.0;

  grid_volume gv = volone(zmax, a);
  structure s(gv, eps, pml(zmax / 6));

  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, 0.25, 3.5, 0.0, 8.0, vec(zmax / 6 + 0.3), 1.0);
  flux_vol *left = f.add_flux_plane(vec(zmax / 3.0), vec(zmax / 3.0));
  flux_vol *right = f.add_flux_plane(vec(zmax * 2.0 / 3.0), vec(zmax * 2.0 / 3.0));

  const double ttot = min(10.0 + 1e5 / zmax, f.last_source_time());

  f.step();
  grid_volume mid = volone(zmax / 3, a);
  mid.set_origin(vec(zmax / 3));
  double flux_left = 0.0, flux_right = 0.0;
  double delta_energy = f.field_energy_in_box(mid.surroundings());
  master_printf("Initial energy is %g\n", f.field_energy_in_box(mid.surroundings()));
  master_printf("Initial electric energy is %g\n", f.electric_energy_in_box(mid.surroundings()));
  while (f.time() < ttot) {
    f.step();
    flux_left += f.dt * left->flux();
    flux_right += f.dt * right->flux();
  }
  delta_energy -= f.field_energy_in_box(mid.surroundings());
  master_printf("Final energy is %g\n", f.field_energy_in_box(mid.surroundings()));
  master_printf("Final electric energy is %g\n", f.electric_energy_in_box(mid.surroundings()));
  const double del = flux_left;
  const double der = flux_right - delta_energy;
  master_printf("  Delta E:\t%g\n  Flux left:\t%g\n  Flux right:\t%g\n  Ratio:\t%g\n", delta_energy,
                del, der, del / der);
  return compare(del, der, 0.06, 0, "Flux");
}

int split_1d(double eps(const vec &), int splitting) {
  const double boxwidth = 5.0, timewait = 1.0;
  const double zmax = 15.0, a = 10.0;

  grid_volume gv = volone(zmax, a);
  structure s1(gv, eps, pml(2.0));
  structure s(gv, eps, pml(2.0), identity(), splitting);

  fields f1(&s1);
  fields f(&s);
  f1.use_real_fields();
  f.use_real_fields();
  f1.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax / 2 + 0.3), 1.0e2);
  f.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax / 2 + 0.3), 1.0e2);
  flux_vol *left1 = f1.add_flux_plane(vec(zmax * .5 - boxwidth), vec(zmax * .5 - boxwidth));
  flux_vol *left = f.add_flux_plane(vec(zmax * .5 - boxwidth), vec(zmax * .5 - boxwidth));
  grid_volume mid = volone(2 * boxwidth, a);
  mid.set_origin(vec(zmax * .5 - boxwidth - 0.25 / a));

  const double ttot = f.last_source_time() + timewait;
  const double tol = sizeof(realnum) == sizeof(float) ? 1e-3 : 1e-9;
  while (f.time() < ttot) {
    f1.step();
    f.step();
    if (!compare(f.dt * left1->flux(), f.dt * left->flux(), tol, tol, "Flux")) return 0;
  }
  return 1;
}

int cavity_1d(const double boxwidth, const double timewait, double eps(const vec &)) {
  const double zmax = 15.0;
  const double a = 10.0;

  grid_volume gv = volone(zmax, a);
  structure s(gv, eps, pml(2.0));

  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, 0.25, 4.5, 0.0, 8.0, vec(zmax / 2 + 0.3), 1.0e2);
  flux_vol *left = f.add_flux_plane(vec(zmax * .5 - boxwidth), vec(zmax * .5 - boxwidth));
  flux_vol *right = f.add_flux_plane(vec(zmax * .5 + boxwidth), vec(zmax * .5 + boxwidth));
  grid_volume mid = volone(2 * boxwidth, a);
  mid.set_origin(vec(zmax * .5 - boxwidth - 0.25 / a));

  while (f.time() < f.last_source_time())
    f.step();
  const double ttot = f.time() + timewait;
  double flux_left = 0.0, flux_right = 0.0;
  const double start_energy = f.field_energy_in_box(mid.surroundings());
  master_printf("  Energy starts at\t%g\n", start_energy);
  while (f.time() < ttot) {
    f.step();
    flux_left += f.dt * left->flux();
    flux_right += f.dt * right->flux();
  }
  const double delta_energy = start_energy - f.field_energy_in_box(mid.surroundings());
  const double defl = flux_right - flux_left;
  master_printf("  Delta E:         \t%g\n  Integrated Flux:\t%g\n", delta_energy, defl);
  master_printf("  Ratio:         \t%g\n", delta_energy / defl);
  master_printf("  Fractional error:\t%g\n", (delta_energy - defl) / start_energy);
  return compare(start_energy - delta_energy, start_energy - defl, (timewait > 50) ? 0.032 : 0.004,
                 0, "Flux");
}

int flux_2d(const double xmax, const double ymax, double eps(const vec &)) {
  const double a = 8.0;

  master_printf("\nFlux_2d(%g,%g) test...\n", xmax, ymax);

  grid_volume gv = voltwo(xmax, ymax, a);
  structure s(gv, eps, pml(0.5));

  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ez, 0.25, 3.5, 0., 8., vec(xmax / 6 + 0.1, ymax / 6 + 0.3), 1.);

  // corners of flux planes and energy box:
  vec lb(vec(xmax / 3, ymax / 3)), rb(vec(2 * xmax / 3, ymax / 3));
  vec lt(vec(xmax / 3, 2 * ymax / 3)), rt(vec(2 * xmax / 3, 2 * ymax / 3));
  volume box(lb, rt);

  flux_vol *left = f.add_flux_plane(lb, lt);
  flux_vol *right = f.add_flux_plane(rb, rt);
  flux_vol *bottom = f.add_flux_plane(lb, rb);
  flux_vol *top = f.add_flux_plane(lt, rt);

  /* measure flux spectra through two concentric flux boxes
     around the source...should be positive and equal */
  volume box1(vec(xmax / 6 - 0.4, ymax / 6 - 0.2), vec(xmax / 6 + 0.6, ymax / 6 + 0.8));
  volume box2(vec(xmax / 6 - 0.9, ymax / 6 - 0.7), vec(xmax / 6 + 1.1, ymax / 6 + 1.3));
  int Nfreq = 10;
  double freq_array[] = {0.230, 0.232, 0.238, 0.241, 0.248, 0.254, 0.256, 0.265, 0.269, 0.270};
  // workaround for C++98 which does not support list initialization
  const std::vector<double> freq(freq_array, freq_array + sizeof(freq_array) / sizeof(double));
  dft_flux flux1 = f.add_dft_flux_box(box1, freq);
  dft_flux flux2 = f.add_dft_flux_box(box2, freq);

  const double ttot = 130;

  /* first check: integral of flux = change in energy of box */
  f.step();
  double init_energy = f.field_energy_in_box(box);
  master_printf("Initial energy is %g\n", init_energy);
  long double fluxL = 0;
  while (f.time() < ttot) {
    f.step();
    fluxL += f.dt * (left->flux() - right->flux() + bottom->flux() - top->flux());
  }
  double flux = fluxL;
  double del_energy = f.field_energy_in_box(box) - init_energy;
  master_printf("Final energy is %g\n", f.field_energy_in_box(box));
  master_printf("  delta E: %g\n  net flux: %g\n  ratio: %g\n", del_energy, flux,
                del_energy / flux);
  if (!compare(del_energy, flux, 0.09, 0, "Flux")) return 0;

  /* second check: flux spectrum is same for two concentric
     boxes containing the source. */
  while (f.time() < ttot * 2) {
    f.step();
  }
  master_printf("  energy after more time is %g\n", f.field_energy_in_box(box));
  master_printf("  and energy in box2 is %g\n", f.field_energy_in_box(box2));
  double *fl1 = flux1.flux();
  double *fl2 = flux2.flux();
  for (int i = 0; i < Nfreq; ++i) {
    master_printf("  flux(%g) = %g vs. %g (rat. = %g)\n", flux1.freq[i], fl1[i], fl2[i],
                  fl1[i] / fl2[i]);
    if (!compare(fl1[i], fl2[i], 0.09, 0, "Flux spectrum")) return 0;
  }
  delete[] fl2;
  delete[] fl1;

  return 1;
}

int flux_cyl(const double rmax, const double zmax, double eps(const vec &), int m) {
  const double a = 8.0;

  master_printf("\nFlux_cyl(%g,%g) test...\n", rmax, zmax);

  grid_volume gv = volcyl(rmax, zmax, a);
  structure s(gv, eps, pml(0.5), identity(), 0, min(0.5, 1 / (abs(m) + 0.5)));

  fields f(&s, m);
  // f.use_real_fields();
  f.add_point_source(Ep, 0.25, 3.5, 0., 8., veccyl(rmax * 5 / 6 + 0.1, zmax / 6 + 0.3), 1.);

  // corners of flux planes and energy box:
  vec lb(veccyl(-rmax / 3, zmax / 3)), rb(veccyl(2 * rmax / 3, zmax / 3));
  vec lt(veccyl(-rmax / 3, 2 * zmax / 3)), rt(veccyl(2 * rmax / 3, 2 * zmax / 3));
  volume box(lb, rt);

  /* measure flux spectra through two concentric flux boxes
     around the source...should be positive and equal */
  volume box1(veccyl(rmax * 5 / 6 - 0.4, zmax / 6 - 0.2),
              veccyl(rmax * 5 / 6 + 0.6, zmax / 6 + 0.8));
  volume box2(veccyl(rmax * 5 / 6 - 0.9, zmax / 6 - 0.7),
              veccyl(rmax * 5 / 6 + 1.1, zmax / 6 + 1.3));
  double fmin = 0.23, fmax = 0.27;
  int Nfreq = 10;
  dft_flux flux1 = f.add_dft_flux_box(box1, fmin, fmax, Nfreq);
  dft_flux flux2 = f.add_dft_flux_box(box2, fmin, fmax, Nfreq);

  flux_vol *left = f.add_flux_plane(lb, lt);
  flux_vol *right = f.add_flux_plane(rb, rt);
  flux_vol *bottom = f.add_flux_plane(lb, rb);
  flux_vol *top = f.add_flux_plane(lt, rt);

  const double ttot = 130;

  f.step();
  double init_energy = f.field_energy_in_box(box);
  master_printf("Initial energy is %g\n", init_energy);
  long double fluxL = 0;
  while (f.time() < ttot) {
    f.step();
    fluxL += f.dt * (left->flux() - right->flux() + bottom->flux() - top->flux());
  }
  double flux = fluxL;
  double del_energy = f.field_energy_in_box(box) - init_energy;
  master_printf("Final energy is %g\n", f.field_energy_in_box(box));
  master_printf("  delta E: %g\n  net flux: %g\n  ratio: %g\n", del_energy, flux,
                del_energy / flux);
  if (!compare(del_energy, flux, 0.08, 0, "Flux")) return 0;

  while (f.time() < ttot * 2) {
    f.step();
  }
  master_printf("  energy after more time is %g\n", f.field_energy_in_box(box));
  master_printf("  and energy in box2 is %g\n", f.field_energy_in_box(box2));
  double *fl1 = flux1.flux();
  double *fl2 = flux2.flux();
  for (int i = 0; i < Nfreq; ++i) {
    master_printf("  flux(%g) = %g vs. %g (rat. = %g)\n", flux1.freq[i], fl1[i], fl2[i],
                  fl1[i] / fl2[i]);
    if (!compare(fl1[i], fl2[i], 0.08, 0, "Flux spectrum")) return 0;
  }
  delete[] fl2;
  delete[] fl1;

  return 1;
}

void attempt(const char *name, int allright) {
  if (allright)
    master_printf("Passed %s\n", name);
  else
    meep::abort("Failed %s!\n", name);
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Trying out the fluxes...\n");

  attempt("Split flux plane split by 7...", split_1d(cavity, 7));

  attempt("Cavity 1D 1.3 73", cavity_1d(1.3, 73.0, cavity));
  attempt("Cavity 1D 5.0   1", cavity_1d(5.0, 1.0, cavity));
  attempt("Cavity 1D 3.85 55", cavity_1d(3.85, 55.0, cavity));

  width = 20.0;
  attempt("Flux 1D 20", flux_1d(100.0, bump));
  width = 10.0;
  attempt("Flux 1D 10", flux_1d(100.0, bump));
  width = 300.0;
  attempt("Flux 1D 300", flux_1d(100, bump));

  width = 5.0;
  attempt("Flux 2D 5", flux_2d(10.0, 10.0, bump2));

  width = 5.0;
  attempt("Flux cylindrical 5", flux_cyl(20.0, 10.0, bump2, 1));

  return 0;
}

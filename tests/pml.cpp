#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
using namespace meep;
using std::complex;
using std::max;
using std::polar;

double Rasymp = 1e-15, stretch = 2.0; // PML parameters

// a simple material with xy offdiagonal terms in the tensors, for testing
class offdiag_material : public material_function {
public:
  offdiag_material(double offdiag) : offdiag(offdiag) {}
  virtual bool has_mu() { return true; }
  virtual void eff_chi1inv_row(component c, double chi1inv_row[3], const volume &v,
                               double tol = DEFAULT_SUBPIXEL_TOL,
                               int maxeval = DEFAULT_SUBPIXEL_MAXEVAL) {
    (void)v;
    (void)tol;
    (void)maxeval; // unused
    // we are returning inv(1+chi1), so we must compute the inverse
    // inv([1+od od; od 1+od]) = [1+od -od; -od 1+od] / (1+2*o)
    double detinv = 1.0 / (1 + 2 * offdiag);
    if (component_direction(c) == X) {
      chi1inv_row[0] = (1 + offdiag) * detinv;
      chi1inv_row[1] = -offdiag * detinv;
      chi1inv_row[2] = 0.0;
    }
    else if (component_direction(c) == Y) {
      chi1inv_row[0] = -offdiag * detinv;
      chi1inv_row[1] = (1 + offdiag) * detinv;
      chi1inv_row[2] = 0.0;
    }
    else {
      chi1inv_row[0] = 0.0;
      chi1inv_row[1] = 0.0;
      chi1inv_row[2] = 1.0;
    }
  }
  double offdiag;
};

static double one(const vec &) { return 1.0; }
static double notone_val = 1.0;
static double notone(const vec &) { return notone_val; }

static complex<double> do_ft(fields &f, component c, const vec &pt, double freq) {
  complex<double> ft = 0.0;
  double emax = 0;
  while (f.time() < f.last_source_time()) {
    complex<double> fpt = f.get_field(c, pt);
    ft += fpt * polar(1.0, 2 * pi * freq * f.time());
    emax = max(emax, abs(fpt));
    f.step();
  }
  do {
    double emaxcur = 0;
    double T = f.time() + 50;
    while (f.time() < T) {
      complex<double> fpt = f.get_field(c, pt);
      ft += fpt * polar(1.0, 2 * pi * freq * f.time());
      double e = abs(fpt);
      emax = max(emax, e);
      emaxcur = max(emaxcur, e);
      f.step();
    }
    if (emaxcur < (sizeof(realnum) == sizeof(float) ? 1e-4 : 1e-6) * emax) break;
    if (T > 500 && emaxcur > 1e-2 * emax) meep::abort("fields do not seem to be decaying");
  } while (1);
  return ft;
}

int check_pml1d(double eps(const vec &), double conductivity) {
  double freq = 1.0, dpml = 1.0;
  double sz = 1.0 + 2 * dpml;
  double sz2 = 1.0 + 2 * dpml * 2;
  complex<double> ft = 0.0, ft2 = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  vec fpt(0.5 * sz - dpml - 0.1);
  master_printf("Checking resolution convergence of 1d PML...\n");
  if (conductivity != 0) master_printf("...with conductivity %g...\n", conductivity);
  notone_val = conductivity;
  for (int i = 0; i < (sizeof(realnum) == sizeof(float) ? 5 : 8); i++) {
    double res = 10.0 + 10.0 * i;
    {
      grid_volume gv = vol1d(sz, res);
      gv.center_origin();
      structure s(gv, eps, pml(dpml, Rasymp, stretch));
      s.set_conductivity(By, notone);
      fields f(&s);
      gaussian_src_time src(freq, freq / 20);
      f.add_point_source(Ex, src, vec(-0.5 * sz + dpml + 0.1));
      ft = do_ft(f, Ex, fpt, freq);
    }
    {
      grid_volume gv = vol1d(sz2, res);
      gv.center_origin();
      structure s(gv, eps, pml(dpml * 2, Rasymp, stretch));
      s.set_conductivity(By, notone);
      fields f(&s);
      gaussian_src_time src(freq, freq / 20);
      f.add_point_source(Ex, src, vec(-0.5 * sz + dpml + 0.1));
      ft2 = do_ft(f, Ex, fpt, freq);
    }
    refl_const = pow(abs(ft - ft2), 2.0) / pow(abs(ft2), 2.0);
    master_printf("refl1d:, %g, %g\n", res, refl_const);
    if (i > 0 && refl_const > prev_refl_const * pow((res - 10) / res, 8.0) * 1.1) return 1;
    prev_refl_const = refl_const;
  }
  master_printf("passed 1d PML check.\n");
  return 0;
}

int check_pml2d(double eps(const vec &), component c, double conductivity, bool dispersion,
                double offdiag) {
  double freq = 1.0, dpml = 1.0, sigma0 = 1.0, omega0 = 1.0, gamma0 = 0.3;
  complex<double> ft = 0.0, ft2 = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  double sxy = 4.0 + 2 * dpml;
  double sxy2 = 4.0 + 2 * dpml * 2;
  double res_step = 6.0;
  vec fpt(0.5 * sxy - dpml - 0.1, 0);
  if (c != Ez && c != Hz) meep::abort("unimplemented component check");
  double symsign = c == Ez ? 1.0 : -1.0;
  master_printf("Checking resolution convergence of 2d %s PML...\n", c == Ez ? "TM" : "TE");
  if (conductivity != 0) master_printf("...with conductivity %g...\n", conductivity);
  if (dispersion) master_printf("...with dispersion\n");
  if (offdiag != 0) master_printf("...with offdiag %g...\n", offdiag);
  offdiag_material mat(offdiag);
  for (int i = 0; i < 4; i++) {
    double res = 10.0 + res_step * i;
    {
      grid_volume gv = vol2d(sxy, sxy, res);
      gv.center_origin();
      const symmetry S =
          offdiag != 0 ? rotate2(Z, gv) : mirror(X, gv) * symsign + mirror(Y, gv) * symsign;
      structure s(gv, eps, pml(dpml, Rasymp, stretch), S);
      if (conductivity != 0) {
        notone_val = conductivity;
        s.set_conductivity(Bx, notone);
        s.set_conductivity(By, notone);
        s.set_conductivity(Bz, notone);
        s.set_conductivity(Dx, notone);
        s.set_conductivity(Dy, notone);
        s.set_conductivity(Dz, notone);
      }
      if (dispersion) {
        notone_val = sigma0;
        s.add_susceptibility(notone, E_stuff, lorentzian_susceptibility(omega0, gamma0));
      }
      if (offdiag != 0) s.set_materials(mat, false);
      fields f(&s);
      f.use_real_fields();
      gaussian_src_time src(freq, freq / 20);
      f.add_point_source(c, src, gv.center());
      ft = do_ft(f, c, fpt, freq);
    }
    {
      grid_volume gv = vol2d(sxy2, sxy2, res);
      gv.center_origin();
      const symmetry S =
          offdiag != 0 ? rotate2(Z, gv) : mirror(X, gv) * symsign + mirror(Y, gv) * symsign;
      structure s(gv, eps, pml(dpml * 2, Rasymp, stretch), S);
      if (conductivity != 0) {
        notone_val = conductivity;
        s.set_conductivity(Bx, notone);
        s.set_conductivity(By, notone);
        s.set_conductivity(Bz, notone);
        s.set_conductivity(Dx, notone);
        s.set_conductivity(Dy, notone);
        s.set_conductivity(Dz, notone);
      }
      if (dispersion) {
        notone_val = sigma0;
        s.add_susceptibility(notone, E_stuff, lorentzian_susceptibility(omega0, gamma0));
      }
      if (offdiag != 0) s.set_materials(mat, false);
      fields f(&s);
      f.use_real_fields();
      gaussian_src_time src(freq, freq / 20);
      f.add_point_source(c, src, gv.center());
      ft2 = do_ft(f, c, fpt, freq);
    }
    refl_const = pow(abs(ft - ft2), 2.0) / pow(abs(ft2), 2.0);
    master_printf("refl2d:, %g, %g\n", res, refl_const);
    if (i > 0 &&
        refl_const > prev_refl_const * pow((res - res_step) / res, offdiag != 0 ? 6.0 : 8.0) * 1.2)
      return 1;
    prev_refl_const = refl_const;
  }
  master_printf("passed 2d %s PML check.\n", c == Ez ? "TM" : "TE");
  return 0;
}

/* The cylindrical case actually shouldn't have a reflection that goes
   to zero with increasing resolution - we implement only a
   "quasi-PML" for cylindrical coordinates, which is basically the PML
   for Cartesian coordinates copied over directly to the cylindrical
   case, rather than doing a proper coordinate stretching of r.  This
   is not a practical issue because, rather than increasing the
   resolution, in practice you increase the PML thickness to eliminate
   reflections, and increasing a quasi-PML thickness makes the
   reflection vanish by the usual adiabatic theorem.

   Because of that, we don't actually run this check as part of the
   test suite, but I'll leave the code here for future study of the
   cylindrical PML. */
int check_pmlcyl(double eps(const vec &)) {
  double freq = 1.0, dpml = 1.0;
  complex<double> ft = 0.0, ft2 = 0.0;
  double refl_const = 0.0;
  double sr = 5.0 + dpml, sz = 1.0 + 2 * dpml;
  double sr2 = 5.0 + dpml * 2, sz2 = 1.0 + 2 * dpml * 2;
  vec fpt = veccyl(sr - dpml - 0.1, 0);
  double res_step = 6.0;
  master_printf("Checking resolution convergence of cylindrical PML...\n");
  for (int i = 0; i < 5; i++) {
    double res = 10.0 + res_step * i;
    master_printf("    checking cylindrical resolution %g...\n", res);
    {
      grid_volume gv = volcyl(sr, sz, res);
      gv.center_origin();
      structure s(gv, eps, pml(dpml, Rasymp, stretch));
      fields f(&s, 0);
      gaussian_src_time src(freq, freq / 20);
      f.add_point_source(Ez, src, veccyl(0.1, 0.1));
      ft = do_ft(f, Ez, fpt, freq);
    }
    {
      grid_volume gv = volcyl(sr2, sz2, res);
      gv.center_origin();
      structure s(gv, eps, pml(dpml * 2, Rasymp, stretch));
      fields f(&s, 0);
      gaussian_src_time src(freq, freq / 20);
      f.add_point_source(Ez, src, veccyl(0.1, 0.1));
      ft2 = do_ft(f, Ez, fpt, freq);
    }
    refl_const = pow(abs(ft - ft2), 2.0) / pow(abs(ft2), 2.0);
    master_printf("reflcyl:, %g, %g\n", res, refl_const);
  }
  master_printf("passed cylindrical PML check.\n");
  return 0;
}

int pml1d_scaling(double eps(const vec &)) {
  double res = 20, freq = 1.0, dpml = 0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double refl_const = 0.0, prev_refl_const = 0.0;
  master_printf("Checking thickness convergence of 1d PML...\n");
  for (int i = 0; i < (sizeof(realnum) == sizeof(float) ? 5 : 7); i++) {
    dpml = pow(2.0, (double)i);
    double sz = 2 * dpml + 10.0 + dpml;
    prev_ft = ft;

    grid_volume gv = vol1d(sz, res);
    structure s(gv, eps,
                (pml(2 * dpml, Z, Low, Rasymp, stretch) + pml(dpml, Z, High, Rasymp, stretch)) *
                    1.5);
    fields f(&s);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ex, src, vec(2 * dpml + 0.1));
    ft = do_ft(f, Ex, vec(sz - dpml - 0.1), freq);

    if (i > 0) {
      refl_const = pow(abs(ft - prev_ft), 2.0) / pow(abs(prev_ft), 2.0);
      master_printf("refl1d:, %g, %g\n", dpml, refl_const);
      if (refl_const > (1e-9) * pow(2 / dpml, 6.0) || refl_const < (1e-10) * pow(2 / dpml, 6.0))
        return 1;
      if (i > 1) {
        master_printf("ratio R(%g)/R(%g) * 2^6 = %g\n", dpml, dpml / 2,
                      (refl_const / prev_refl_const) * 64.0);
        if ((refl_const / prev_refl_const) * 64.0 > 1.1) return 1;
      }
      prev_refl_const = refl_const;
    }
  }
  master_printf("pml scales correctly with length.\n");
  return 0;
}

int pmlcyl_scaling(double eps(const vec &), int m) {
  double res = 10, freq = 1.0, dpml = 0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double refl_const = 0.0, prev_refl_const = 0.0;
  master_printf("Checking thickness convergence of cylindrical PML for m=%d...\n", m);
  for (int i = 0; i < 3; i++) {
    dpml = pow(2.0, (double)i);
    double sr = 5.0 + dpml, sz = dpml + 5.0 + dpml;
    prev_ft = ft;

    grid_volume gv = volcyl(sr, sz, res);
    gv.center_origin();
    structure s(gv, eps, pml(dpml, Rasymp, stretch));
    fields f(&s, m);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ez, src, veccyl(0.5 * (sr - dpml), 0.1));
    ft = do_ft(f, Ez, veccyl(sr - dpml - 0.1, 0), freq);

    if (i > 0) {
      refl_const = pow(abs(ft - prev_ft), 2.0) / pow(abs(prev_ft), 2.0);
      master_printf("reflcyl:, %g, %g\n", dpml, refl_const);
      if (refl_const > (1e-5) * pow(2 / dpml, 6.0) || refl_const < (1e-8) * pow(2 / dpml, 6.0))
        return 1;
      if (i > 1) {
        master_printf("ratio R(%g)/R(%g) * 2^6 = %g\n", dpml, dpml / 2,
                      (refl_const / prev_refl_const) * 64.0);
        if ((refl_const / prev_refl_const) * 64.0 > 1.1) return 1;
      }
      prev_refl_const = refl_const;
    }
  }
  master_printf("pml scales correctly with length.\n");
  return 0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  master_printf("Running PML tests...\n");
  // if (check_pml1d(one, 0)) meep::abort("not a pml in 1d.");
  if (check_pml1d(one, 10.0)) meep::abort("not a pml in 1d + conductivity.");
  if (check_pml2d(one, Hz, 1, true, 0.5))
    meep::abort("not a pml in 2d TE + conduct. + dispersion + offdiag");
  // if (check_pml2d(one, Ez, 0, false, 0)) meep::abort("not a pml in 2d TM.");
  // if (check_pml2d(one,Ez,1,false,0)) meep::abort("not a pml in 2d TM + conduct.");
  // if (check_pml2d(one,Hz,0,false,0)) meep::abort("not a pml in 2d TE.");
  // if (check_pml2d(one, Hz, 1, false, 0)) meep::abort("not a pml in 2d TE + conduct.");
  //  if (check_pml2d(one,Ez,0,true,0)) meep::abort("not a pml in 2d TM + dispersion.");
  // if (check_pml2d(one, Hz, 0, true, 0)) meep::abort("not a pml in 2d TE + dispersion.");
  // if (check_pml2d(one, Ez, 0, false, 0.5)) meep::abort("not a pml in 2d TM + offdiag.");
  // if (check_pml2d(one, Hz, 0, false, 0.5)) meep::abort("not a pml in 2d TE + offdiag.");
  // if (check_pmlcyl(one)) meep::abort("not a pml in cylincrical co-ordinates.");
  if (sizeof(realnum) == sizeof(double)) {
    if (pml1d_scaling(one)) meep::abort("pml doesn't scale properly with length.");
  }
  if (pmlcyl_scaling(one, 0))
    meep::abort("m=0 cylindrical pml doesn't scale properly with length.");
  if (pmlcyl_scaling(one, 1))
    meep::abort("m=1 cylindrical pml doesn't scale properly with length.");
  if (pmlcyl_scaling(one, 2))
    meep::abort("m=2 cylindrical pml doesn't scale properly with length.");

  return 0;
}

#include <vector>
#include <stdio.h>
#include <meep.hpp>
using namespace meep;
#include "config.h"
using std::complex;

double eps(const vec &pt) { return ((pt.r() < 0.5 + 1e-6) ? 9.0 : 1.0); }

#define MINRES 10
#define MAXRES 25
#define RESSTEP 3 // should be odd

int find_exponent(double a_mean, double a_meansqr, double a2_mean, double a2_meansqr,
                  const char *name) {
  // Verdict on convergence
  double a_sigma, a2_sigma;
  a_sigma = sqrt(a_meansqr - a_mean * a_mean);
  a2_sigma = sqrt(a2_meansqr - a2_mean * a2_mean);
  master_printf("%s a's: ", name);
  if (a2_sigma / a2_mean < 0.15) {
    master_printf("converged as %3.1e / (a*a)\n", a_mean);
    return 2;
  }
  else if (a_sigma / a_mean < 0.15) {
    master_printf("converged as %3.1e / a\n", a_mean);
    return 1;
  }
  else {
    master_printf("Not clear if it converges...\n");
    return 0;
  }
}

void test_convergence_without_averaging() {
  double w0 = 0.2858964; // exact to last digit

  int n[2] = {0, 0};
  double a_mean[2] = {0, 0}, a_meansqr[2] = {0, 0}, a2_mean[2] = {0, 0}, a2_meansqr[2] = {0, 0};

  for (int a = MINRES; a <= MAXRES; a += RESSTEP) {
    grid_volume vol = volcyl(1.0, 0.0, a);
    structure s(vol, eps);
    fields f(&s, 1);
    f.use_bloch(0.1);
    f.set_boundary(High, R, Metallic);
    f.add_point_source(Hr, w0, 2.0, 0.0, 5.0, veccyl(0.2, 0.0));
    while (f.time() < f.last_source_time())
      f.step();
    int t_harminv_max = 2500; // try increasing this in case of failure
    std::vector<complex<double> > mon_data(t_harminv_max);
    int t = 0;
    monitor_point mp;
    while (t < t_harminv_max) {
      f.step();
      f.get_point(&mp, veccyl(0.2, 0.0));
      mon_data[t] = mp.get_component(Er);
      t++;
    }
    const int maxbands = 10;
    std::vector<complex<double> > amps(maxbands);
    std::vector<double> freq_re(maxbands), freq_im(maxbands), errors(maxbands);
    int nfreq = do_harminv(mon_data.data(), t_harminv_max - 1, f.dt, 0.10, 0.50, maxbands,
                           amps.data(), freq_re.data(), freq_im.data(), errors.data());
    double w = 0.0;
    for (int jf = 0; jf < nfreq; jf++)
      if (abs(freq_re[jf] - w0) < abs(w - w0)) w = freq_re[jf];
    double e = -(w - w0) / w0, ea = e * a, ea2 = e * a * a; //  to check 1/a and 1/(a*a) convergence
    // master_printf("Using a = %d ...\n", a);
    // master_printf("a = %3d\tw = %g \t(w-w0)/w0*a = %4.2e \t(w-w0)/w0*a*a = %4.2e\n", a, w, ea,
    // ea2);
    master_printf("noavg:, %d, %g, %g\n", a, w, fabs(e));

    // Statistical analysis
    int index = (2 * (a / 2) == a) ? 0 : 1; // even / odd
    a_mean[index] += ea;
    a_meansqr[index] += ea * ea;
    a2_mean[index] += ea2;
    a2_meansqr[index] += ea2 * ea2;
    n[index]++;
  }
  for (int i = 0; i < 2; i++)
    a_mean[i] /= n[i];
  for (int i = 0; i < 2; i++)
    a_meansqr[i] /= n[i];
  for (int i = 0; i < 2; i++)
    a2_mean[i] /= n[i];
  for (int i = 0; i < 2; i++)
    a2_meansqr[i] /= n[i];

  /* Note: in older versions of Meep, even with "no averaging" there
     was some funny averaging that happened to give quadratic convergence
     for the even-resolution cylindrical case here.  We no longer do this
     -- "no averaging" really means no averaging now. */
  if (find_exponent(a_mean[0], a_meansqr[0], a2_mean[0], a2_meansqr[0], "Even") != 1)
    meep::abort("Failed even convergence test with no fancy averaging!\n");
  if (find_exponent(a_mean[1], a_meansqr[1], a2_mean[1], a2_meansqr[1], "Odd") != 1)
    meep::abort("Failed odd convergence test with no fancy averaging!\n");
  master_printf("Passed convergence test with no fancy averaging!\n");
}

void test_convergence_with_averaging() {
  double w0 = 0.2858964; // exact to last digit

  int n[2] = {0, 0};
  double a_mean[2] = {0, 0}, a_meansqr[2] = {0, 0}, a2_mean[2] = {0, 0}, a2_meansqr[2] = {0, 0};

  for (int a = MINRES; a <= MAXRES; a += RESSTEP) {
    grid_volume vol = volcyl(1.0, 0.0, a);
    structure s(vol, eps);
    s.set_epsilon(eps);

    fields f(&s, 1);
    f.use_bloch(0.1);
    f.set_boundary(High, R, Metallic);
    f.add_point_source(Hr, w0, 2.0, 0.0, 5.0, veccyl(0.2, 0.0));
    while (f.time() < f.last_source_time())
      f.step();
    int t_harminv_max = 2500; // try increasing this in case of failure
    std::vector<complex<double> > mon_data(t_harminv_max);
    int t = 0;
    monitor_point mp;
    while (t < t_harminv_max) {
      f.step();
      f.get_point(&mp, veccyl(0.2, 0.0));
      mon_data[t] = mp.get_component(Er);
      t++;
    }
    const int maxbands = 10;
    std::vector<complex<double> > amps(maxbands);
    std::vector<double> freq_re(maxbands), freq_im(maxbands), errors(maxbands);
    int nfreq = do_harminv(mon_data.data(), t_harminv_max - 1, f.dt, 0.10, 0.50, maxbands,
                           amps.data(), freq_re.data(), freq_im.data(), errors.data());
    double w = 0.0;
    for (int jf = 0; jf < nfreq; jf++)
      if (abs(freq_re[jf] - w0) < abs(w - w0)) w = freq_re[jf];
    double e = -(w - w0) / w0, ea = e * a, ea2 = e * a * a; //  to check 1/a and 1/(a*a) convergence
    // master_printf("Using a = %d ...\n", a);
    // master_printf("a = %3d\tw = %g \t(w-w0)/w0*a = %4.2e \t(w-w0)/w0*a*a = %4.2e\n", a, w, ea,
    // ea2);
    master_printf("avg:, %d, %g, %g\n", a, w, fabs(e));

    // Statistical analysis
    int index = (2 * (a / 2) == a) ? 0 : 1; // even / odd
    a_mean[index] += ea;
    a_meansqr[index] += ea * ea;
    a2_mean[index] += ea2;
    a2_meansqr[index] += ea2 * ea2;
    n[index]++;
  }
  for (int i = 0; i < 2; i++)
    a_mean[i] /= n[i];
  for (int i = 0; i < 2; i++)
    a_meansqr[i] /= n[i];
  for (int i = 0; i < 2; i++)
    a2_mean[i] /= n[i];
  for (int i = 0; i < 2; i++)
    a2_meansqr[i] /= n[i];

  if (find_exponent(a_mean[0], a_meansqr[0], a2_mean[0], a2_meansqr[0], "Even") != 2)
    meep::abort("Failed convergence test with anisotropic dielectric averaging!\n");
  if (find_exponent(a_mean[1], a_meansqr[1], a2_mean[1], a2_meansqr[1], "Odd") != 2)
    meep::abort("Failed convergence test with anisotropic dielectric averaging!\n");
  master_printf("Passed convergence test with anisotropic dielectric averaging!\n");
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
#ifdef HAVE_HARMINV
  master_printf("Testing convergence of a waveguide mode frequency...\n");
  test_convergence_without_averaging();
  test_convergence_with_averaging();
#endif
  return 0;
}

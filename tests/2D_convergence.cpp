#include <meep.hpp>
using namespace meep;
#include "config.h"
using std::complex;

const double diameter = 0.8;
const double r = diameter * 0.5;

double holey_2d(const vec &xx) {
  const grid_volume gv = vol2d(2.0, 1.0, 100.0);
  vec p = xx - gv.center();
  while (p.x() < -0.5)
    p += vec(1.0, 0);
  while (p.x() > 0.5)
    p -= vec(1.0, 0);
  while (p.y() < -0.5)
    p += vec(0, 1.0);
  while (p.y() > 0.5)
    p -= vec(0, 1.0);
  if (fabs(p & p) < r * r - 1e-12) return 1.0;
  return 12.0;
}

double holey_shifted_2d(const vec &xx) { return holey_2d(xx + vec(pi * 0.01, 3 - pi) * 0.5); }

double get_the_freq(monitor_point *p, component c) {
  complex<double> *amp, *freqs;
  int num;
  p->harminv(c, &amp, &freqs, &num, 0.15, 0.20, 8);
  if (!num) return 0.0;
  double best_amp = abs(amp[0]), best_freq = fabs(real(freqs[0]));
  for (int i = 1; i < num; i++)
    if (abs(amp[i]) > best_amp && fabs(imag(freqs[i]) / real(freqs[i])) < 0.002) {
      best_amp = abs(amp[i]);
      best_freq = fabs(real(freqs[i]));
    }
  delete[] freqs;
  delete[] amp;
  return best_freq;
}

double freq_at_resolution(double e(const vec &), double a, component c, double beta) {
  const grid_volume gv = vol2d(2.0, 1.0, a);
  structure s(gv, e);
  s.set_epsilon(e);

  fields f(&s, 0, beta);
  f.use_real_fields();
  f.use_bloch(vec(0, 0));
  f.add_point_source(c, 0.18, 2.5, 0.0, 6.0, vec(0.5, 0.5), 1.0);
  f.add_point_source(c, 0.18, 2.5, 0.0, 6.0, vec(1.5, 0.5), -1.0);

  while (f.time() <= f.last_source_time() + 10.0)
    f.step();
  const double fourier_timesteps = 3000.0;
  const double ttot = fourier_timesteps / a + f.time();
  monitor_point *p = NULL;
  while (f.time() <= ttot) {
    f.step();
    p = f.get_new_point(vec(0.52, 0.97), p);
  }
  const double freq = get_the_freq(p, c);
  delete p;
  return freq;
}

void check_convergence(component c, double best_guess, double beta) {
  const double amin = 5.0, amax = 30.0, adelta = 5.0;

  master_printf("Checking convergence for %s field...\n", component_name(c));
  if (beta != 0) master_printf("... using exp(i beta z) z-dependence with beta=%g\n", beta);
  if (best_guess) master_printf("(The correct frequency should be %g.)\n", best_guess);

  for (double a = amax; a >= amin; a -= adelta) {
    const double freq = freq_at_resolution(holey_2d, a, c, beta);
    const double freq_shifted = freq_at_resolution(holey_shifted_2d, a, c, beta);

    // Initialize best guess at the correct freq.
    if (!best_guess) {
      best_guess = freq + 0.5 * (freq_shifted - freq);
      master_printf("The frequency is approximately %g\n", best_guess);
    }
    else {
      master_printf("frequency for a=%g is %g, %g (shifted), %g (mean)\n", a, freq, freq_shifted,
                    0.5 * (freq + freq_shifted));
      master_printf("Unshifted freq error is %g/%g/%g\n", (freq - best_guess) * a * a, a, a);
      if (fabs(freq - best_guess) * a * a > 0.4)
        meep::abort("Frequency doesn't converge properly with a.\n");
      master_printf("Shifted freq error is %g/%g/%g\n", (freq_shifted - best_guess) * a * a, a, a);
      if (fabs(freq_shifted - best_guess) * a * a > 0.4)
        meep::abort("Frequency doesn't converge properly with a.\n");
    }

    // Check frequency difference...
    master_printf("Frequency difference with a of %g is %g/%g/%g\n", a,
                  (freq - freq_shifted) * a * a, a, a);
    if (fabs(freq - freq_shifted) * a * a > 0.4)
      meep::abort("Frequency difference = doesn't converge properly with a.\n");
  }
  master_printf("Passed 2D resolution convergence test for %s!\n", component_name(c));
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
#ifdef HAVE_HARMINV
  master_printf("Running holes square-lattice resolution convergence test.\n");
  check_convergence(Ey, 0.179944, 0); // from MPB; correct to >= 4 dec. places
  // check_convergence(Ez, 0.166998, 0);  // from MPB; correct to >= 4 dec. places
  check_convergence(Ez, 0.173605, .1); // from MPB; correct to >= 4 dec. places
#endif
  return 0;
}

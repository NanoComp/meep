#include "meep.h"

const double diameter = 0.8;
const double r = diameter*0.5;

double holey_2d(const vec &xx) {
  const volume v = vol2d(2.0,1.0,100.0);
  vec p = xx - v.center();
  while (p.x() <-0.5) p += vec2d(1.0,0);
  while (p.x() > 0.5) p -= vec2d(1.0,0);
  while (p.y() <-0.5) p += vec2d(0,1.0);
  while (p.y() > 0.5) p -= vec2d(0,1.0);
  if (fabs(p & p) < r*r) return 1.0;
  return 12.0;
}

double holey_shifted_2d(const vec &xx) {
  return holey_2d(xx + vec2d(pi*0.01, 3 - pi)*0.5);
}

double get_the_freq(monitor_point *p) {
  complex<double> *amp, *freqs;
  int num;
  p->harminv(Ey, &amp, &freqs, &num, 0.15, 0.20, 8);
  double best_amp = 0.0, best_freq = 0.0;
  for (int i=0;i<num;i++)
    if (abs(amp[i]) > best_amp &&
        fabs(imag(freqs[i])/real(freqs[i])) < 0.002) {
      best_amp = abs(amp[i]);
      best_freq = fabs(real(freqs[i]));
    }
  delete[] freqs;
  delete[] amp;
  return best_freq;
}

double freq_at_resolution(double e(const vec &), double a) {
  const volume v = vol2d(2.0,1.0,a);
  mat ma(v, e);
  ma.set_epsilon(e, 0.0, true);

  fields f(&ma);
  f.use_bloch(vec2d(0,0));
  f.add_point_source(Ey, 0.18, 2.5, 0.0, 6.0, vec2d(0.5,0.5), 1.0);
  f.add_point_source(Ey, 0.18, 2.5, 0.0, 6.0, vec2d(1.5,0.5),-1.0);
  f.use_real_fields();

  while (f.time() <= f.find_last_source() + 10.0 && !interrupt) f.step();
  const double fourier_timesteps = 2000.0;
  const double ttot = fourier_timesteps/a + f.time();
  monitor_point *p = NULL, *p_shifted = NULL;
  while (f.time() <= ttot) {
    f.step();
    p = f.get_new_point(vec2d(0.5,0.5), p);
  }
  const double freq = get_the_freq(p);
  delete p;
  return freq;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  master_printf("Running square holes resolution convergence test.\n");
  const double amin = 5.0, amax = 30.0, adelta = 5.0;
  double best_guess = 0.0;
  for (double a=amax; a >= amin; a-=adelta) {
    const double freq = freq_at_resolution(holey_2d, a);
    const double freq_shifted = freq_at_resolution(holey_shifted_2d, a);

    // Initialize best guess at the correct freq.
    if (!best_guess) {
      best_guess = freq + 0.5*(freq_shifted - freq);
      master_printf("The frequency is approximately %lg\n", best_guess);
    } else {
      master_printf("Unshifted freq error is %lg/%lg/%lg\n",
                    (freq - best_guess)*a*a, a, a);
      if (fabs(freq - best_guess)*a*a > 0.3)
        abort("Frequency doesn't converge properly with a.\n");
      master_printf("Shifted freq error is %lg/%lg/%lg\n",
                    (freq_shifted - best_guess)*a*a, a, a);
      if (fabs(freq_shifted - best_guess)*a*a > 0.3)
        abort("Frequency doesn't converge properly with a.\n");
    }    

    // Check frequency difference...
    master_printf("Frequency difference with a of %lg is %lg/%lg/%lg\n",
                  a, (freq - freq_shifted)*a*a, a, a);
    if (fabs(freq - freq_shifted)*a*a > 0.3)
      abort("Frequency difference = doesn't converge properly with a.\n");
  }
  master_printf("Passed 2D resolution convergence test!\n");
  return 0;
}

#include <meep.hpp>
using namespace meep;

const double sx = 5.0;
const double sy = 3.0;
const double dpml = 1.0;
const double d = 0.35;
const double sw = 1.0;
const double res = 20;

double two_waveguides(const vec &p) {
  if ((fabs(p.x()) >= 0.5 * d) && (fabs(p.x()) <= 0.5 * d + sw) && (p.y() <= 0.5 * sw) &&
      (p.y() >= -0.5 * sw))
    return 11.9;
  else
    return 1.0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;
  grid_volume gv = vol3d(sx + 2 * dpml, sy + 2 * dpml, 0, res);
  gv.center_origin();
  const symmetry S = mirror(X, gv) - mirror(Y, gv);
  structure s(gv, two_waveguides, pml(dpml, X) + pml(dpml, Y), S);
  s.set_epsilon(two_waveguides);
  fields f(&s);
  f.use_bloch(vec(0, 0, 0.5));
  f.add_point_source(Ey, 0.22, 0.06, 0.0, 4.0, vec(0.5 * (d + sw), 0, 0));
  f.add_point_source(Ey, 0.22, 0.06, 0.0, 4.0, vec(-0.5 * (d + sw), 0, 0));

#if 0
  double T = f.last_source_time();
  int iT = T / f.dt;
  while (f.t < iT) {
    if (f.t % (iT / 10) == 0)
      master_printf("%g%% done with source\n", f.time()/T * 100);
    f.step();
  }
  double T2 = 300;
  int iT2 = T2 / f.dt;
  complex<double> *vals = new complex<double>[iT2];
  while (f.t - iT < iT2) {
    if ((f.t - iT) % (iT2 / 10) == 0)
      master_printf("%g%% done with harminv\n", (f.t - iT) * 100.0 / iT2);
    f.step();
    vals[f.t - iT - 1] = f.get_field(Ey, vec(0.5*(d+sw),0.,0.));
  }
  complex<double> amps[8];
  double freqs_re[8], freqs_im[8];

  master_printf("done with timestepping, running harminv...\n");
  int num = do_harminv(vals, iT2, f.dt, 0.19, 0.25, 8,
  		       amps, freqs_re, freqs_im);
  master_printf("harminv found %d modes\n",num);
  for (int i=0;i<num;i++)
    master_printf("freq %d is %0.6g, %0.6g\n",i,freqs_re[i],freqs_im[i]);
  double freq = freqs_re[0];
  double T3 = T + T2;
#else
  double freq = 0.220847; // hard-code since we are testing force, not harminv
  double T3 = f.last_source_time() + 300;
#endif

  double dpad = 0.1;
  volume box(vec(0.5 * d - dpad, -0.5 * sw - dpad, 0),
             vec(0.5 * d + sw + dpad, 0.5 * sw + dpad, 0));
  dft_flux fluxR = f.add_dft_flux_plane(box, freq, freq, 1);
  dpad = 0.02;
  volume_list line(volume(vec(0.5 * d - dpad, -0.5 * sy, 0), vec(0.5 * d - dpad, 0.5 * sy, 0)), Sx);
  dft_force forceR = f.add_dft_force(&line, freq, freq, 1);
  f.zero_fields();
  f.t = 0;
  while (f.time() < T3) {
    if (f.t % 2000 == 0) master_printf("%g%% done with flux and force\n", f.time() / T3 * 100);
    f.step();
  }
  double *fl = fluxR.flux();
  double *fr = forceR.force();
  master_printf("flux is %0.8g, force is %0.8g, F/P = %0.8g\n", fl[0], fr[0], -0.5 * fr[0] / fl[0]);
  double FoverP = -0.5 * fr[0] / fl[0];
  delete[] fl;
  delete[] fr;
  return fabs(FoverP + 0.33628872) / 0.33628872 > 0.1;
  // MPB: -0.33628872
}

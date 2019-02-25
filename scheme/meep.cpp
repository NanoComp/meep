#include "meep-ctl.hpp"

using namespace meep;
using namespace std;

/**************************************************************************/

/* The following are hook functions called from main() when
   starting the program and just before exiting.  */

static initialize *meep_init = 0;

void ctl_start_hook(int *argc, char ***argv) {
  meep_init = new initialize(*argc, *argv);
#ifdef HAVE_LIBCTL_QUIET
  extern int libctl_quiet;
  libctl_quiet = !am_master();
#endif
}

void ctl_stop_hook(void) { delete meep_init; }

extern "C" void SWIG_init();

void ctl_export_hook(void) { SWIG_init(); }

/**************************************************************************/

ctlio::cvector3_list do_harminv(ctlio::cnumber_list vals, double dt, double fmin, double fmax,
                                int maxbands, double spectral_density, double Q_thresh,
                                double rel_err_thresh, double err_thresh, double rel_amp_thresh,
                                double amp_thresh) {
  complex<double> *amp = new complex<double>[maxbands];
  double *freq_re = new double[maxbands];
  double *freq_im = new double[maxbands];
  double *freq_err = new double[maxbands];
  maxbands = do_harminv(reinterpret_cast<complex<double> *>(vals.items), vals.num_items, dt, fmin,
                        fmax, maxbands, amp, freq_re, freq_im, freq_err, spectral_density, Q_thresh,
                        rel_err_thresh, err_thresh, rel_amp_thresh, amp_thresh);
  ctlio::cvector3_list res;
  res.num_items = maxbands;
  res.items = new cvector3[maxbands];
  for (int i = 0; i < maxbands; ++i) {
    res.items[i].x.re = freq_re[i];
    res.items[i].x.im = freq_im[i];
    res.items[i].y.re = real(amp[i]);
    res.items[i].y.im = imag(amp[i]);
    res.items[i].z.re = freq_err[i];
    res.items[i].z.im = 0;
  }
  delete[] freq_err;
  delete[] freq_im;
  delete[] freq_re;
  delete[] amp;
  return res;
}

kpoint_list do_get_eigenmode_coefficients(fields *f, dft_flux flux, const volume &eig_vol,
                                          int *bands, int num_bands, int parity,
                                          std::complex<double> *coeffs, double *vgrp,
                                          double eig_resolution, double eigensolver_tol,
                                          meep::kpoint_func user_kpoint_func,
                                          void *user_kpoint_data) {
  size_t num_kpoints = num_bands * flux.Nfreq;
  meep::vec *kpoints = new meep::vec[num_kpoints];
  meep::vec *kdom = new meep::vec[num_kpoints];

  f->get_eigenmode_coefficients(flux, eig_vol, bands, num_bands, parity, eig_resolution,
                                eigensolver_tol, coeffs, vgrp, user_kpoint_func, user_kpoint_data,
                                kpoints, kdom);

  kpoint_list res = {kpoints, num_kpoints, kdom, num_kpoints};

  return res;
}

/**************************************************************************/

/* This is a wrapper function to fool SWIG...since our list constructor
   takes ownership of the next pointer, we have to make sure that SWIG
   does not garbage-collect volume_list objects.  We do
   this by wrapping a "helper" function around the constructor which
   does not have the %newobject SWIG attribute.   Note that we then
   need to deallocate the list explicitly in Scheme. */
volume_list *make_volume_list(const volume &v, int c, complex<double> weight, volume_list *next) {
  return new volume_list(v, c, weight, next);
}

/***************************************************************************/

ctlio::number_list dft_flux_flux(dft_flux *f) {
  ctlio::number_list res;
  res.num_items = f->Nfreq;
  res.items = f->flux();
  return res;
}

ctlio::number_list dft_energy_electric(dft_energy *f)
{
  ctlio::number_list res;
  res.num_items = f->Nfreq;
  res.items = f->electric();
  return res;
}

ctlio::number_list dft_energy_magnetic(dft_energy *f)
{
  ctlio::number_list res;
  res.num_items = f->Nfreq;
  res.items = f->magnetic();
  return res;
}

ctlio::number_list dft_energy_total(dft_energy *f)
{
  ctlio::number_list res;
  res.num_items = f->Nfreq;
  res.items = f->total();
  return res;
}

ctlio::number_list dft_force_force(dft_force *f) {
  ctlio::number_list res;
  res.num_items = f->Nfreq;
  res.items = f->force();
  return res;
}

ctlio::number_list dft_ldos_ldos(dft_ldos *f) {
  ctlio::number_list res;
  res.num_items = f->Nomega;
  res.items = f->ldos();
  return res;
}

ctlio::cnumber_list dft_ldos_F(dft_ldos *f) {
  ctlio::cnumber_list res;
  res.num_items = f->Nomega;
  res.items = (cnumber *)f->F();
  return res;
}

ctlio::cnumber_list dft_ldos_J(dft_ldos *f) {
  ctlio::cnumber_list res;
  res.num_items = f->Nomega;
  res.items = (cnumber *)f->J();
  return res;
}

ctlio::cnumber_list dft_near2far_farfield(dft_near2far *f, const vec &x) {
  ctlio::cnumber_list res;
  res.num_items = f->Nfreq * 6;
  res.items = (cnumber *)f->farfield(x);
  return res;
}

ctlio::number_list dft_near2far_flux(dft_near2far *f, direction df, const volume &where,
                                     double resolution) {
  ctlio::number_list res;
  res.num_items = f->Nfreq;
  res.items = f->flux(df, where, resolution);
  return res;
}

/***************************************************************************/

ctlio::cnumber_list make_casimir_g(double T, double dt, double sigma, meep::field_type ft,
                                   complex<double> (*eps_func)(complex<double> omega),
                                   double Tfft) {
  ctlio::cnumber_list res;
  res.num_items = int(ceil(T / dt));
  res.items = new cnumber[res.num_items];
  complex<double> *g = meep::make_casimir_gfunc(T, dt, sigma, ft, eps_func, Tfft);
  for (int i = 0; i < res.num_items; ++i) {
    res.items[i].re = real(g[i]);
    res.items[i].im = imag(g[i]);
  }
  delete[] g;
  return res;
}

ctlio::cnumber_list make_casimir_g_kz(double T, double dt, double sigma, meep::field_type ft) {
  ctlio::cnumber_list res;
  res.num_items = int(ceil(T / dt));
  res.items = new cnumber[res.num_items];
  complex<double> *g = meep::make_casimir_gfunc_kz(T, dt, sigma, ft);
  for (int i = 0; i < res.num_items; ++i) {
    res.items[i].re = real(g[i]);
    res.items[i].im = imag(g[i]);
  }
  delete[] g;
  return res;
}

#include <stdio.h>
#include <meep.h>
using namespace meep;

double guided_eps(const vec &v) { return ((v.r() < 0.5+1e-6) ? 9.0 : 1.0); }

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  master_printf("Testing convergence of a waveguide mode frequency...\n");
  double w0 = 0.2858964; // exact to last digit 

  int n[2] = {0,0};
  double a_mean[2] = {0,0}, a_meansqr[2] = {0,0}, a2_mean[2] = {0,0}, a2_meansqr[2] = {0,0}; 

  for (int a=10; a <= 25; a+=3) {
  volume vol = volcyl(1.0,0.0,a);  
  structure s(vol, guided_eps, 0);
  
  fields f(&s, 1);
  f.use_bloch(0.1);
  f.set_boundary(High, R, Metallic);
  f.add_point_source(Hr, w0, 2.0, 0.0, 5.0, vec(0.2,0.0));
  while (f.time() < f.find_last_source()) f.step();
  int t_harminv_max = 2500; // try increasing this in case of failure
  complex<double> *mon_data = new complex<double>[t_harminv_max];
  int t = 0;
  monitor_point mp;
  while (t < t_harminv_max) {
      f.step();
      f.get_point(&mp,  vec(0.2,0.0));
      mon_data[t] = mp.get_component(Er);
      t++;
  }
  int maxbands = 10, nfreq;
  complex<double> *amps = new complex<double>[maxbands]; ;
  double *freq_re = new double[maxbands], *freq_im = new double[maxbands], *errors  = new double[maxbands];
  nfreq = do_harminv(mon_data, t_harminv_max - 1, 1, a, 0.10, 0.50, maxbands, amps, freq_re, freq_im, errors);
  double w = 0.0;
  for (int jf = 0; jf < nfreq; jf++) 
    if (abs(freq_re[jf] - w0) < 0.03)
      w = freq_re[jf];
  double e = -(w-w0)/w0, ea = e*a, ea2=e*a*a; //  to check 1/a and 1/(a*a) convergence
  master_printf("Using a = %d ...\n", a);
  //master_printf("a = %3d\tw = %lg \t(w-w0)/w0*a = %4.2e \t(w-w0)/w0*a*a = %4.2e\n", a, w, ea, ea2);

    // Statistical analysis
    int index = (2*(a/2)==a) ? 0 : 1; // even / odd
    a_mean[index]     = (n[index]*a_mean[index] + ea)          / (n[index] + 1);
    a_meansqr[index]  = (n[index]*a_meansqr[index] + ea*ea)    / (n[index] + 1);
    a2_mean[index]    = (n[index]*a2_mean[index] + ea2)        / (n[index] + 1);
    a2_meansqr[index] = (n[index]*a2_meansqr[index] + ea2*ea2) / (n[index] + 1);
    n[index]++;
  }
  
  // Verdict on convergence
  double a_sigma[2], a2_sigma[2];
  int convergence_exponent[2];
  for (int index = 0; index <=1; index++) { 
    a_sigma[index] = sqrt(a_meansqr[index] - a_mean[index]*a_mean[index]);
    a2_sigma[index] = sqrt(a2_meansqr[index] - a2_mean[index]*a2_mean[index]);
    master_printf("%s a's: ", (index==0) ? "Even" : "Odd");
    if (a2_sigma[index]/a2_mean[index] < 0.15) {
      master_printf("converged as %3.1e / (a*a)\n", a_mean[index]);
      convergence_exponent[index] = 2;
    }
    else if (a_sigma[index]/a_mean[index] < 0.15) {
      master_printf("converged as %3.1e / a\n", a_mean[index]);
      convergence_exponent[index] = 1;
    }
    else {
      master_printf("Not clear if it converges...\n"); 
      convergence_exponent[index] = 0;
    }
  }
  if (convergence_exponent[0] != 2  || convergence_exponent[1] != 1)
    abort("Failed convergence test!\n");
  else
    master_printf("Passed convergence test!\n");
}

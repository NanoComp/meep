#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
using namespace meep;

double one(const vec &) { return 1.0; }

int check_pml1d(double eps(const vec &)) {
  double freq = 1.0, dpml = 1.0;
  double sz = 10.0 + 2*dpml;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  master_printf("Checking resolution convergence of 1d PML...\n");
  for (int i=0; i<15; i++) {
    double res = 50.0 + 10.0*i;
    master_printf("    checking 1d resolution %g...\n", res);
    volume v = vol1d(sz,res);
    v.center_origin();
    structure s(v, eps, pml(dpml));
    fields f(&s);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ex, src, vec(-0.5*sz+dpml+0.1));
    vec fpt(0.5*sz - dpml - 0.1);
    if (i > 1) prev_refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
    prev_ft = ft;
    ft = 0.0;
    double emax = 0;
    while (f.time() < f.last_source_time()) {
      ft += f.get_field(Ex, fpt) * polar(1.0, 2*pi*freq * f.time());
      emax = max(emax, abs(f.get_field(Ex, fpt)));
      f.step();
    }
    do {
      double emaxcur = 0;
      double T = f.time() + 50;
      while (f.time() < T) {
	ft += f.get_field(Ex, fpt) * polar(1.0, 2*pi*freq * f.time());
	double e = abs(f.get_field(Ex, fpt));
	emax = max(emax, e);
	emaxcur = max(emaxcur, e);
	f.step();
      }
      if (emaxcur < 1e-6 * emax) break;
    } while(1);
    if (i > 1) {
      refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
      if (refl_const > prev_refl_const) return 1;
    }
  }
  master_printf("passed 1d PML check.\n");
  return 0;
}

int check_pml2dTM(double eps(const vec &)) {
  double freq = 1.0, dpml = 1.0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  double sxy = 5.0 + 2*dpml;
  vec fpt(0.5*sxy - dpml - 0.1,0);
  master_printf("Checking resolution convergence of 2d TM PML...\n");
  for (int i=0; i<5; i++) {
    double res = 10.0 + 5.0*i;
    master_printf("    checking 2d TM resolution %g...\n", res);
    volume v = vol2d(sxy,sxy,res);
    v.center_origin();
    const symmetry S = mirror(X,v) + mirror(Y,v);
    structure s(v, eps, pml(dpml), S);
    fields f(&s);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ez, src, v.center());
    if (i > 1) 
      prev_refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
    prev_ft = ft;
    ft = 0.0;
    double emax = 0;
    while (f.time() < f.last_source_time()) {
      ft += f.get_field(Ez, fpt) * polar(1.0, 2*pi*freq * f.time());
      emax = max(emax, abs(f.get_field(Ez, fpt)));
      f.step();
    }
    do {
      double emaxcur = 0;
      double T = f.time() + 50;
      while (f.time() < T) {
	ft += f.get_field(Ez, fpt) * polar(1.0, 2*pi*freq * f.time());
	double e = abs(f.get_field(Ez, fpt));
	emax = max(emax, e);
	emaxcur = max(emaxcur, e);
	f.step();
      }
      if (emaxcur < 1e-6 * emax) break;
    } while(1);
    if (i > 1) {
      refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
      if (refl_const > prev_refl_const) return 1;      
    }
  }
  master_printf("passed 2d TM PML check.\n");
  return 0;
}

int check_pml2dTE(double eps(const vec &)) {
  double freq = 1.0, dpml = 1.0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  double sxy = 5.0 + 2*dpml;
  vec fpt(0.5*sxy - dpml - 0.1,0);
  master_printf("Checking resolution convergence of 2d TE PML...\n");
  for (int i=0; i<5; i++) {
    double res = 10.0 + 5.0*i;
    master_printf("    checking 2d TE resolution %g...\n", res);
    volume v = vol2d(sxy,sxy,res);
    v.center_origin();
    const symmetry S = -mirror(X,v) - mirror(Y,v);
    structure s(v, eps, pml(dpml), S);
    fields f(&s);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Hz, src, v.center());
    if (i > 1) 
      prev_refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
    prev_ft = ft;
    ft = 0.0;
    double emax = 0;
    while (f.time() < f.last_source_time()) {
      ft += f.get_field(Hz, fpt) * polar(1.0, 2*pi*freq * f.time());
      emax = max(emax, abs(f.get_field(Hz, fpt)));
      f.step();
    }
    do {
      double emaxcur = 0;
      double T = f.time() + 50;
      while (f.time() < T) {
	ft += f.get_field(Hz, fpt) * polar(1.0, 2*pi*freq * f.time());
	double e = abs(f.get_field(Hz, fpt));
	emax = max(emax, e);
	emaxcur = max(emaxcur, e);
	f.step();
      }
      if (emaxcur < 1e-6 * emax) break;
    } while(1);
    if (i > 1) {
      refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
      if (refl_const > prev_refl_const) return 1;      
    }
  }
  master_printf("passed 2d TE PML check.\n");
  return 0;
}

int check_pml3d(double eps(const vec &)) {
  double freq = 1.0, dpml = 1.0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  double sxyz = 3.0 + 2*dpml;
  vec fpt(0.5*sxyz - dpml - 0.1,0);
  master_printf("Checking resolution convergence of 3d PML...\n");
  for (int i=0; i<5; i++) {
    double res = 6.0 + 2.0*i;
    master_printf("    checking 3d resolution %g...\n", res);
    volume v = vol3d(sxyz,sxyz,sxyz,res);
    v.center_origin();
    const symmetry S = mirror(X,v) + mirror(Y,v) - mirror(Z, v);
    structure s(v, eps, pml(dpml), S);
    fields f(&s);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ez, src, vec(0,0,0));
    if (i > 1) prev_refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
    prev_ft = ft;
    ft = 0.0;
    double emax = 0;
    while (f.time() < f.last_source_time()) {
      ft += f.get_field(Hz, fpt) * polar(1.0, 2*pi*freq * f.time());
      emax = max(emax, abs(f.get_field(Hz, fpt)));
      f.step();
    }
    do {
      double emaxcur = 0;
      double T = f.time() + 50;
      while (f.time() < T) {
	ft += f.get_field(Hz, fpt) * polar(1.0, 2*pi*freq * f.time());
	double e = abs(f.get_field(Hz, fpt));
	emax = max(emax, e);
	emaxcur = max(emaxcur, e);
	f.step();
      }
      if (emaxcur < 1e-6 * emax) break;
    } while(1);
    if (i > 1) {
      refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
      if (refl_const > prev_refl_const) return 1;
    }
  }
  return 0;
}

int check_pmlcyl(double eps(const vec &)) {
  double freq = 1.0, dpml = 1.0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double prev_refl_const = 0.0, refl_const = 0.0;
  double sr = 5.0 + dpml, sz = 1.0;
  vec fpt = veccyl(sr - dpml - 0.1,0);
  master_printf("Checking resolution convergence of cylindrical PML...\n");
  for (int i=0; i<5; i++) {
    double res = 10.0 + 5.0*i;
    master_printf("    checking cylindrical resolution %g...\n", res);
    volume v = volcyl(sr,sz,res);
    structure s(v, eps, pml(dpml));
    fields f(&s, 0);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ez, src, veccyl(0.1,0.1));
    if (i > 1) 
      prev_refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
    prev_ft = ft;
    ft = 0.0;
    double emax = 0;
    while (f.time() < f.last_source_time()) {
      ft += f.get_field(Ez, fpt) * polar(1.0, 2*pi*freq * f.time());
      emax = max(emax, abs(f.get_field(Ez, fpt)));
      f.step();
    }
    do {
      double emaxcur = 0;
      double T = f.time() + 50;
      while (f.time() < T) {
	ft += f.get_field(Ez, fpt) * polar(1.0, 2*pi*freq * f.time());
	double e = abs(f.get_field(Ez, fpt));
	emax = max(emax, e);
	emaxcur = max(emaxcur, e);
	f.step();
      }
      if (emaxcur < 1e-6 * emax) break;
    } while(1);
    if (i > 1) {
      refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
      if (refl_const > prev_refl_const) return 1;      
    }
  }
  master_printf("passed cylindrical PML check.\n");
  return 0;
}

int pml1d_scaling(double eps(const vec &)) {
  double res = 20, freq = 1.0, dpml = 0;
  complex<double> prev_ft = 0.0, ft = 0.0;
  double refl_const = 0.0;
  master_printf("Checking thickness convergence of 1d PML...\n");
  for (int i=0; i<10; i++) {
    dpml = pow(2.0,(double)i);
    master_printf("    checking 1d thickness %g...\n", dpml);
    double sz = 2*dpml + 10.0 + dpml;
    volume v = vol1d(sz,res);
    structure s(v, eps, (pml(2*dpml,Z,Low)*2.0 + pml(dpml,Z,High))*2);
    fields f(&s);
    gaussian_src_time src(freq, freq / 20);
    f.add_point_source(Ex, src, vec(2*dpml+0.1));
    vec fpt(sz - dpml - 0.1);
    prev_ft = ft;
    ft = 0.0;
    double emax = 0;
    while (f.time() < f.last_source_time()) {
      ft += f.get_field(Ex, fpt) * polar(1.0, 2*pi*freq * f.time());
      emax = max(emax, abs(f.get_field(Ex, fpt)));
      f.step();
    }
    do {
      double emaxcur = 0;
      double T = f.time() + 50;
      while (f.time() < T) {
	ft += f.get_field(Ex, fpt) * polar(1.0, 2*pi*freq * f.time());
	double e = abs(f.get_field(Ex, fpt));
	emax = max(emax, e);
	emaxcur = max(emaxcur, e);
	f.step();
      }
      if (emaxcur < 1e-6 * emax) break;
    } while(1);
    if (i > 0) {
      refl_const = pow(abs(ft - prev_ft),2.0) / pow(abs(prev_ft),2.0);
      if (refl_const > (1e-6)*pow(1/dpml,6.0)) return 1;
    }      
  }
  master_printf("pml scales correctly with length.\n");
  return 0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  quiet = true;
  master_printf("Running PML tests...\n");
  if (check_pml1d(one)) abort("not a pml in 1d.");
  if (check_pml2dTM(one)) abort("not a pml in 2d TM.");
  if (check_pml2dTE(one)) abort("not a pml in 2d TE."); 
  if (check_pmlcyl(one)) abort("not a pml in cylincrical co-ordinates.");
  // if (check_pml3d(one)) abort("not a pml in 3d.");
  if (pml1d_scaling(one)) abort("pml doesn't scale properly with length.");
  return 0;
}

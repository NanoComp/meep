/* Copyright (C) 2003 Massachusetts Institute of Technology
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
#include <math.h>

#include "dactyl.h"
#include "dactyl_internals.h"

mat::mat(const volume &thev, double eps(const vec &)) {
  num_chunks = 1;
  v = thev;
  chunks = new (mat_chunk *)[num_chunks];
  chunks[0] = new mat_chunk(v,eps);
}

mat::mat(const mat *m) {
  num_chunks = m->num_chunks;
  v = m->v;
  chunks = new (mat_chunk *)[num_chunks];
  for (int i=0;i<num_chunks;i++) chunks[i] = new mat_chunk(m->chunks[i]);
}

mat::~mat() {
  for (int i=0;i<num_chunks;i++) {
    delete chunks[i];
  }
  delete[] chunks;
}

void mat::make_average_eps() {
  for (int i=0;i<num_chunks;i++) chunks[i]->make_average_eps();
}

void mat::use_pml_left(double dx) { // FIXME
  for (int i=0;i<num_chunks;i++) chunks[i]->use_pml_left(dx);
}
void mat::use_pml_right(double dx) { // FIXME
  for (int i=0;i<num_chunks;i++) chunks[i]->use_pml_right(dx);
}
void mat::use_pml_radial(double dx) { // FIXME
  for (int i=0;i<num_chunks;i++) chunks[i]->use_pml_radial(dx);
}
void mat::set_output_directory(const char *name) {
  outdir = name;
  for (int i=0;i<num_chunks;i++) chunks[i]->set_output_directory(name);
}
void mat::mix_with(const mat *oth, double f) {
  if (num_chunks != oth->num_chunks) {
    printf("You can't phase materials with different chunk topologies...\n");
    exit(1);
  }
  for (int i=0;i<num_chunks;i++) chunks[i]->mix_with(oth->chunks[i], f);
}

mat_chunk::~mat_chunk() {
  for (int c=0;c<10;c++) delete[] inveps[c];
  delete[] eps;

  for (int c=0;c<10;c++) delete[] Cmain[c];
  for (int c=0;c<10;c++) delete[] Cother[c];
  if (pb) delete pb;
}

static double sig(double r, double power);

static double minimize_badness(double sig[], int thickness, double eps, double fmin, int i);
inline void reverse(double sig[], int l) {
  for (int i=0;i<l/2;i++) {
    double temp = sig[i];
    sig[i] = sig[l-1-i];
    sig[l-1-i] = temp;
  }
}

static double badness(double sig[], int thickness, double epsilon, double fmin) {
  if (thickness < 1) return 1;
  const double A = .0001/fmin*.1/fmin, K = 6.0/epsilon*2.25/epsilon;
  double sofar = 1.0;
  for (int i=0;i<thickness-1;i++) {
    double first_trans = exp(-K*sig[i+1]);
    double refl = A*fabs(sig[i]-sig[i+1])*fabs(sig[i]-sig[i+1]);
    double total_trans = exp(-K*sig[i])*first_trans;
    sofar = refl + (1-refl)*total_trans*sofar;
    if (sofar > 1.0) sofar = 1.0;
  }
  double last_refl = A*fabs(sig[thickness-1]);
  sofar = last_refl + (1-last_refl)*sofar;
  return sofar;
}

static double minimize_badness(double sig[], int thickness,
                               double epsilon, double fmin, int i) {
  double behind_reflection = badness(sig, i-1, epsilon, fmin);
  

  double now = badness(sig, thickness, epsilon, fmin);
  double tried = now;
  do {
    now = tried;
    sig[i] *= 1.001;
    tried = badness(sig, thickness, epsilon, fmin);
  } while (tried < now);
  sig[i] /= 1.001;
  tried = now = badness(sig, thickness, epsilon, fmin);
  do {
    now = tried;
    sig[i] /= 1.001;
    tried = badness(sig, thickness, epsilon, fmin);
  } while (tried < now);
  sig[i] *= 1.001;
  return badness(sig, thickness, epsilon, fmin);
}

static double sig(double r, double power) {
  return pow(r, power);
}

void mat_chunk::mix_with(const mat_chunk *n, double f) {
  for (int i=0;i<v.ntot();i++)
    eps[i] = 1.0/(1.0/eps[i] + f*(1.0/n->eps[i]-1.0/eps[i]));
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c))
      for (int i=0;i<v.ntot();i++)
        inveps[c][i] += f*(n->inveps[c][i] - inveps[c][i]);
  // Mix in the polarizability...
  polarizability *po = pb, *pn = n->pb;
  while (po && pn) {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_electric((component)c))
        for (int i=0;i<v.ntot();i++)
          po->s[c][i] += f*(pn->s[c][i] - po->s[c][i]);
    po = po->next;
    pn = pn->next;
  }
}

void mat_chunk::make_average_eps() {
  double meaneps = 0;
  for (int i=0;i<v.ntot();i++) {
    meaneps += eps[i]; // This isn't quite correct as some points are counted twice...
  }
  meaneps /= v.ntot();
  for (int i=0;i<v.ntot();i++)
    eps[i] = meaneps;
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c))
      for (int i=0;i<v.ntot();i++)
        inveps[c][i] = 1/meaneps;
}

const double Cmax = 0.5;

void mat_chunk::use_pml_right(double dx) {
  if (v.dim == d1) {
    const double border = v.nz()/v.a + v.origin.z() - dx;
    const double prefac = Cmax/(dx*dx);
    if (!Cmain[Hy]) {
      Cmain[Hy] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Hy][i] = 0.0;
    }
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(Hy,i).z();
      if (x > border) Cmain[Hy][i] = prefac*(x-border)*(x-border);
    }
    if (!Cmain[Ex]) {
      Cmain[Ex] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Ex][i] = 0.0;
    }
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(Ex,i).z();
      if (x > border) Cmain[Ex][i] = prefac*(x-border)*(x-border);
    }
  } else if (v.dim == dcyl) {
    const double border = v.nz()/v.a + v.origin.z() - dx;
    const double prefac = Cmax/(dx*dx);
    if (!Cmain[Ep]) {
      Cmain[Ep] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Ep][i] = 0.0;
    }
    if (!Cmain[Hp]) {
      Cmain[Hp] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Hp][i] = 0.0;
    }
    if (!Cother[Er]) {
      Cother[Er] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Er][i] = 0.0;
    }
    if (!Cother[Hr]) {
      Cother[Hr] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Hr][i] = 0.0;
    }
    component m=Er;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x > border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Hr;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x > border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Ep;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x > border) Cmain[m][i] = prefac*(x-border)*(x-border);
    }
    m=Hp;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x > border) Cmain[m][i] = prefac*(x-border)*(x-border);
    }
  } else {
    printf("Unsupported dimension?!\n");
    exit(1);
  }
}

void mat_chunk::use_pml_left(double dx) {
  if (v.dim == d1) {
    const double border = dx + v.origin.z();
    const double prefac = Cmax/(dx*dx);
    if (!Cmain[Hy]) {
      Cmain[Hy] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Hy][i] = 0.0;
    }
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(Hy,i).z();
      if (x < border) Cmain[Hy][i] = prefac*(x-border)*(x-border);
    }
    if (!Cmain[Ex]) {
      Cmain[Ex] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Ex][i] = 0.0;
    }
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(Ex,i).z();
      if (x < border) Cmain[Ex][i] = prefac*(x-border)*(x-border);
    }
  } else if (v.dim == dcyl) {
    const double border = dx + v.origin.z();
    const double prefac = Cmax/(dx*dx);
    if (!Cmain[Ep]) {
      Cmain[Ep] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Ep][i] = 0.0;
    }
    if (!Cmain[Hp]) {
      Cmain[Hp] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cmain[Hp][i] = 0.0;
    }
    if (!Cother[Er]) {
      Cother[Er] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Er][i] = 0.0;
    }
    if (!Cother[Hr]) {
      Cother[Hr] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Hr][i] = 0.0;
    }
    component m=Er;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x < border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Hr;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x < border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Ep;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x < border) Cmain[m][i] = prefac*(x-border)*(x-border);
    }
    m=Hp;
    for (int i=0;i<v.ntot();i++) {
      double x = v.loc(m,i).z();
      if (x < border) Cmain[m][i] = prefac*(x-border)*(x-border);
    }
  } else {
    printf("Unsupported dimension?!\n");
    exit(1);
  }
}

void mat_chunk::use_pml_radial(double dx) {
  if (v.dim == dcyl) {
    const double border = dx + v.origin.z();
    const double prefac = Cmax/(dx*dx);
    if (!Cother[Ep]) {
      Cother[Ep] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Ep][i] = 0.0;
    }
    if (!Cother[Hp]) {
      Cother[Hp] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Hp][i] = 0.0;
    }
    if (!Cother[Ez]) {
      Cother[Ez] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Ez][i] = 0.0;
    }
    if (!Cother[Hz]) {
      Cother[Hz] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) Cother[Hz][i] = 0.0;
    }
    component m=Ez;
    for (int i=0;i<v.ntot();i++) {
      double x = v.origin.r() + v.nr()/v.a - v.loc(m,i).r();
      if (x < border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Hz;
    for (int i=0;i<v.ntot();i++) {
      double x = v.origin.r() + v.nr()/v.a - v.loc(m,i).r();
      if (x < border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Ep;
    for (int i=0;i<v.ntot();i++) {
      double x = v.origin.r() + v.nr()/v.a - v.loc(m,i).r();
      if (x < border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
    m=Hp;
    for (int i=0;i<v.ntot();i++) {
      double x = v.origin.r() + v.nr()/v.a - v.loc(m,i).r();
      if (x < border) Cother[m][i] = prefac*(x-border)*(x-border);
    }
  } else {
    printf("Radial pml only works in cylindrical coordinates!\n");
    exit(1);
  }
}

mat_chunk::mat_chunk(const mat_chunk *o) {
  outdir = o->outdir;
  if (o->pb) pb = new polarizability(o->pb);
  else pb = NULL;
  a = o->a;
  v = o->v;
  eps = new double[v.ntot()];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }
  for (int i=0;i<v.ntot();i++) eps[i] = o->eps[i];
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c)) {
      inveps[c] = new double[v.ntot()];
      for (int i=0;i<v.ntot();i++) inveps[c][i] = o->inveps[c][i];
    } else {
      inveps[c] = NULL;
    }
  // Allocate the conductivity arrays:
  for (int c=0;c<10;c++) {
    Cmain[c] = NULL;
    Cother[c] = NULL;
  }
  // Copy over the conductivity arrays:
  for (int c=0;c<10;c++)
    if (v.has_field((component)c)) {
      if (o->Cmain[c]) {
        Cmain[c] = new double[v.ntot()];
        for (int i=0;i<v.ntot();i++) Cmain[c][i] = o->Cmain[c][i];
      }
      if (o->Cother[c]) {
        Cother[c] = new double[v.ntot()];
        for (int i=0;i<v.ntot();i++) Cother[c][i] = o->Cother[c][i];
      }
    }
}

mat_chunk::mat_chunk(const volume &thev, double feps(const vec &)) {
  pml_fmin = 0.2;
  outdir = ".";
  pb = NULL;
  v = thev;
  a = thev.a;
  eps = new double[v.ntot()];
  if (eps == NULL) {
    printf("Out of memory!\n");
    exit(1);
  }

  if (v.dim == dcyl) {
    for (int i=0;i<v.ntot();i++) eps[i] = feps(v.loc(Hp,i));
  } else if (v.dim == d1) {
    for (int i=0;i<v.ntot();i++) eps[i] = feps(v.loc(Ex,i));
  } else {
    printf("Unsupported symmetry!\n");
    exit(1);
  } 
  for (int c=0;c<10;c++)
    if (v.has_field((component)c) && is_electric((component)c)) {
      inveps[c] = new double[v.ntot()];
      // Initialize eps to 1;
      for (int i=0;i<v.ntot();i++) inveps[c][i] = 1;
    } else {
      inveps[c] = NULL;
    }

  if (v.dim == dcyl) {
    const vec dr = v.dr()*0.5; // The distance between Yee field components
    const vec dz = v.dz()*0.5; // The distance between Yee field components
    for (int r=1;r<v.nr();r++) {
      const int ir = r*(v.nz()+1);
      const int irm1 = (r-1)*(v.nz()+1);
      for (int z=1;z<=v.nz();z++) {
        inveps[Er][z + ir] = 2./(eps[z+ir] + eps[z+ir-1]);
        inveps[Ep][z + ir] = 4./(eps[z+ir] + eps[z+ir-1] +
                               eps[z+irm1] + eps[z+irm1-1]);
        inveps[Ez][z + ir] = 2./(eps[z+ir] + eps[z+irm1]);
      }
    }
    for (int r=0;r<v.nr();r++) {
      const int ir = r*(v.nz()+1);
      const vec here = v.loc(Ep,ir);
      inveps[Er][ir] = 2./(feps(here+dr+dz) + feps(here+dr-dz));
      inveps[Ep][ir] = 4./(feps(here+dr+dz) + feps(here-dr+dz) +
                         feps(here+dr-dz) + feps(here-dr-dz));
      inveps[Ez][ir] = 2./(feps(here+dr+dz) + feps(here-dr+dz));
    }
    for (int z=0;z<v.nz();z++) {
      const vec here = v.loc(Ep,z);
      inveps[Er][z] = 2./(feps(here+dr+dz) + feps(here+dr-dz));
      inveps[Ep][z] = 4./(feps(here+dr+dz) + feps(here-dr+dz) +
                        feps(here+dr-dz) + feps(here-dr-dz));
      inveps[Ez][z] = 2./(feps(here+dr+dz) + feps(here-dr+dz));
    }
  } else if (v.dim == d1) {
    for (int i=0;i<v.ntot();i++) inveps[Ex][i] = 1.0/eps[i];
  } else {
    printf("Unsupported symmetry!\n");
    exit(1);
  }
  // Allocate the conductivity arrays:
  for (int c=0;c<10;c++) {
    Cmain[c] = NULL;
    Cother[c] = NULL;
  }
}

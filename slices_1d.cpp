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

#include "tidod.h"
#include "tidod_internals.h"

static double get_phase(double *f[2], int nz) {
  complex<double> mean=0;
  for (int z=0;z<nz;z++) {
    mean += complex<double>(RE(f,z), IM(f,z));
  }
  complex<double> meanr=0;
  for (int z=0;z<nz;z++) {
    if (RE(f,z) > 0) meanr += complex<double>(RE(f,z), IM(f,z));
    else meanr += complex<double>(-RE(f,z), -IM(f,z));
  }
  complex<double> meani=0;
  for (int z=0;z<nz;z++) {
    if (IM(f,z) > 0) meani += complex<double>(RE(f,z), IM(f,z));
    else meani += complex<double>(-RE(f,z), -IM(f,z));
  }
  if (abs(mean) > abs(meanr) && abs(mean) > abs(meani)) return arg(mean);
  if (abs(meanr) > abs(meani)) return arg(meanr);
  if (abs(meani) > 0.0) return arg(meani);
  return 0;
}

static void output_complex_slice(double *f[2], int nz, const char *name) {
  FILE *out = fopen(name, "w");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  double phase = get_phase(f, nz);
  double c = cos(phase), s = sin(phase);
  for (int z=0;z<nz;z++) {
    fprintf(out, "%d\t0\t%lg\n", z, c*RE(f,z)+s*IM(f,z));
  }
  fclose(out);
}

static void output_slice(const double *f, int nz, const char *name) {
  FILE *out = fopen(name, "w");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  for (int z=0;z<nz;z++) {
    fprintf(out, "%d\t0\t%lg\n", z, MA(f,z));
  }
  fclose(out);
}

void mat_1d::output_sigma_slice(const char *filename) {
  double *sigma = new double[nz+1];
  for (int z=0;z<nz;z++) MA(sigma,z) = 0.0;
  if (npmlz) {
    for (int z=0;z<npmlz;z++) MA(sigma,z) = Czex[z];
    for (int z=0;z<npmlz;z++) MA(sigma,nz-z-1) = Czex[z];
  }
  output_slice(sigma, nz, filename);
  delete[] sigma;
}

void mat_1d::output_slices(const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  if (!n) {
    printf("Allocation failure!\n");
    exit(1);
  }
  snprintf(n, buflen, "%s/%sepsilon.sli", outdir, nname);
  output_slice(eps, nz, n);
  snprintf(n, buflen, "%s/%ssigma.sli", outdir, nname);
  output_sigma_slice(n);
  free(n);
}

void fields_1d::output_real_imaginary_slices(const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) {
    printf("Allocation failure!\n");
    exit(1);
  }
  char *r_or_i = "-re";
  for (int cmp=0;cmp<2;cmp++) {
    if (a == 1) snprintf(n, buflen, "%s/%shy%s-%08.0f.sli", outdir, nname, r_or_i, t*inva*c);
    else snprintf(n, buflen, "%s/%shy%s-%09.2f.sli", outdir, nname, r_or_i, t*inva*c);
    output_slice(hy[cmp], nz, n);
    if (a == 1) snprintf(n, buflen, "%s/%sex%s-%08.0f.sli", outdir, nname, r_or_i, t*inva*c);
    else snprintf(n, buflen, "%s/%sex%s-%09.2f.sli", outdir, nname, r_or_i, t*inva*c);
    output_slice(ex[cmp], nz, n);
    r_or_i = "-im";
  }

  free(n);
}

void fields_1d::output_slices(const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) {
    printf("Allocation failure!\n");
    exit(1);
  }
  char time_step_string[buflen];
  if (a == 1)
    snprintf(time_step_string, buflen, "%08.0f", time());
  else
    snprintf(time_step_string, buflen, "%09.2f", time());

  {
    polarization_1d *p = olpol;
    int polnum = 0;
    while (p) {
      snprintf(n, buflen, "%s/%sp%d-%s.sli", outdir, nname, polnum, time_step_string);
      output_complex_slice(p->Px, nz, n);
      polnum++;
      p = p->next;
    }
  }
  snprintf(n, buflen, "%s/%shy-%s.sli", outdir, nname, time_step_string);
  output_complex_slice(hy, nz, n);
  snprintf(n, buflen, "%s/%sex-%s.sli", outdir, nname, time_step_string);
  output_complex_slice(ex, nz, n);
  
  if (new_ma) {
    snprintf(n, buflen, "%s/%sepsilon-%s.sli", outdir, nname, time_step_string);
    output_slice(ma->eps, nz, n);
  }
  
  free(n);
}

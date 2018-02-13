/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
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

/* Check of Green's functions (analytical vs. numerical)
   and near-to-far-field transformation. */

#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
using namespace meep;
using namespace std;

double two(const vec &) { return 2.0; }

const int EHcomp[10] = {0,1,0,1,2, 3,4,3,4,5};

int check_2d_3d(ndim dim, const double xmax, double a, component c0) {
  const double dpml = 1;
  if (dim != D2 && dim != D3) abort("2d or 3d required");
  grid_volume gv = dim == D2 ? vol2d(xmax + 2*dpml,xmax + 2*dpml,a) : vol3d(xmax + 2*dpml,xmax + 2*dpml,xmax + 2*dpml,a);
  gv.center_origin();

  if (!gv.has_field(c0)) return 1;
  master_printf("TESTING %s AT RESOLUTION %g FOR %s SOURCE...\n",
                dim == D2 ? "2D" : "3D", a, component_name(c0));

  structure s(gv, two, pml(dpml));
  fields f(&s);
  double w = 0.30;
  continuous_src_time src(w);
  f.add_point_source(c0, src, zero_vec(dim));
  f.solve_cw(1e-6);

  FOR_E_AND_H(c) if (gv.has_field(c)) {
          const int N = 20;
          double dx = 0.75 * (xmax/4) / N;
          complex<double> F[N], F0[N], EH[6];
          double diff = 0.0, dot0 = 0.0, dot = 0.0;
          complex<double> phase = polar(1.0, (4*w*f.dt)*pi);
          vec x0 = zero_vec(dim);
          for (int i = 0; i < N; ++i) {
              double s = xmax/4 + dx*i;
              vec x = dim == D2 ? vec(s, 0.5*s) : vec(s, 0.5*s, 0.3*s);
              F[i] = f.get_field(c, x) * phase;
              (dim == D2 ? green2d : green3d)(EH, x, w, 2.0, 1.0, x0, c0, 1.0);
              F0[i] = EH[EHcomp[c]];
              double d = abs(F0[i] - F[i]);
              double f = abs(F[i]);
              double f0 = abs(F0[i]);
              diff += d*d;
              dot += f*f;
              dot0 += f0*f0;
          }
          if (dot0 == 0) continue; /* zero field component */
          double relerr = sqrt(diff) / sqrt(dot0);

          master_printf("  GREEN: %s -> %s, resolution %g: relerr = %g\n",
                        component_name(c0), component_name(c), a, relerr);

          if (relerr > 0.05 * 30/a) {
              for (int i = 0; i < N; ++i)
                  master_printf("%g, %g,%g, %g,%g\n",
                                xmax/4 + dx*i,
                                real(F[i]), imag(F[i]),
                                real(F0[i]), imag(F0[i]));
              return 0;
          }
  }

  const double L = xmax/4;
  volume_list vl =
     dim == D2 ?
     volume_list(                volume(vec(-L,+L),vec(+L,+L)), Sy, 1.0,
                 new volume_list(volume(vec(+L,+L),vec(+L,-L)), Sx, 1.0,
                 new volume_list(volume(vec(-L,-L),vec(+L,-L)), Sy, -1.0,
                 new volume_list(volume(vec(-L,-L),vec(-L,+L)), Sx, -1.0)))) :
    volume_list(                volume(vec(+L,-L,-L),vec(+L,+L,+L)), Sx, +1.,
                new volume_list(volume(vec(-L,-L,-L),vec(-L,+L,+L)), Sx, -1.,
                new volume_list(volume(vec(-L,+L,-L),vec(+L,+L,+L)), Sy, +1.,
                new volume_list(volume(vec(-L,-L,-L),vec(+L,-L,+L)), Sy, -1.,
                new volume_list(volume(vec(-L,-L,+L),vec(+L,+L,+L)), Sz, +1.,
                new volume_list(volume(vec(-L,-L,-L),vec(+L,+L,-L)), Sz, -1.
                                ))))));

  dft_near2far n2f = f.add_dft_near2far(&vl, w, w, 1);
  f.update_dfts();
  n2f.scale_dfts(sqrt(2*pi)/f.dt); // cancel time-integration factor

  FOR_E_AND_H(c) if (gv.has_field(c)) {
          const int N = 20;
          double dx = 0.75 * (xmax/4) / N;
          complex<double> F[N], F0[N], EH_[6], EH[6];
          double diff = 0.0, dot = 0.0;
          complex<double> phase = polar(1.0, (4*w*f.dt)*pi);
          vec x0 = zero_vec(dim);
          for (int i = 0; i < N; ++i) {
              double s = xmax + dx*i;
              vec x = dim == D2 ? vec(s, 0.5*s) : vec(s, 0.5*s, 0.3*s);
              n2f.farfield_lowlevel(EH_, x);
              sum_to_all(EH_, EH, 6);
              F[i] = EH[EHcomp[c]] * phase;
              (dim == D2 ? green2d : green3d)(EH, x, w, 2.0, 1.0, x0, c0, 1.0);
              F0[i] = EH[EHcomp[c]];
              double d = abs(F0[i] - F[i]);
              double f0 = abs(F0[i]);
              diff += d*d;
              dot += f0*f0;
          }
          if (dot == 0) continue; /* zero field component */
          double relerr = sqrt(diff) / sqrt(dot);

          master_printf("  NEAR2FAR: %s -> %s, resolution %g: relerr = %g\n",
                        component_name(c0), component_name(c), a, relerr);

          if (relerr > 0.05 * 30/a) {
              for (int i = 0; i < N; ++i)
                  master_printf("%g, %g,%g, %g,%g\n",
                                xmax + dx*i,
                                real(F[i]), imag(F[i]),
                                real(F0[i]), imag(F0[i]));
              return 0;
          }
  }

  return 1;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);

  const double a2d = argc > 1 ? atof(argv[1]) : 20, a3d = argc > 1 ? a2d : 10;

  if (!check_2d_3d(D3, 4, a3d, Ez)) return 1;
  if (!check_2d_3d(D3, 4, a3d, Hz)) return 1;
#ifdef HAVE_LIBGSL
  FOR_E_AND_H(c0) if (!check_2d_3d(D2, 8, a2d, c0)) return 1;
#endif

  return 0;
}

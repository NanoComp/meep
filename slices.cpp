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

static double get_phase(component c, double *f[2],
                        const volume &v, const volume &what) {
  complex<double> mean=0;
  for (int i=0;i<v.ntot();i++)
    mean += v.dv(c,i)*complex<double>(f[0][i],f[1][i]);
  complex<double> meanr=0;
  for (int i=0;i<v.ntot();i++) {
    complex<double> val = v.dv(c,i)*complex<double>(f[0][i],f[1][i]);
    if (f[0][i] > 0) meanr += val;
    else meanr -= val;
  }
  complex<double> meani=0;
  for (int i=0;i<v.ntot();i++) {
    complex<double> val = v.dv(c,i)*complex<double>(f[0][i],f[1][i]);
    if (f[1][i] > 0) meani += val;
    else meani -= val;
  }
  if (abs(mean) > abs(meanr) && abs(mean) > abs(meani)) return arg(mean);
  if (abs(meanr) > abs(meani)) return arg(meanr);
  if (abs(meani) > 0.0) return arg(meani);
  return 0;
}

static void output_complex_slice(component m, double *f[2], const volume &v,
                                 const volume &what, const char *name) {
  if (!f[0] || ! f[1]) return; // Field doesn't exist...
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  double phase = get_phase(m, f, v, what);
  double c = cos(phase), s = sin(phase);
  for (int i=0;i<v.ntot();i++) {
    if (what.contains(v.loc(m,i))) {
      v.loc(m,i).print(out);
      fprintf(out, "\t%lg\n", c*f[0][i]+s*f[1][i]);
    }
  }
  fclose(out);
}

static void output_slice(component m, const double *f, const volume &v,
                         const volume &what, const char *name) {
  if (!f) return; // Field doesn't exist...
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  for (int i=0;i<v.ntot();i++) {
    if (what.contains(v.loc(m,i))) {
      v.loc(m,i).print(out);
      fprintf(out, "\t%.17lg\n", f[i]);
    }
  }
  fclose(out);
}

const double default_eps_resolution = 20.0;
const double default_eps_size = 1000.0;

static void eps_header(double xmin, double ymin, double xmax, double ymax,
                       double fmax, double a, FILE *out, const char *name) {
  fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
  double size = xmax - xmin + ymax - ymin;
  fprintf(out, "%%%%BoundingBox: %lg %lg %lg %lg\n",
          xmin*default_eps_size/size, ymin*default_eps_size/size,
          xmax*default_eps_size/size, ymax*default_eps_size/size);
  fprintf(out, "gsave\n");
  fprintf(out, "%lg %lg scale\n", default_eps_size/size, default_eps_size/size);
  fprintf(out, "/Times-Roman findfont 20 scalefont setfont\n");
  fprintf(out, "newpath 140 280 moveto (%s) show\n", name);
  fprintf(out, "/max %lg def\n", fmax);
  const double dx = 1.0/min(a, default_eps_resolution);
  fprintf(out, "/dx %lg def\n", dx);
  fprintf(out, "/hdx %lg def\n", dx*0.5);
  fprintf(out, "dx 10 div setlinewidth\n");
  fprintf(out, "1 setlinecap\n");
  fprintf(out, "/P {\n\
    max div\n\
    dup 0 lt {\n\
        1 add\n\
        dup 1\n\
    }{\n\
        neg 1 add\n\
        dup 1 3 1 roll\n\
    } ifelse\n\
    setrgbcolor\n\
    newpath\n\
    moveto\n");
  fprintf(out, "    %lg %lg rmoveto\n", dx*0.5, dx*0.5);
  fprintf(out, "    0 %lg rlineto\n", -dx);
  fprintf(out, "    %lg 0 rlineto\n", -dx);
  fprintf(out, "    0 dx rlineto\n", dx);
  fprintf(out, "    dx 0 rlineto\n", dx);
  fprintf(out, "    gsave\n\
    fill\n\
    grestore\n");
  fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  fprintf(out, "    stroke\n\
} def\n\
/LV {\n\
    0 0 0 setrgbcolor\n\
    moveto\n");
  fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  fprintf(out, "    0 %lg rmoveto\n", dx*0.5);
  fprintf(out, "    0 %lg rlineto\n", -dx);
  fprintf(out, "    stroke\n\
} def\n\
/LH {\n");
  fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  fprintf(out, "    0 0 0 setrgbcolor\n\
    moveto\n");
  fprintf(out, "    %lg 0 rmoveto\n", dx*0.5);
  fprintf(out, "    %lg 0 rlineto\n", -dx);
  fprintf(out, "    stroke\n\
} def\n");
  fprintf(out, "    /DV { [0 %lg] 0 setdash LV } def\n", dx/4);
  fprintf(out, "    /DH { [0 %lg] 0 setdash LH } def\n", dx/4);
  fprintf(out, "    /D { moveto\n\
    0 1 0 setrgbcolor [0 %lg] 0 setdash\n\
    lineto stroke } def\n", dx);
}

static void eps_1d_header(double xmin, double ymin, double xmax, double ymax,
                          double fmax, double a, FILE *out, const char *name) {
  fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
  const double size = xmax - xmin;
  const double fsize = (5.0 < 0.2*size)?0.2*size:5.0;
  const double dx = 1.0/a;
  fprintf(out, "%%%%BoundingBox: %lg %lg %lg %lg\n",
          xmin*default_eps_size/size, -default_eps_size*0.5*fsize/size,
          xmax*default_eps_size/size, default_eps_size*0.5*fsize/size);
  fprintf(out, "gsave\n");
  fprintf(out, "%lg 0 moveto %lg 0 lineto %lg setlinewidth stroke\n",
          xmin*default_eps_size/size, xmax*default_eps_size/size, dx*0.1); 
  fprintf(out, "%lg %lg scale\n", default_eps_size/size, default_eps_size/size);
  fprintf(out, "/height %lg def\n", (default_eps_size*0.5*fsize)/size);
  fprintf(out, "/Times-Roman findfont 12 scalefont setfont\n");
  fprintf(out, "newpath 220 height 0.75 mul moveto (%s) show\n", name);
  fprintf(out, "/max %lg def\n", fmax);
  fprintf(out, "/fscale %lg def\n", 2.2*fmax/fsize);
  fprintf(out, "/dotrad %lg def\n", dx);
  fprintf(out, "/dx %lg def\n", dx);
  fprintf(out, "/hdx %lg def\n", dx*0.5);
  fprintf(out, "dx 10 div setlinewidth\n");
  fprintf(out, "1 setlinecap\n");
  fprintf(out, "/P {\n\
    fscale div exch pop \n\
    newpath dotrad 0 360 arc fill \n\
} def\n");
  fprintf(out, "/LV {\n\
    0.8 0.8 0 setrgbcolor\n\
    pop dup height moveto height neg lineto\n");
  fprintf(out, "    %lg setlinewidth\n", dx*0.5);
  fprintf(out, "    stroke\n\
} def\n\
/LH {\n");
  fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  fprintf(out, "    0 0 0 setrgbcolor\n\
    moveto\n");
  fprintf(out, "    %lg 0 rmoveto\n", 10*dx*0.5);
  fprintf(out, "    %lg 0 rlineto\n", -10*dx);
  fprintf(out, "    stroke\n\
} def\n");
  fprintf(out, "    /DV { [0 %lg] 0 setdash LV } def\n", dx/4);
  fprintf(out, "    /D { moveto\n\
    0 1 0 setrgbcolor [0 %lg] 0 setdash\n\
    lineto stroke } def\n", dx);
}

static void eps_trailer(FILE *out) {
  fprintf(out, "grestore\n");
  fprintf(out, "showpage\n");
  fprintf(out, "%%%%Trailer\n");
  fprintf(out, "%%%%EOF\n");
}

static void eps_dotted(FILE *out, component m, const double *f, const volume &v,
                       const volume &what) {
  if (!f) return; // Field doesn't exist...
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i)))
      switch (v.dim) {
      case Dcyl:
        {
          ivec next = v.iloc(m,i)+ivec(2,0);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDH\n", v[next].z(), v[next].r() - 0.5/v.a);
          next = v.iloc(m,i)+ivec(0,2);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDV\n", v[next].z() - 0.5/v.a, v[next].r());
          break;
        }
      case D2:
        {
          ivec next = v.iloc(m,i)+ivec2d(2,0);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDH\n", v[next].x() - 0.5/v.a, v[next].y());
          next = v.iloc(m,i)+ivec2d(0,2);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDV\n", v[next].x(), v[next].y() - 0.5/v.a);
          break;
        }
      }
}

void fields::outline_chunks(const char *name) {
  if (v.dim == D1) return;
  if (my_rank()) return;
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for epsilon border output.\n", name);
    return;
  }
  for (int i=0;i<num_chunks;i++) {
    double x1, y1, x2, y2, x3, y3, x4, y4;
    switch (v.dim) {
    case Dcyl:
      x1 = chunks[i]->v.zmin(); y1 = chunks[i]->v.rmin();
      x2 = chunks[i]->v.zmin(); y2 = chunks[i]->v.rmax();
      x3 = chunks[i]->v.zmax(); y3 = chunks[i]->v.rmax();
      x4 = chunks[i]->v.zmin(); y4 = chunks[i]->v.rmin();
    break;
    case D2:
      x1 = chunks[i]->v.xmin(); y1 = chunks[i]->v.ymin();
      x2 = chunks[i]->v.xmin(); y2 = chunks[i]->v.ymax();
      x3 = chunks[i]->v.xmax(); y3 = chunks[i]->v.ymax();
      x4 = chunks[i]->v.xmax(); y4 = chunks[i]->v.ymin();
    }
    fprintf(out, "%lg\t%lg\t%lg\t%lg\tD\n", x1, y1, x2, y2);
    fprintf(out, "%lg\t%lg\t%lg\t%lg\tD\n", x2, y2, x3, y3);
    fprintf(out, "%lg\t%lg\t%lg\t%lg\tD\n", x3, y3, x4, y4);
    fprintf(out, "%lg\t%lg\t%lg\t%lg\tD\n", x4, y4, x1, y1);
  }
  fclose(out);
}

static void eps_outline(component m, const double *f,
                        const volume &v, const volume &what,
                        symmetry S, int symnum, const char *name) {
  if (!f) return; // Field doesn't exist...
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for epsilon border output.\n", name);
    return;
  }
  for (int i=0;i<v.ntot();i++) {
    const vec here = S.transform(v.loc(m,i),symnum);
    if (what.contains(here))
      switch (v.dim) {
      case Dcyl: {
        ivec next = v.iloc(m,i)+ivec(2,0);
        vec nextrot = v[S.transform(next - ivec(1,0),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          fprintf(out, "%lg\t%lg\tLH\n", nextrot.z(), nextrot.r());
        next = v.iloc(m,i)+ivec(0,2);
        nextrot = v[S.transform(next - ivec(0,1),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          fprintf(out, "%lg\t%lg\tLV\n", nextrot.z(), nextrot.r());
        break;
      }
      case D2: {
        ivec next = v.iloc(m,i)+ivec2d(0,2);
        vec nextrot = v[(S.transform(next - ivec2d(0,1),symnum))];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLH\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)-ivec2d(0,2);
        nextrot = v[S.transform(next + ivec2d(0,1),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLH\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)+ivec2d(2,0);
        nextrot = v[S.transform(next - ivec2d(1,0),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLV\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)-ivec2d(2,0);
        nextrot = v[S.transform(next + ivec2d(1,0),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLV\n", nextrot.x(), nextrot.y());
        break;
      }
      case D1: {
        const ivec next = v.iloc(m,i)+ivec(2);
        const vec nextrot = v[S.transform(next - ivec(1),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          fprintf(out, "%lg\t%lg\tLV\n", nextrot.z(), 0.0);
        break;
      }
      }
  }
  fclose(out);
}

static void output_complex_eps_body(component m, double *f[2], const volume &v,
                                    symmetry S, int symnum,
                                    const volume &what, const char *name) {
  if (!f[0] || !f[1]) return; // Field doesn't exist...
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  const complex<double> ph = S.phase_shift(m, symnum);
  for (int i=0;i<v.ntot();i++) {
    const vec here = S.transform(v.loc(m,i),symnum);
    if (what.contains(here)) {
      double x = 0, y = 0;
      switch (v.dim) {
      case Dcyl: x = here.z(); y = here.r(); break;
      case D1: x = here.z(); break;
      case D2: x = here.x(); y = here.y(); break;
      }
      if (f[1]) fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y,
                        real(ph)*f[0][i] - imag(ph)*f[1][i]);
      else fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y, real(ph)*f[0][i]);
    }
  }
  fclose(out);
}

void fields_chunk::output_eps_body(component c, const symmetry &S, int sn,
                                   const volume &what, const char *n) {
  if (f[c][0]) {
    if (v.a <= default_eps_resolution) {
      output_complex_eps_body(c, f[c], v, S, sn, what, n);
    } else {
      const complex<double> ph = S.phase_shift(c, sn);
      FILE *out = fopen(n, "a");
      if (!out) {
        printf("Unable to open file '%s' for slice output.\n", n);
        return;
      }
      const int n = v.ntot_at_resolution(default_eps_resolution);
      for (int i=0;i<n;i++) {
        const vec here = v.loc_at_resolution(i,default_eps_resolution);
        const vec there = S.transform(here, sn);
        if (what.contains(here)) {
          double x = 0, y = 0;
          switch (v.dim) {
          case Dcyl: x = there.z(); y = there.r(); break;
          case D1: x = there.z(); break;
          case D2: x = there.x(); y = there.y(); break;
          }
          complex<double> val[8];
          for (int j=0;j<8;j++) val[j] = 0.0;
          interpolate_field_private(c, here, val, ph);
          double fhere = 0.0;
          for (int j=0;j<8;j++) fhere += real(val[j]);
          fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y, fhere);
        }
      }
      fclose(out);
    }
  }
}

static void output_complex_eps_header(component m, double fmax, const volume &v,
                                      const volume &what, const char *name,
                                      component om = Hx) {
  double xmin = 1e300, ymin = 1e300, xmax = -1e300, ymax = -1e300;
  switch (v.dim) {
  case Dcyl:
    xmin = what.origin.z();
    xmax = what.origin.z() + what.nz()*what.inva;
    ymin = what.origin.r();
    ymax = what.origin.r() + what.nr()*what.inva;
    break;
  case D2:
    xmin = what.origin.x();
    xmax = what.origin.x() + what.nx()*what.inva;
    ymin = what.origin.y();
    ymax = what.origin.y() + what.ny()*what.inva;
    break;
  case D1:
    ymin = -v.inva*0.5;
    ymax = v.inva*0.5;
    xmin = what.origin.z();
    xmax = what.origin.z() + what.nz()*what.inva;
    break;
  }
  if (ymax == ymin) ymax = ymin + 1.0/v.a;
  if (fmax == 0.0) fmax = 0.0001;
  if (v.dim == D1) {
    // Make a 1D line plot!
    FILE *out = fopen(name, "w");
    if (!out) {
      printf("Unable to open file '%s' for slice output.\n", name);
      return;
    }
    eps_1d_header(xmin, ymin, xmax, ymax, fmax, v.a, out, name);
    fclose(out);
  } else {
    // Make a 2D color plot!
    FILE *out = fopen(name, "w");
    if (!out) {
      printf("Unable to open file '%s' for slice output.\n", name);
      return;
    }
    eps_header(xmin, ymin, xmax, ymax, fmax, v.a, out, name);
    fclose(out);
  }
}

static void output_complex_eps_tail(const char *name) {
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  eps_trailer(out);
  fclose(out);
}

void mat::output_slices(const char *name) const {
  output_slices(v,name);
}
void mat::output_slices(const volume &what, const char *name) const {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = new char[buflen];
  if (!n) abort("Allocation failure!\n");
  snprintf(n, buflen, "%s/%sepsilon.sli", outdir, nname);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      output_slice(v.eps_component(), chunks[i]->eps, chunks[i]->v, what, n);
  //snprintf(n, buflen, "%s/%ssigma.sli", outdir, nname);
  //output_sigma_slice(n);
  delete[] n;
}

void fields::output_real_imaginary_slices(const char *name) {
  output_real_imaginary_slices(v,name);
}

void fields::output_real_imaginary_slices(const volume &what,
                                          const char *name) {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) abort("Allocation failure!\n");
  char *r_or_i = "-re";
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c)) {
        snprintf(n, buflen, "%s/%s%s%s-%09.2f.sli",
                 outdir, nname, component_name((component)c),
                 r_or_i, time());
        for (int i=0;i<num_chunks;i++)
          if (chunks[i]->is_mine())
            output_slice((component)c, chunks[i]->f[c][cmp], chunks[i]->v, what, n);
      }
    r_or_i = "-im";
  }

  free(n);
}

void fields::output_slices(const char *name) {
  output_slices(v, name);
}
void fields::output_slices(const volume &what, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  if (a == 1)
    snprintf(time_step_string, buflen, "%08.0f", time());
  else
    snprintf(time_step_string, buflen, "%09.2f", time());
  for (int c=0;c<10;c++)
    if (v.has_field((component)c)) {
      snprintf(n, buflen, "%s/%s%s-%s.sli", outdir, nname,
               component_name((component)c), time_step_string);
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          output_complex_slice((component)c, chunks[i]->f[c], chunks[i]->v, what, n);
    }
  //if (new_ma) {
  //  snprintf(n, buflen, "%s/%sepsilon-%s.sli", outdir, nname, time_step_string);
  //  output_slice(v.eps_component(), ma->eps, v, what, n);
  //}
  free(n);
  finished_working();
}

void fields::eps_slices(const char *name) {
  if (v.dim == Dcyl || v.dim == D1 || v.dim == D2)
    eps_slices(user_volume, name);
}

void fields::eps_slices(const volume &what, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());
  for (int c=0;c<10;c++)
    if (v.has_field((component)c)) {
      snprintf(n, buflen, "%s/%s%s-%s.eps", outdir, nname,
               component_name((component)c), time_step_string);
      const double fmax = maxfieldmag_to_master((component)c);
      if (am_master())
        output_complex_eps_header((component)c, fmax,
                                  user_volume, what,
                                  n, v.eps_component());
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            for (int otherc=0;otherc<10;otherc++)
              if (S.transform((component)otherc,sn) == c)
                chunks[i]->output_eps_body((component)otherc,
                                           S, sn, what, n);
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            for (int otherc=0;otherc<10;otherc++)
              if (S.transform((component)otherc,sn) == c)
                eps_outline(v.eps_component(), chunks[i]->ma->eps,
                            chunks[i]->v, what, S, sn, n);
      all_wait();
      outline_chunks(n);
      all_wait();
      if (am_master()) output_complex_eps_tail(n);
    }
  free(n);
  finished_working();
}

double fields::maxfieldmag_to_master(component c) const {
  double themax = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      themax = max(themax,chunks[i]->maxfieldmag(c));
  return max_to_master(themax);
}

double fields_chunk::maxfieldmag(component c) const {
  double themax = 0.0;
  if (f[c][0])
    for (int i=0;i<v.ntot();i++) {
      double norm = 0.0;
      DOCMP norm += f[c][cmp][i]*f[c][cmp][i];
      themax = max(themax,sqrt(norm));
    }
  return themax;
}

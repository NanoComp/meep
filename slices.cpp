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

static void eps_header(double xmin, double ymin, double xmax, double ymax,
                       double fmax, double a, FILE *out, const char *name) {
  fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
  double size = xmax - xmin + ymax - ymin;
  fprintf(out, "%%%%BoundingBox: %lg %lg %lg %lg\n",
          xmin*500/size, ymin*500/size, xmax*500/size, ymax*500/size);
  fprintf(out, "gsave\n");
  fprintf(out, "%lg %lg scale\n", 500/size, 500/size);
  fprintf(out, "/Times-Roman findfont 20 scalefont setfont\n");
  fprintf(out, "newpath 140 280 moveto (%s) show\n", name);
  fprintf(out, "/max %lg def\n", fmax);
  const double dx = 1.0/a;
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
}

static void eps_1d_header(double xmin, double ymin, double xmax, double ymax,
                          double fmax, double a, FILE *out, const char *name) {
  fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
  const double size = xmax - xmin;
  const double fsize = (5.0 < 0.2*size)?0.2*size:5.0;
  const double dx = 1.0/a;
  fprintf(out, "%%%%BoundingBox: %lg %lg %lg %lg\n",
          xmin*500/size, -250*fsize/size, xmax*500/size, 250*fsize/size);
  fprintf(out, "gsave\n");
  fprintf(out, "%lg 0 moveto %lg 0 lineto %lg setlinewidth stroke\n",
          xmin*500/size, xmax*500/size, dx*0.1); 
  fprintf(out, "%lg %lg scale\n", 500/size, 500/size);
  fprintf(out, "/height %lg def\n", (250*fsize)/size);
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
}

static void eps_trailer(FILE *out) {
  fprintf(out, "grestore\n");
  fprintf(out, "%%%%Trailer\n");
  fprintf(out, "%%%%EOF\n");
}

static void eps_dotted(FILE *out, component m, const double *f, const volume &v,
                       const volume &what) {
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i)))
      switch (v.dim) {
      case dcyl:
        {
          vec next = v.loc(m,i)+v.dr();
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDH\n", next.z(), next.r() - 0.5/v.a);
          next = v.loc(m,i)+v.dz();
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDV\n", next.z() - 0.5/v.a, next.r());
          break;
        }
      case d2:
        {
          vec next = v.loc(m,i)+v.dx();
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDH\n", next.x() - 0.5/v.a, next.y());
          next = v.loc(m,i)+v.dy();
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            fprintf(out, "%lg\t%lg\tDV\n", next.x(), next.y() - 0.5/v.a);
          break;
        }
      }
}

static void eps_outline(FILE *out, component m, const double *f, const volume &v,
                        const volume &what) {
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i)))
      switch (v.dim) {
      case dcyl: {
        vec next = v.loc(m,i)+v.dr();
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          fprintf(out, "%lg\t%lg\tLH\n", next.z(), next.r() - 0.5/v.a);
        next = v.loc(m,i)+v.dz();
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          fprintf(out, "%lg\t%lg\tLV\n", next.z() - 0.5/v.a, next.r());
        break;
      }
      case d2: {
        vec next = v.loc(m,i)+v.dy();
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLH\n", next.x(), next.y() - 0.5/v.a);
        next = v.loc(m,i)-v.dy();
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLH\n", next.x(), next.y() + 0.5/v.a);
        next = v.loc(m,i)+v.dx();
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLV\n", next.x() - 0.5/v.a, next.y());
        next = v.loc(m,i)-v.dx();
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            fprintf(out, "%lg\t%lg\tLV\n", next.x() + 0.5/v.a, next.y());
        break;
      }
      case d1: {
        const vec next = v.loc(m,i)+v.dz();
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          fprintf(out, "%lg\t%lg\tLV\n", next.z() - 0.5/v.a, 0.0);
        break;
      }
      }
}

static void output_eps(component m, const double *f, const volume &v,
                       const volume &what, const char *name,
                       component om = Hx, const double *overlay = NULL,
                       const double *dashed = NULL) {
  double xmin = 1e300, ymin = 1e300, xmax = -1e300, ymax = -1e300, fmax = 0.0;
  switch (v.dim) {
  case dcyl:
    xmin = what.origin.z();
    xmax = what.origin.z() + what.nz()*what.inva;
    ymin = ymax = what.origin.r();
    ymax = what.origin.r() + what.nr()*what.inva;
    break;
  case d1:
    ymin = -v.inva*0.5;
    ymax = v.inva*0.5;
    xmin = what.origin.z();
    xmax = what.origin.z() + what.nz()*what.inva;
    break;
  case d2:
    xmin = what.origin.x();
    xmax = what.origin.x() + what.nx()*what.inva;
    ymin = ymax = what.origin.y();
    ymax = what.origin.y() + what.ny()*what.inva;
  }
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i)))
      fmax = max(fmax, fabs(f[i]));
  if (fmax == 0.0) fmax = 0.0001;
  if (v.dim == d1) {
    // Make a 1D line plot!
    grace g(name);
    int skipnum = 1 + v.ntot()/1000;
    for (int i=0;i<v.ntot();i+=skipnum)
      if (what.contains(v.loc(m,i)))
        g.output_point(v.loc(m,i).z(), f[i]);
    //for (int i=1;i<v.ntot()-1;i++)
    //  if (what.contains(v.loc(m,i)) && f[i] > f[i-1] && f[i] > f[i+1])
    //    g.output_point(v.loc(m,i).z(), f[i]);
  } else {
    // Make a 2D grayscale plot!
    FILE *out = fopen(name, "w");
    if (!out) {
      printf("Unable to open file '%s' for slice output.\n", name);
      return;
    }
    eps_header(xmin, ymin, xmax, ymax, fmax, v.a, out, name);
    for (int i=0;i<v.ntot();i++) {
      if (what.contains(v.loc(m,i))) {
        double x = 0, y = 0;
        switch (v.dim) {
        case dcyl: x = v.loc(m,i).z(); y = v.loc(m,i).r(); break;
        case d1: x = v.loc(m,i).z(); break;
        case d2: x = v.loc(m,i).x(); y = v.loc(m,i).y(); break;
        }
        fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y, f[i]);
      }
    }
    if (overlay) eps_outline(out, om, overlay, v, what);
    if (dashed) eps_dotted(out, om, dashed, v, what);
    eps_trailer(out);
    fclose(out);
  }
}

static void output_complex_eps_body(component m, double *f[2], const volume &v,
                                    const volume &what, const char *name,
                                    component om = Hx, const double *overlay = NULL,
                                    const double *dashed = NULL) {
  FILE *out = fopen(name, "a");
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", name);
    return;
  }
  for (int i=0;i<v.ntot();i++) {
    if (what.contains(v.loc(m,i))) {
      double x = 0, y = 0;
      switch (v.dim) {
      case dcyl: x = v.loc(m,i).z(); y = v.loc(m,i).r(); break;
      case d1: x = v.loc(m,i).z(); break;
      case d2: x = v.loc(m,i).x(); y = v.loc(m,i).y(); break;
      }
      fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y, f[0][i]);
    }
  }
  if (overlay) eps_outline(out, om, overlay, v, what);
  if (dashed) eps_dotted(out, om, dashed, v, what);
  fclose(out);
}

static void output_complex_eps_header(component m, double fmax, const volume &v,
                                      const volume &what, const char *name,
                                      component om = Hx) {
  double xmin = 1e300, ymin = 1e300, xmax = -1e300, ymax = -1e300;
  switch (v.dim) {
  case dcyl:
    xmin = what.origin.z();
    xmax = what.origin.z() + what.nz()*what.inva;
    ymin = what.origin.r();
    ymax = what.origin.r() + what.nr()*what.inva;
    break;
  case d2:
    xmin = what.origin.x();
    xmax = what.origin.x() + what.nx()*what.inva;
    ymin = what.origin.y();
    ymax = what.origin.y() + what.ny()*what.inva;
    break;
  case d1:
    ymin = -v.inva*0.5;
    ymax = v.inva*0.5;
    xmin = what.origin.z();
    xmax = what.origin.z() + what.nz()*what.inva;
    break;
  }
  if (ymax == ymin) ymax = ymin + 1.0/v.a;
  if (fmax == 0.0) fmax = 0.0001;
  if (v.dim == d1) {
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

static void output_complex_eps_tail(component m, const volume &v,
                                    const volume &what, const char *name) {
  if (v.dim == d1 && 0) {
    // Make a 1D line plot!
  } else {
    // Make a 2D color plot!
    FILE *out = fopen(name, "a");
    if (!out) {
      printf("Unable to open file '%s' for slice output.\n", name);
      return;
    }
    eps_trailer(out);
    fclose(out);
  }
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

void fields::output_real_imaginary_slices(const char *name) const {
  output_real_imaginary_slices(v,name);
}

void fields::output_real_imaginary_slices(const volume &what, const char *name) const {
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

void fields::output_slices(const char *name) const {
  output_slices(v, name);
}
void fields::output_slices(const volume &what, const char *name) const {
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
  {
    //polarization *p = pol;
    //int polnum = 0;
    //while (p) {
    //  for (int c=0;c<10;c++)
    //    if (v.has_field((component)c) && is_electric((component)c)) {
    //      snprintf(n, buflen, "%s/%sp%d%s-%s.sli", outdir, nname, polnum,
    //               component_name((component)c), time_step_string);
    //      output_complex_slice((component)c, p->P[c], v, what, n);
    //    }
    //  polnum++;
    //  p = p->next;
    //}
  }
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
}

void fields::eps_slices(const char *name) const {
  if (v.dim == dcyl || v.dim == d1 || v.dim == d2)
    eps_slices(v, name);
}

void fields::eps_slices(const volume &what, const char *name) const {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  int i;
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());
  /*{
    polarization *p = pol;
    int polnum = 0;
    while (p) {
      for (int c=0;c<10;c++)
        if (v.has_field((component)c) && is_electric((component)c)) {
          snprintf(n, buflen, "%s/%sp%d%s-%s.eps", outdir, nname, polnum,
                   component_name((component)c), time_step_string);
          output_complex_eps((component)c, p->P[c], v, what, n);
          snprintf(n, buflen, "%s/%senergy%d%s-%s.eps", outdir, nname, polnum,
                   component_name((component)c), time_step_string);
          output_eps((component)c, p->energy[c], v, what, n);
        }
      polnum++;
      p = p->next;
    }
  }*/
  for (int c=0;c<10;c++)
    if (v.has_field((component)c)) {
      snprintf(n, buflen, "%s/%s%s-%s.eps", outdir, nname,
               component_name((component)c), time_step_string);
      const double fmax = maxfieldmag_to_master((component)c);
      if (am_master())
        output_complex_eps_header((component)c, fmax,
                                  chunks[0]->v, what, n, v.eps_component());
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          output_complex_eps_body((component)c, chunks[i]->f[c], chunks[i]->v, what, n,
                                  v.eps_component(), chunks[i]->ma->eps);
      all_wait();
      if (am_master())
        output_complex_eps_tail((component)c, v, what, n);
    }
  free(n);
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

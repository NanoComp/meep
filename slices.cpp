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

#include "meep.h"
#include "meep_internals.h"

complex<double> fields::optimal_phase_shift(component c) const {
  complex<double> mean = field_mean(c, false, false);
  complex<double> meanr= field_mean(c, true , false);
  complex<double> meani= field_mean(c, false, true );
  if (abs(mean) > abs(meanr) && abs(mean) > abs(meani)) return abs(mean)/mean;
  if (abs(meanr) > abs(meani)) return abs(meanr)/meanr;
  if (abs(meani) > 0.0) return abs(meani)/meani;
  return 1.0;
}

complex<double> fields::field_mean(component c, bool abs_real,
                                   bool abs_imag) const {
  complex<double> themean = 0.0;
  for (int i=0;i<num_chunks;i++)
    themean += chunks[i]->field_mean(c, abs_real, abs_imag);
  return sum_to_all(themean);
}

complex<double> fields_chunk::field_mean(component c, bool abs_real,
                                         bool abs_imag) const {
  complex<double> themean = 0.0;
  if (f[c][0]) for (int i=0;i<v.ntot();i++)
    themean += (v.dV(c,i) & gv).full_volume()*
      (abs_real)?fabs(f[c][0][i]):f[c][0][i];
  if (f[c][1]) for (int i=0;i<v.ntot();i++)
    themean += (v.dV(c,i) & gv).full_volume()*
      complex<double>(0,abs_imag?fabs(f[c][1][i]):f[c][1][i]);
  return themean;
}

static void output_complex_slice(component m, double *f[2],
                                 complex<double> phshift, const volume &v,
                                 const geometric_volume &what, file *out) {
  if (!f[0] || ! f[1]) return; // Field doesn't exist...
  double c = real(phshift), s = imag(phshift);
  for (int i=0;i<v.ntot();i++) {
    if (what.contains(v.loc(m,i))) {
      v.loc(m,i).print(out);
      i_fprintf(out, "\t%lg\n", c*f[0][i]+s*f[1][i]);
    }
  }
}

static void output_slice(component m, const double *f, const volume &v,
                         const geometric_volume &what, file *out) {
  if (!f) return; // Field doesn't exist...
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i))) {
      v.loc(m,i).print(out);
      i_fprintf(out, "\t%.18lg\n", f[i]);
    }
}

static const double default_eps_resolution = 20.0;
static const double default_eps_size = 1000.0;

static void eps_header(double xmin, double ymin, double xmax, double ymax,
                       double fmax, double a, file *out, const char *name) {
  i_fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
  double size = xmax - xmin + ymax - ymin;
  i_fprintf(out, "%%%%BoundingBox: %lg %lg %lg %lg\n",
           xmin*default_eps_size/size, ymin*default_eps_size/size,
           xmax*default_eps_size/size, ymax*default_eps_size/size);
  i_fprintf(out, "gsave\n");
  i_fprintf(out, "%lg %lg scale\n", default_eps_size/size, default_eps_size/size);
  i_fprintf(out, "/Times-Roman findfont 20 scalefont setfont\n");
  i_fprintf(out, "newpath 140 280 moveto (%s) show\n", name);
  i_fprintf(out, "/max %lg def\n", fmax);
  const double dx = 1.0/min(a, default_eps_resolution);
  i_fprintf(out, "/dx %lg def\n", dx);
  i_fprintf(out, "/hdx %lg def\n", dx*0.5);
  i_fprintf(out, "dx 10 div setlinewidth\n");
  i_fprintf(out, "/dotrad %lg def\n", dx*0.25);
  i_fprintf(out, "1 setlinecap\n");
  i_fprintf(out, "/P {\n\
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
  i_fprintf(out, "    %lg %lg rmoveto\n", dx*0.5, dx*0.5);
  i_fprintf(out, "    0 %lg rlineto\n", -dx);
  i_fprintf(out, "    %lg 0 rlineto\n", -dx);
  i_fprintf(out, "    0 dx rlineto\n", dx);
  i_fprintf(out, "    dx 0 rlineto\n", dx);
  i_fprintf(out, "    gsave\n\
    fill\n\
    grestore\n");
  i_fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  i_fprintf(out, "    stroke\n\
} def\n\
/LV {\n\
    0 0 0 setrgbcolor\n\
    moveto\n");
  i_fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  i_fprintf(out, "    0 %lg rmoveto\n", dx*0.5);
  i_fprintf(out, "    0 %lg rlineto\n", -dx);
  i_fprintf(out, "    stroke\n\
} def\n\
/LH {\n");
  i_fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  i_fprintf(out, "    0 0 0 setrgbcolor\n\
    moveto\n");
  i_fprintf(out, "    %lg 0 rmoveto\n", dx*0.5);
  i_fprintf(out, "    %lg 0 rlineto\n", -dx);
  i_fprintf(out, "    stroke\n\
} def\n");
  i_fprintf(out, "    /DV { [0 %lg] 0 setdash LV } def\n", dx/4);
  i_fprintf(out, "    /DH { [0 %lg] 0 setdash LH } def\n", dx/4);
  i_fprintf(out, "    /D { moveto\n\
    %lg setlinewidth\n\
    0 0.8 0 setrgbcolor [0 %lg] 0 setdash\n\
    lineto stroke } def\n", 0.6*dx, 3*dx);
  // B for Boundary...
  i_fprintf(out, "/B {\n\
    0 0 0 setrgbcolor \n\
    newpath dotrad 0 360 arc fill \n\
} def\n");
}

static void eps_1d_header(double xmin, double ymin, double xmax, double ymax,
                          double fmax, double a, file *out, const char *name) {
  i_fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
  const double size = xmax - xmin;
  const double fsize = (5.0 < 0.2*size)?0.2*size:5.0;
  const double dx = 1.0/a;
  i_fprintf(out, "%%%%BoundingBox: %lg %lg %lg %lg\n",
           xmin*default_eps_size/size, -default_eps_size*0.5*fsize/size,
           xmax*default_eps_size/size, default_eps_size*0.5*fsize/size);
  i_fprintf(out, "gsave\n");
  i_fprintf(out, "%lg 0 moveto %lg 0 lineto %lg setlinewidth stroke\n",
           xmin*default_eps_size/size, xmax*default_eps_size/size, dx*0.1); 
  i_fprintf(out, "%lg %lg scale\n", default_eps_size/size, default_eps_size/size);
  i_fprintf(out, "/height %lg def\n", (default_eps_size*0.5*fsize)/size);
  i_fprintf(out, "/Times-Roman findfont 12 scalefont setfont\n");
  i_fprintf(out, "newpath 220 height 0.75 mul moveto (%s) show\n", name);
  i_fprintf(out, "/max %lg def\n", fmax);
  i_fprintf(out, "/fscale %lg def\n", 2.2*fmax/fsize);
  i_fprintf(out, "/dotrad %lg def\n", dx);
  i_fprintf(out, "/dx %lg def\n", dx);
  i_fprintf(out, "/hdx %lg def\n", dx*0.5);
  i_fprintf(out, "dx 10 div setlinewidth\n");
  i_fprintf(out, "1 setlinecap\n");
  i_fprintf(out, "/P {\n\
    fscale div exch pop \n\
    newpath dotrad 0 360 arc fill \n\
} def\n");
  i_fprintf(out, "/LV {\n\
    0.8 0.8 0 setrgbcolor\n\
    pop dup height moveto height neg lineto\n");
  i_fprintf(out, "    %lg setlinewidth\n", dx*0.5);
  i_fprintf(out, "    stroke\n\
} def\n\
/LH {\n");
  i_fprintf(out, "    %lg setlinewidth\n", dx*0.1);
  i_fprintf(out, "    0 0 0 setrgbcolor\n\
    moveto\n");
  i_fprintf(out, "    %lg 0 rmoveto\n", 10*dx*0.5);
  i_fprintf(out, "    %lg 0 rlineto\n", -10*dx);
  i_fprintf(out, "    stroke\n\
} def\n");
  i_fprintf(out, "    /DV { [0 %lg] 0 setdash LV } def\n", dx/4);
  i_fprintf(out, "    /D { moveto\n\
    %lg setlinewidth\n\
    0 0.8 0 setrgbcolor [0 %lg] 0 setdash\n\
    lineto stroke } def\n", 0.6*dx, 3*dx);
}

static void eps_trailer(file *out) {
  i_fprintf(out, "grestore\n");
  i_fprintf(out, "showpage\n");
  i_fprintf(out, "%%%%Trailer\n");
  i_fprintf(out, "%%%%EOF\n");
}

static void eps_dotted(file *out, component m, const double *f, const volume &v,
                       const geometric_volume &what) {
  if (!f) return; // Field doesn't exist...
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i)))
      switch (v.dim) {
      case Dcyl:
        {
          ivec next = v.iloc(m,i)+ivec(2,0);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            i_fprintf(out, "%lg\t%lg\tDH\n", v[next].z(), v[next].r() - 0.5/v.a);
          next = v.iloc(m,i)+ivec(0,2);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            i_fprintf(out, "%lg\t%lg\tDV\n", v[next].z() - 0.5/v.a, v[next].r());
          break;
        }
      case D2:
        {
          ivec next = v.iloc(m,i)+ivec2d(2,0);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            i_fprintf(out, "%lg\t%lg\tDH\n", v[next].x() - 0.5/v.a, v[next].y());
          next = v.iloc(m,i)+ivec2d(0,2);
          if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
              f[i]*f[v.index(m,next)] == 0.0)
            i_fprintf(out, "%lg\t%lg\tDV\n", v[next].x(), v[next].y() - 0.5/v.a);
          break;
        }
      }
}

void fields::outline_chunks(file *out) {
  if (v.dim == D1) return;
  if (my_rank()) return;
  for (int i=0;i<num_chunks;i++) {
    double xlo, xhi, ylo, yhi;
    switch (v.dim) {
    case Dcyl:
      xlo = chunks[i]->v.boundary_location(Low,Z);
      xhi = chunks[i]->v.boundary_location(High,Z);
      ylo = chunks[i]->v.boundary_location(Low,R);
      yhi = chunks[i]->v.boundary_location(High,R);
    break;
    case D2:
      xlo = chunks[i]->v.boundary_location(Low,X);
      ylo = chunks[i]->v.boundary_location(Low,Y);
      xhi = chunks[i]->v.boundary_location(High,X);
      yhi = chunks[i]->v.boundary_location(High,Y);
    break;
    case D3: // FIXME make this smart about plane direction.
      xlo = chunks[i]->v.boundary_location(Low,X);
      ylo = chunks[i]->v.boundary_location(Low,Y);
      xhi = chunks[i]->v.boundary_location(High,X);
      yhi = chunks[i]->v.boundary_location(High,Y);
    }
    i_fprintf(out, "%lg\t%lg\t%lg\t%lg\tD\n", xlo, yhi, xhi, yhi);
    i_fprintf(out, "%lg\t%lg\t%lg\t%lg\tD\n", xhi, yhi, xhi, ylo);
  }
}

static void eps_outline(component m, const double *f,
                        const volume &v, const geometric_volume &what,
                        symmetry S, int symnum, file *out) {
  if (!f) return; // Field doesn't exist...
  for (int i=0;i<v.ntot();i++) {
    const vec here = S.transform(v.loc(m,i),symnum);
    if (what.contains(here))
      switch (v.dim) {
      case Dcyl: {
        ivec next = v.iloc(m,i)+ivec(2,0);
        vec nextrot = v[S.transform(next - ivec(1,0),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          i_fprintf(out, "%lg\t%lg\tLH\n", nextrot.z(), nextrot.r());
        next = v.iloc(m,i)+ivec(0,2);
        nextrot = v[S.transform(next - ivec(0,1),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          i_fprintf(out, "%lg\t%lg\tLV\n", nextrot.z(), nextrot.r());
        break;
      }
      case D2: {
        ivec next = v.iloc(m,i)+ivec2d(0,2);
        vec nextrot = v[(S.transform(next - ivec2d(0,1),symnum))];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            i_fprintf(out, "%lg\t%lg\tLH\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)-ivec2d(0,2);
        nextrot = v[S.transform(next + ivec2d(0,1),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            i_fprintf(out, "%lg\t%lg\tLH\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)+ivec2d(2,0);
        nextrot = v[S.transform(next - ivec2d(1,0),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            i_fprintf(out, "%lg\t%lg\tLV\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)-ivec2d(2,0);
        nextrot = v[S.transform(next + ivec2d(1,0),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            i_fprintf(out, "%lg\t%lg\tLV\n", nextrot.x(), nextrot.y());
        break;
      }
      case D1: {
        const ivec next = v.iloc(m,i)+ivec(2);
        const vec nextrot = v[S.transform(next - ivec(1),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          i_fprintf(out, "%lg\t%lg\tLV\n", nextrot.z(), 0.0);
        break;
      }
      }
  }
}

static void output_complex_eps_body(component m, double *f[2], const volume &v,
                                    symmetry S, int symnum,
                                    const geometric_volume &what, file *out) {
  if (!f[0] || !f[1]) return; // Field doesn't exist...
  const complex<double> ph = S.phase_shift(m, symnum);
  for (int i=0;i<v.ntot();i++) {
    const vec here = S.transform(v.loc(m,i),symnum);
    if (what.contains(here)) {
      double x = 0, y = 0;
      switch (v.dim) {
      case Dcyl: x = here.z(); y = here.r(); break;
      case D1: x = here.z(); break;
      case D2: x = here.x(); y = here.y(); break;
      case D3: x = here.x(); y = here.y(); break; // FIXME use right directions!
      }
      if (f[1]) i_fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y,
                         real(ph)*f[0][i] - imag(ph)*f[1][i]);
      else i_fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y, real(ph)*f[0][i]);
    }
  }
}

void fields_chunk::output_eps_body(component c, const symmetry &S, int sn,
                                   const geometric_volume &what, file *out) {
  if (f[c][0]) {
    if (v.a <= default_eps_resolution) {
      output_complex_eps_body(c, f[c], v, S, sn, what, out);
    } else {
      const complex<double> ph = S.phase_shift(c, sn);
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
          ivec ilocs[8];
          double w[8];
          double fhere = 0.0;
          v.interpolate(c, here, ilocs, w);
          for (int i=0;i<8&&w[i];i++)
            if (v.contains(S.transform(ilocs[i],sn)))
              fhere += w[i]*real(get_field(c,ilocs[i]));
          i_fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y, fhere);
        }
      }
    }
  }
}

static void output_complex_eps_header(component m, double fmax, const volume &v,
                                      const geometric_volume &what, file *out,
                                      const char *name, component om = Hx) {
  double xmin = 1e300, ymin = 1e300, xmax = -1e300, ymax = -1e300;
  direction xdir, ydir;
  switch (v.dim) {
  case Dcyl: xdir = Z; ydir = R; break;
  case D3: xdir = X; ydir = Y; break; // FIXME: check the thin direction of what.
  case D2: xdir = X; ydir = Y; break;
  case D1: xdir = Z; ydir = Z; break; // a tad ugly...
  }
  xmin = what.in_direction_min(xdir);
  xmax = what.in_direction_max(xdir);
  ymin = what.in_direction_min(ydir);
  ymax = what.in_direction_max(ydir);
  if (v.dim == D1) { ymin = -v.inva*0.5; ymax = -ymin; }
  if (ymax == ymin) ymax = ymin + v.inva;
  if (fmax == 0.0) fmax = 0.0001;
  if (v.dim == D1) {
    // Make a 1D line plot!
    eps_1d_header(xmin, ymin, xmax, ymax, fmax, v.a, out, name);
  } else {
    // Make a 2D color plot!
    eps_header(xmin, ymin, xmax, ymax, fmax, v.a, out, name);
  }
}

static void output_eps_header(double fmax, double dx,
                              const double xmin, const double xmax,
                              const double ymin, const double ymax,
                              file *out, const char *name) {
  if (fmax == 0.0) fmax = 0.0001;
  // Make a 2D color plot!
  eps_header(xmin, ymin, xmax, ymax, fmax, 1.0/dx, out, name);
}

static void output_complex_eps_tail(file *out) {
  eps_trailer(out);
}

void mat::output_slices(const char *name) const {
  output_slices(v.surroundings(),name);
}
void mat::output_slices(const geometric_volume &what, const char *name) const {
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = new char[buflen];
  if (!n) abort("Allocation failure!\n");
  snprintf(n, buflen, "%s/%sepsilon.sli", outdir, nname);
  file *out = everyone_open_write(n);
  if (!out) {
    printf("Unable to open file '%s' for slice output.\n", n);
    return;
  }
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      output_slice(v.eps_component(), chunks[i]->eps, chunks[i]->v, what, out);
  everyone_close(out);

  FOR_DIRECTIONS(d)
    FOR_COMPONENTS(c) {
      snprintf(n, buflen, "%s/%ssigma-C-%s-%s.sli", outdir, nname, component_name(c), direction_name(d));
      file *out = everyone_open_write(n);
      if (!out) {
	printf("Unable to open file '%s' for slice output.\n", n);
	return;
      }
      for (int i=0;i<num_chunks;i++)
	if (chunks[i]->is_mine())
	  output_slice(c, chunks[i]->C/*decay*/[d][c]/*[component_direction(c)]*/, chunks[i]->v, what, out);
      everyone_close(out);
    }
  delete[] n;
}

void fields::output_real_imaginary_slices(const char *name) {
  output_real_imaginary_slices(v.surroundings(),name);
}

void fields::output_real_imaginary_slices(const geometric_volume &what,
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
        file *out = everyone_open_write(n);
        if (!out) {
          printf("Unable to open file '%s' for slice output.\n", n);
          return;
        }
        for (int i=0;i<num_chunks;i++)
          if (chunks[i]->is_mine())
            output_slice((component)c, chunks[i]->f[c][cmp],
                         chunks[i]->v, what, out);
        everyone_close(out);
      }
    r_or_i = "-im";
  }

  free(n);
}

void fields::output_slices(const char *name) {
  output_slices(v.surroundings(), name);
}
void fields::output_slices(const geometric_volume &what, const char *name) {
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
  FOR_COMPONENTS(c)
    if (v.has_field(c)) {
      snprintf(n, buflen, "%s/%s%s-%s.sli", outdir, nname,
               component_name(c), time_step_string);
      file *out = everyone_open_write(n);
      if (!out) {
        printf("Unable to open file '%s' for slice output.\n", n);
        return;
      }
      complex<double> phshift = optimal_phase_shift(c);
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          output_complex_slice(c, chunks[i]->f[c], phshift,
                               chunks[i]->v, what, out);
      everyone_close(out);
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
    eps_slices(user_volume.surroundings(), name);
}

void fields::eps_slices(const vec &origin, const vec &xside, const vec &yside,
                        const double dx, const char *name) {
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
  // Define some convenient grid-related values:
  const vec xhat = xside / abs(xside);
  const vec yhat = yside / abs(yside);
  const double xlen = abs(xside);
  const double ylen = abs(yside);
  const double xmin = origin & xhat;
  const double ymin = origin & yhat;
  // Now start outputing pretty pictures.
  FOR_COMPONENTS(c)
    if (v.has_field(c)) {
      snprintf(n, buflen, "%s/%s%s-%s.eps", outdir, nname,
               component_name(c), time_step_string);
      const double fmax = maxfieldmag_to_master(c);
      file *out = everyone_open_write(n);
      if (!out) {
        printf("Unable to open file '%s' for slice output.\n", n);
        return;
      }
      if (am_master())
        output_eps_header(fmax, dx, xmin, xmin + xlen, ymin, ymin + ylen, out, n);
      complex<double> phshift = optimal_phase_shift(c);
      for (double x = xmin; x <= xmin + xlen + dx; x += dx)
        for (double y = ymin; y <= ymin + ylen + dx; y += dx)
          master_fprintf(out, "%lg\t%lg\t%lg\tP\n", x, y,
                         real(phshift*
                              get_field(c, origin + xhat*(x-xmin) + yhat*(y-ymin))));
      for (double x = xmin; x <= xmin + xlen + dx; x += 1.0/v.a)
        for (double y = ymin; y <= ymin + ylen + dx; y += 1.0/v.a) {
          vec loc = origin + xhat*(x-xmin) + yhat*(y-ymin);
          if (has_eps_interface(&loc))
            master_fprintf(out, "%lg\t%lg\tB\n", loc & xhat, loc & yhat);
        }
      //outline_chunks(out);
      //all_wait();
      if (am_master()) output_complex_eps_tail(out);
      everyone_close(out);
    }
  free(n);
  finished_working();
}

bool fields::has_eps_interface(vec *loc) const {
  ivec ilocs[8];
  double w[8];
  double val[8];
  for (int i=0;i<8;i++) val[i] = 0.0;
  v.interpolate(v.eps_component(), *loc, ilocs, w);
  for (int argh=0;argh<8&&w[argh];argh++)
    val[argh] = get_eps(ilocs[argh]);
  double epshere = 0.0, epsmin = 1e100, epsmax = -1e100;
  for (int argh=0;argh<8&&w[argh];argh++) {
    epsmin = min(epsmin, val[argh]);
    epsmax = max(epsmax, val[argh]);
    epshere += w[argh] * val[argh];
  }
  vec grad = zero_vec(v.dim);
  vec center = zero_vec(v.dim);
  double norm = 0.0, mean = 0.0;
  for (int argh=0;argh<8&&w[argh];argh++) {
    center += v[ilocs[argh]];
    mean += val[argh];
    norm++;
  }
  center = center/norm;
  mean *= 1.0/norm;
  for (int argh=0;argh<8&&w[argh];argh++) {
    const vec dist = v[ilocs[argh]] - center;
    grad += dist*(val[argh] - mean)/(dist & dist);
  }
  for (int argh=0;argh<8&&w[argh];argh++)
    if (val[argh] != val[0]){
      *loc = *loc + grad*(0.5*(epsmax - epshere) - mean)/abs(grad & grad);
      return true;
    }
  return false;
}

void fields::eps_slices(const geometric_volume &what, const char *name) {
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
      file *out = everyone_open_write(n);
      if (!out) {
        printf("Unable to open file '%s' for slice output.\n", n);
        return;
      }
      if (am_master())
        output_complex_eps_header((component)c, fmax,
                                  user_volume, what,
                                  out, n, v.eps_component());
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            for (int otherc=0;otherc<10;otherc++)
              if (S.transform((component)otherc,sn) == c)
                chunks[i]->output_eps_body((component)otherc,
                                           S, sn, what, out);
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            for (int otherc=0;otherc<10;otherc++)
              if (S.transform((component)otherc,sn) == c)
                eps_outline(v.eps_component(), chunks[i]->ma->eps,
                            chunks[i]->v, what, S, sn, out);
      all_wait();
      outline_chunks(out);
      all_wait();
      if (am_master()) output_complex_eps_tail(out);
      everyone_close(out);
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

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
#include <stdarg.h>

#include "meep.h"
#include "meep_internals.h"

namespace meep {

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

const int bufsize = 1<<16;

class bufprint {
public:
  bufprint(file *out);
  ~bufprint();
  void flushme();
  void printf(char *fmt, ...);
private:
  file *f;
  char *buf;
  int inbuf;
};

bufprint::bufprint(file *out) {
  f = out;
  buf = new char[bufsize];
  if (!buf) abort("not enough memory for a buffer!");
  inbuf = 0;
}

bufprint::~bufprint() {
  flushme();
  delete[] buf;
}

void bufprint::flushme() {
  if (inbuf) {
    buf[inbuf] = 0;
    i_fprintf(f, "%s", buf);
    inbuf = 0;
  }
}

void bufprint::printf(char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  int written = vsnprintf(buf + inbuf, bufsize - inbuf - 1, fmt, ap);
  if (written <= 0 || written > bufsize - inbuf - 1) {
    flushme();
    written = vsnprintf(buf + inbuf, bufsize - inbuf - 1, fmt, ap);
  }
  inbuf += written;
  if (inbuf > bufsize/2) flushme();
}

static void output_complex_slice(component m, double *f[2],
                                 complex<double> phshift, const volume &v,
                                 const geometric_volume &what, file *out) {
  if (!f[0] || ! f[1]) return; // Field doesn't exist...
  double c = real(phshift), s = imag(phshift);
  bufprint buf(out);
  for (int i=0;i<v.ntot();i++) {
    if (what.contains(v.loc(m,i))) {
      v.loc(m,i).print(out);
      buf.printf("\t%g\n", c*f[0][i]+s*f[1][i]);
    }
  }  
}

static void output_slice(component m, const double *f, const volume &v,
                         const geometric_volume &what, file *out) {
  if (!f) return; // Field doesn't exist...
  bufprint buf(out);
  for (int i=0;i<v.ntot();i++)
    if (what.contains(v.loc(m,i))) {
      v.loc(m,i).print(out);
      buf.printf("\t%.18lg\n", f[i]);
    }
}

static const double default_eps_resolution = 20.0;
static const double default_eps_size = 1000.0;

static void eps_header(double xmin, double ymin, double xmax, double ymax,
                       double fmax, double a, file *out, const char *name) {
  bufprint buf(out);
  buf.printf("%%!PS-Adobe-3.0 EPSF\n");
  double size = xmax - xmin + ymax - ymin;
  buf.printf("%%%%BoundingBox: 0 0 %g %g\n",
           (xmax-xmin)*default_eps_size/size, (ymax-ymin)*default_eps_size/size);
  buf.printf("gsave\n");
  buf.printf("gsave\n");
  buf.printf("/title (%s) def\n", name);
  buf.printf("%g %g scale\n", default_eps_size/size, default_eps_size/size);
  buf.printf("%g %g translate\n", -xmin, -ymin);
  buf.printf("/max %g def\n", fmax);
  const double dx = 1.0/min(a, default_eps_resolution);
  buf.printf("/dx %g def\n", dx);
  buf.printf("/hdx %g def\n", dx*0.5);
  buf.printf("dx 10 div setlinewidth\n");
  buf.printf("/dotrad %g def\n", dx*0.25);
  buf.printf("1 setlinecap\n");
  buf.printf("/P {\n\
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
  buf.printf("    %g %g rmoveto\n", dx*0.5, dx*0.5);
  buf.printf("    0 %g rlineto\n", -dx);
  buf.printf("    %g 0 rlineto\n", -dx);
  buf.printf("    0 dx rlineto\n", dx);
  buf.printf("    dx 0 rlineto\n", dx);
  buf.printf("    gsave\n\
    fill\n\
    grestore\n");
  buf.printf("    %g setlinewidth\n", dx*0.1);
  buf.printf("    stroke\n\
} def\n\
/LV {\n\
    0 0 0 setrgbcolor\n\
    moveto\n");
  buf.printf("    %g setlinewidth\n", dx*0.1);
  buf.printf("    0 %g rmoveto\n", dx*0.5);
  buf.printf("    0 %g rlineto\n", -dx);
  buf.printf("    stroke\n\
} def\n\
/LH {\n");
  buf.printf("    %g setlinewidth\n", dx*0.1);
  buf.printf("    0 0 0 setrgbcolor\n\
    moveto\n");
  buf.printf("    %g 0 rmoveto\n", dx*0.5);
  buf.printf("    %g 0 rlineto\n", -dx);
  buf.printf("    stroke\n\
} def\n");
  buf.printf("    /DV { [0 %g] 0 setdash LV } def\n", dx/4);
  buf.printf("    /DH { [0 %g] 0 setdash LH } def\n", dx/4);
  buf.printf("    /D { moveto\n\
    %g setlinewidth\n\
    0 0.8 0 setrgbcolor [0 %g] 0 setdash\n\
    lineto stroke } def\n", 0.6*dx, 3*dx);
  // B for Boundary...
  buf.printf("/B {\n\
    0 0 0 setrgbcolor \n\
    newpath dotrad 0 360 arc fill \n\
} def\n");
}

static void eps_1d_header(double xmin, double ymin, double xmax, double ymax,
                          double fmax, double a, file *out, const char *name) {
  bufprint buf(out);
  buf.printf("%%!PS-Adobe-3.0 EPSF\n");
  const double size = xmax - xmin;
  const double fsize = (5.0 < 0.2*size)?0.2*size:5.0;
  const double dx = 1.0/a;
  buf.printf("%%%%BoundingBox: 0 0 %g %g\n",
           (xmax-xmin)*default_eps_size/size, default_eps_size*fsize/size);
  buf.printf("gsave\n");
  buf.printf("gsave\n");
  buf.printf("/title (%s) def\n", name);
  buf.printf("%g %g scale\n", default_eps_size/size, default_eps_size/size);
  buf.printf("%g %g translate\n", -xmin, 0.5*fsize);
  buf.printf("%g 0 moveto %g 0 lineto %g setlinewidth stroke\n",
            xmin, xmax, dx*0.1); 
  buf.printf("/height %g def\n", (default_eps_size*0.5*fsize)/size);
  buf.printf("/max %g def\n", fmax);
  buf.printf("/fscale %g def\n", 2.2*fmax/fsize);
  buf.printf("/dotrad %g def\n", dx);
  buf.printf("/dx %g def\n", dx);
  buf.printf("/hdx %g def\n", dx*0.5);
  buf.printf("dx 10 div setlinewidth\n");
  buf.printf("1 setlinecap\n");
  buf.printf("/P {\n\
    fscale div exch pop \n\
    newpath dotrad 0 360 arc fill \n\
} def\n");
  buf.printf("/LV {\n\
    0.8 0.8 0 setrgbcolor\n\
    pop dup height moveto height neg lineto\n");
  buf.printf("    %g setlinewidth\n", dx*0.5);
  buf.printf("    stroke\n\
} def\n\
/LH {\n");
  buf.printf("    %g setlinewidth\n", dx*0.1);
  buf.printf("    0 0 0 setrgbcolor\n\
    moveto\n");
  buf.printf("    %g 0 rmoveto\n", 10*dx*0.5);
  buf.printf("    %g 0 rlineto\n", -10*dx);
  buf.printf("    stroke\n\
} def\n");
  buf.printf("    /DV { [0 %g] 0 setdash LV } def\n", dx/4);
  buf.printf("    /D { moveto\n\
    %g setlinewidth\n\
    0 0.8 0 setrgbcolor [0 %g] 0 setdash\n\
    lineto stroke } def\n", 0.6*dx, 3*dx);
}

static void eps_trailer(file *out) {
  bufprint buf(out);
  buf.printf("grestore\n");
  buf.printf(" 0.25 0 0.5 setrgbcolor\n");
  buf.printf("/Times-Roman findfont 16 scalefont setfont\n");
  buf.printf("newpath 5 5 moveto title show\n");
  buf.printf("grestore\n");
  buf.printf("showpage\n");
  buf.printf("%%%%Trailer\n");
  buf.printf("%%%%EOF\n");
}

// static void eps_dotted(file *out, component m, const double *f, const volume &v,
//                        const geometric_volume &what) {
//   bufprint buf(out);
//   if (!f) return; // Field doesn't exist...
//   for (int i=0;i<v.ntot();i++)
//     if (what.contains(v.loc(m,i)))
//       switch (v.dim) {
//       case Dcyl:
//         {
//           ivec next = v.iloc(m,i)+ivec(2,0);
//           if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
//               f[i]*f[v.index(m,next)] == 0.0)
//             buf.printf("%g\t%g\tDH\n", v[next].z(), v[next].r() - 0.5/v.a);
//           next = v.iloc(m,i)+ivec(0,2);
//           if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
//               f[i]*f[v.index(m,next)] == 0.0)
//             buf.printf("%g\t%g\tDV\n", v[next].z() - 0.5/v.a, v[next].r());
//           break;
//         }
//       case D2:
//         {
//           ivec next = v.iloc(m,i)+ivec2d(2,0);
//           if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
//               f[i]*f[v.index(m,next)] == 0.0)
//             buf.printf("%g\t%g\tDH\n", v[next].x() - 0.5/v.a, v[next].y());
//           next = v.iloc(m,i)+ivec2d(0,2);
//           if (v.contains(next) && f[i] + f[v.index(m,next)] != 0.0 &&
//               f[i]*f[v.index(m,next)] == 0.0)
//             buf.printf("%g\t%g\tDV\n", v[next].x(), v[next].y() - 0.5/v.a);
//           break;
//         }
//       case D1: case D3: break;
//       }
// }

void fields::outline_chunks(file *out) {
  bufprint buf(out);
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
    case D1: abort("Error in outline chunks 1D.\n"); break;
    }
    buf.printf("%g\t%g\t%g\t%g\tD\n", xlo, yhi, xhi, yhi);
    buf.printf("%g\t%g\t%g\t%g\tD\n", xhi, yhi, xhi, ylo);
  }
}

static void eps_outline(component m, const double *f,
                        const volume &v, const geometric_volume &what,
                        symmetry S, int symnum, file *out) {
  bufprint buf(out);
  if (!f) return; // Field doesn't exist...
  for (int i=0;i<v.ntot();i++) {
    const vec here = S.transform(v.loc(m,i),symnum);
    if (what.contains(here))
      switch (v.dim) {
      case Dcyl: {
        ivec next = v.iloc(m,i)+ivec(2,0);
        vec nextrot = v[S.transform(next - ivec(1,0),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          buf.printf("%g\t%g\tLH\n", nextrot.z(), nextrot.r());
        next = v.iloc(m,i)+ivec(0,2);
        nextrot = v[S.transform(next - ivec(0,1),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          buf.printf("%g\t%g\tLV\n", nextrot.z(), nextrot.r());
        break;
      }
      case D2: {
        ivec next = v.iloc(m,i)+ivec2d(0,2);
        vec nextrot = v[(S.transform(next - ivec2d(0,1),symnum))];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            buf.printf("%g\t%g\tLH\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)-ivec2d(0,2);
        nextrot = v[S.transform(next + ivec2d(0,1),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            buf.printf("%g\t%g\tLH\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)+ivec2d(2,0);
        nextrot = v[S.transform(next - ivec2d(1,0),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            buf.printf("%g\t%g\tLV\n", nextrot.x(), nextrot.y());
        next = v.iloc(m,i)-ivec2d(2,0);
        nextrot = v[S.transform(next + ivec2d(1,0),symnum)];
        if (v.owns(next))
          if (f[i] != f[v.index(m,next)])
            buf.printf("%g\t%g\tLV\n", nextrot.x(), nextrot.y());
        break;
      }
      case D1: {
        const ivec next = v.iloc(m,i)+ivec(2);
        const vec nextrot = v[S.transform(next - ivec(1),symnum)];
        if (v.contains(next) && f[i] != f[v.index(m,next)])
          buf.printf("%g\t%g\tLV\n", nextrot.z(), 0.0);
        break;
      }
      case D3: abort("Error in eps_outline.\n"); break;
      }
  }
}

static void output_complex_eps_body(component m, double *f[2], const volume &v,
                                    symmetry S, int symnum,
                                    const geometric_volume &what, file *out) {
  bufprint buf(out);
  if (!f[0]) return; // Field doesn't exist...
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
      if (f[1]) buf.printf("%g\t%g\t%g\tP\n", x, y,
                         real(ph)*f[0][i] - imag(ph)*f[1][i]);
      else buf.printf("%g\t%g\t%g\tP\n", x, y, real(ph)*f[0][i]);
    }
  }
}

void fields_chunk::output_eps_body(component c, const symmetry &S, int sn,
                                   const geometric_volume &what, file *out,
                                   complex<double> phshift) {
  bufprint buf(out);
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
          case D3: abort("Don't support 3D o_e_b\n"); break;
          }
          ivec ilocs[8];
          double w[8];
          double fhere = 0.0;
          v.interpolate(c, here, ilocs, w);
          for (int i=0;i<8&&w[i];i++)
            if (v.contains(S.transform(ilocs[i],sn)))
              fhere += w[i]*real(phshift*get_field(c,ilocs[i]));
          if (fhere != 0.0) // save space by leaving out blanks.
            buf.printf("%g\t%g\t%g\tP\n", x, y, fhere);
        }
      }
    }
  }
}

void fields_chunk::output_eps_body(const polarizability_identifier &p,
                                   component c, const symmetry &S, int sn,
                                   const geometric_volume &what, file *out,
                                   complex<double> phshift) {
  bufprint buf(out);
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
          case D3: abort("Don't support 3D o_e_b\n"); break;
          }
          ivec ilocs[8];
          double w[8];
          double fhere = 0.0;
          v.interpolate(c, here, ilocs, w);
          for (int i=0;i<8&&w[i];i++)
            if (v.contains(S.transform(ilocs[i],sn)))
              fhere += w[i]*real(phshift*get_polarization_field(p,c,ilocs[i]));
          if (fhere != 0.0) // save space by leaving out blanks.
            buf.printf("%g\t%g\t%g\tP\n", x, y, fhere);
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

void structure::output_slices(const char *name) const {
  output_slices(v.surroundings(),name);
}
void structure::output_slices(const geometric_volume &what, const char *name) const {
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
  if (!n) abort("Allocation failure!\n");
  char *r_or_i = "-re";
  DOCMP {
    FOR_COMPONENTS(c)
      if (v.has_field(c)) {
        snprintf(n, buflen, "%s/%s%s%s-%09.2f.sli",
                 outdir, nname, component_name(c),
                 r_or_i, time());
        file *out = everyone_open_write(n);
        if (!out) {
          printf("Unable to open file '%s' for slice output.\n", n);
          return;
        }
        for (int i=0;i<num_chunks;i++)
          if (chunks[i]->is_mine())
            output_slice(c, chunks[i]->f[c][cmp],
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
  //if (new_s) {
  //  snprintf(n, buflen, "%s/%sepsilon-%s.sli", outdir, nname, time_step_string);
  //  output_slice(v.eps_component(), s->eps, v, what, n);
  //}
  free(n);
  finished_working();
}

void fields::eps_energy_slice(const char *name) {
  if (v.dim == D3) abort("FIXME need to support 3D energy slices...\n");
  geometric_volume what = user_volume.surroundings();
  if (v.dim == Dcyl) what.set_direction_min(R,-what.in_direction_max(R));
  eps_energy_slice(what,name);
}

void fields::eps_energy_slice(const geometric_volume &what, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = new char[buflen];
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());
  
  snprintf(n, buflen, "%s/%senergy-%s.eps", outdir, nname, time_step_string);
  file *out = everyone_open_write(n);
  if (!out) abort("Unable to open file '%s' for slice output.\n", n);
  const double fmax = max(maxpolenergy_to_master(), -minpolenergy_to_master());
  if (am_master())
    output_complex_eps_header(v.eps_component(), fmax, user_volume,
                              what, out, n, v.eps_component());
  if (v.dim != D1) abort("Still only works in 1D.  :( \n");
  {
    bufprint buf(out);
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->is_mine())
        for (int sn=0;sn<S.multiplicity();sn++)
          for (int n=0;n<chunks[i]->v.ntot();n++) {
            const ivec here = chunks[i]->v.iloc(v.eps_component(), n);
            buf.printf("%g\t0\t%g\tP\n", chunks[i]->v[S.transform(here,sn)].z(),
                       chunks[i]->get_polarization_energy(here));
          }
  }
  all_wait();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      for (int sn=0;sn<S.multiplicity();sn++)
        for (int otherc=0;otherc<10;otherc++)
          if (S.transform((component)otherc,sn) == c)
            eps_outline(v.eps_component(), chunks[i]->s->eps,
                        chunks[i]->v, what, S, sn, out);
  all_wait();
  if (am_master()) output_complex_eps_tail(out);
  everyone_close(out);
  delete[] n;
  finished_working();
}

void fields::eps_energy_slice(const polarizability_identifier &p, const char *name) {
  if (v.dim == D3) abort("FIXME need to support 3D energy slices...\n");
  geometric_volume what = user_volume.surroundings();
  if (v.dim == Dcyl) what.set_direction_min(R,-what.in_direction_max(R));
  eps_energy_slice(p, what,name);
}

void fields::eps_energy_slice(const polarizability_identifier &p,
                              const geometric_volume &what, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = new char[buflen];
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());
  
  snprintf(n, buflen, "%s/%senergy-%s.eps", outdir, nname, time_step_string);
  file *out = everyone_open_write(n);
  if (!out) abort("Unable to open file '%s' for slice output.\n", n);
  const double fmax = max(maxpolenergy_to_master(), -minpolenergy_to_master());
  if (am_master())
    output_complex_eps_header(v.eps_component(), fmax, user_volume,
                              what, out, n, v.eps_component());
  if (v.dim != D1) abort("Still only works in 1D.  :( \n");
  {
    bufprint buf(out);
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->is_mine())
        for (int sn=0;sn<S.multiplicity();sn++)
          for (int n=0;n<chunks[i]->v.ntot();n++) {
            const ivec here = chunks[i]->v.iloc(v.eps_component(), n);
            buf.printf("%g\t0\t%g\tP\n", chunks[i]->v[S.transform(here,sn)].z(),
                       chunks[i]->get_polarization_energy(p, here));
          }
  }
  all_wait();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      for (int sn=0;sn<S.multiplicity();sn++)
        for (int otherc=0;otherc<10;otherc++)
          if (S.transform((component)otherc,sn) == c)
            eps_outline(v.eps_component(), chunks[i]->s->eps,
                        chunks[i]->v, what, S, sn, out);
  all_wait();
  if (am_master()) output_complex_eps_tail(out);
  everyone_close(out);
  delete[] n;
  finished_working();
}

void fields::eps_polarization_slice(const polarizability_identifier &p, const char *name) {
  if (v.dim == D3) abort("FIXME need to support 3D polarization slices...\n");
  geometric_volume what = user_volume.surroundings();
  if (v.dim == Dcyl) what.set_direction_min(R,-what.in_direction_max(R));
  eps_polarization_slice(p, what,name);
}

void fields::eps_polarization_slice(const polarizability_identifier &p,
                                    const geometric_volume &what, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = new char[buflen];
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());

  FOR_ELECTRIC_COMPONENTS(c)
    if (v.has_field(c)) {
      snprintf(n, buflen, "%s/%sP-%s-%s.eps", outdir, nname,
               component_name(c), time_step_string);
      const double fmax = maxfieldmag_to_master(c);
      file *out = everyone_open_write(n);
      if (!out) {
        printf("Unable to open file '%s' for slice output.\n", n);
        return;
      }
      if (am_master())
        output_complex_eps_header(c, fmax, user_volume, what,
                                  out, n, v.eps_component());
      complex<double> phshift = optimal_phase_shift(c);
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            FOR_COMPONENTS(otherc)
              if (S.transform(otherc,sn) == c)
                chunks[i]->output_eps_body(p, otherc,
                                           S, sn, what, out, phshift);
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            FOR_COMPONENTS(otherc)
              if (S.transform(otherc,sn) == c)
                eps_outline(v.eps_component(), chunks[i]->s->eps,
                            chunks[i]->v, what, S, sn, out);
      all_wait();
      outline_chunks(out);
      all_wait();
      if (am_master()) output_complex_eps_tail(out);
      everyone_close(out);
    }
  delete[] n;
  finished_working();
}

void fields::eps_envelope(const char *name) {
  if (v.dim != D1) abort("Envelope only supported in 1D so far.\n");
  if (v.dim == D3) abort("Specify directions for EPS slices in 3D\n");
  geometric_volume what = user_volume.surroundings();
  if (v.dim == Dcyl) what.set_direction_min(R,-what.in_direction_max(R));
  eps_envelope(what,name);
}

void fields::eps_envelope(const geometric_volume &what, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());
  FOR_COMPONENTS(c)
    if (v.has_field(c)) {
      snprintf(n, buflen, "%s/%s%s-%s.eps", outdir, nname,
               component_name((component)c), time_step_string);
      const double fmax = maxfieldmag_to_master(c);
      file *out = everyone_open_write(n);
      if (!out) {
        printf("Unable to open file '%s' for slice output.\n", n);
        return;
      }
      if (am_master())
        output_complex_eps_header(c, fmax,
                                  user_volume, what,
                                  out, n, v.eps_component());
      for (double z = 0.0 + inva; z < user_volume.ntot(); z += inva)
        if (what.contains(vec(z))) {
          const double fhere = real(get_field(c, vec(z)));
          if (fhere > 0.0 && fhere > real(get_field(c, vec(z)-inva)) &&
              fhere > real(get_field(c, vec(z)+inva)))
            master_fprintf(out, "%g\t0\t%g\tP\n", z, fhere);
        }
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            FOR_COMPONENTS(otherc)
              if (S.transform(otherc,sn) == c)
                eps_outline(v.eps_component(), chunks[i]->s->eps,
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

void fields::eps_slices(const char *name) {
  if (v.dim == D3) abort("Specify directions for EPS slices in 3D\n");
  geometric_volume what = user_volume.surroundings();
  if (v.dim == Dcyl) what.set_direction_min(R,-what.in_direction_max(R));
  eps_slices(what,name);
}

void fields::eps_slices(const vec &origin, const vec &xside, const vec &yside,
                        const double dx, const char *name) {
  am_now_working_on(Slicing);
  const int buflen = 1024;
  char nname[buflen];
  if (*name) snprintf(nname, buflen, "%s-", name);
  else *nname = 0; // No additional name!
  char *n = (char *)malloc(buflen);
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
          master_fprintf(out, "%g\t%g\t%g\tP\n", x, y,
                         real(phshift*
                              get_field(c, origin + xhat*(x-xmin) + yhat*(y-ymin))));
      for (double x = xmin; x <= xmin + xlen + dx; x += 1.0/v.a)
        for (double y = ymin; y <= ymin + ylen + dx; y += 1.0/v.a) {
          vec loc = origin + xhat*(x-xmin) + yhat*(y-ymin);
          if (has_eps_interface(&loc))
            master_fprintf(out, "%g\t%g\tB\n", loc & xhat, loc & yhat);
        }
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
  if (!n) abort("Allocation failure!\n");
  char time_step_string[buflen];
  snprintf(time_step_string, buflen, "%09.2f", time());
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
        output_complex_eps_header(c, fmax, user_volume, what,
                                  out, n, v.eps_component());
      complex<double> phshift = optimal_phase_shift(c);
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            FOR_COMPONENTS(otherc)
              if (S.transform(otherc,sn) == c)
                chunks[i]->output_eps_body(otherc,
                                           S, sn, what, out, phshift);
      all_wait();
      for (int i=0;i<num_chunks;i++)
        if (chunks[i]->is_mine())
          for (int sn=0;sn<S.multiplicity();sn++)
            FOR_COMPONENTS(otherc)
              if (S.transform(otherc,sn) == c)
                eps_outline(v.eps_component(), chunks[i]->s->eps,
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

double fields::minpolenergy_to_master() const {
  double themin = 1e300;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      themin = min(themin,chunks[i]->minpolenergy());
  return -max_to_master(-themin);
}

double fields::maxpolenergy_to_master() const {
  double themax = -1e300;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      themax = max(themax,chunks[i]->maxpolenergy());
  return max_to_master(themax);
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

double fields_chunk::minpolenergy() const {
  double themin = 1e300;
  const component c = v.eps_component();
  for (int i=0;i<v.ntot();i++)
    themin = min(themin,my_polarization_energy(v.iloc(c, i)));
  return themin;
}

double fields_chunk::maxpolenergy() const {
  double themax = -1e300;
  const component c = v.eps_component();
  for (int i=0;i<v.ntot();i++)
    themax = max(themax,my_polarization_energy(v.iloc(c, i)));
  return themax;
}

} // namespace meep

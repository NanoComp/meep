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

/* Below are some routines to output to a grace file. */

const char grace_header[] = "# Grace project file
#
@page size 792, 612
@page scroll 5%
@page inout 5%
@default symbol size 0.330000\n
";

class grace_point {
public:
  int n;
  double x, y, dy;
  grace_point *next;
};

grace::grace(const char *fname, const char *dirname) {
  fn = fname;
  dn = dirname;
  char buf[300];
  snprintf(buf,300,"%s/%s", dirname, fname);
  f = fopen(buf, "w");
  if (!f) {
    printf("Unable to open file %f\n", buf);
    exit(1);
  }
  set_num = 0;
  sn = 0;
  pts = NULL;
  fprintf(f,"%s", grace_header);
}

grace::~grace() {
  flush_pts();
  fclose(f);
  char gracecmd[500];
  snprintf(gracecmd, 500, "gracebat -hdevice EPS -printfile %s/%s.eps -hardcopy %s/%s",
           dn, fn, dn, fn);
  system(gracecmd);
}

void grace::new_set(grace_type pt) {
  flush_pts();
  set_num++;
  sn++;
  fprintf(f, "@    s%d line color %d\n", sn, set_num);
  fprintf(f, "@    s%d symbol color %d\n", sn, set_num);
  fprintf(f, "@    s%d errorbar color %d\n", sn, set_num);
  fprintf(f, "@    target G0.S%d\n", sn);
  if (pt == ERROR_BARS) fprintf(f, "@    type xydy\n");
  else fprintf(f, "@    type xy\n");
}

void grace::set_legend(const char *l) {
  fprintf(f, "@    s%d legend  \"%s\"\n", sn, l);
}

void grace::new_curve() {
  sn++;
  fprintf(f, "@    s%d line color %d\n", sn, set_num);
  fprintf(f, "@    s%d symbol color %d\n", sn, set_num);
  fprintf(f, "@    s%d errorbar color %d\n", sn, set_num);
  if (!set_num) new_set();
  fprintf(f, "\n");
}

void grace::output_point(double x, double y, double dy) {
  if (dy >= 0) fprintf(f, "%lg\t%lg\t%lg\n", x, y, dy);
  else fprintf(f, "%lg\t%lg\n", x, y);
  fflush(f);
}

void grace::output_out_of_order(int n, double x, double y, double dy) {
  grace_point *gp = new grace_point;
  gp->n = n;
  gp->x = x;
  gp->y = y;
  gp->dy = dy;
  gp->next = pts;
  pts = gp;
}

void grace::flush_pts() {
  int first_time = 1;
  while (pts) {
    grace_point *p = pts;
    int num_seen = 0;
    while (p) {
      if (p->n <= 0) num_seen++;
      p = p->next;
    }
    if (num_seen && !first_time) new_curve();
    first_time = 0;
    p = pts;
    grace_point **last = &pts;
    while (p) {
      if (p->n <= 0) {
        *last = p->next;
        output_point(p->x,p->y,p->dy);
        delete p;
        p = *last;
      } else {
        p->n -= 1;
        last = &p->next;
        p = p->next;
      }
    }
  }
}

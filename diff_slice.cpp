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
#include <stdarg.h>
#include <dirent.h>
#include <string.h>

void usage() {
  printf("Usage:  diff_slice <name1> <name2> <small float>\n");
  exit(0);
}

int verbose = 1;

void ack(const char *format, ...);
void hey(const char *format, ...);
void read_slice(const char *fn, double *&a, int &nr, int &nz);
int compare_files(const char *f1, const char *f2, double small);

int main(int argc, const char *argv[]) {
  char *n1 = new char[1024], *n2 = new char[1024];
  if (argc != 4) usage();
  const char *name1 = argv[1];
  const char *name2 = argv[2];
  char time[7] = "000000";
  const double little = atof(argv[3]);
  printf("Comparing %s with %s, small number being %lg.\n",
         name1, name2, little);
  // First compare epsilon...
  sprintf(n1, "%s-epsilon.sli", name1);
  sprintf(n2, "%s-epsilon.sli", name2);
  if (compare_files(n1, n2, little)) ack("Epsilon differs!");
  // Then compare sigma!
  sprintf(n1, "%s-sigma.sli", name1);
  sprintf(n2, "%s-sigma.sli", name2);
  if (compare_files(n1, n2, little)) ack("Sigma differs!");

  char *lookfor = new char[strlen(name1)+1+50];
  strcpy(lookfor,name1);
  strcat(lookfor,"-hz-re-");

  struct dirent **namelist;
  int nfiles;
  nfiles = scandir(".", &namelist, 0, alphasort);
  if (nfiles <= 0) ack("Problem reading directory.");
  else {
    for (int n=0;n<nfiles;n++) {
      char *name = namelist[n]->d_name;
      if (strcmp(name+strlen(name)-4,".sli") == 0 &&
          strncmp(name, lookfor, strlen(lookfor)) == 0) {
        strncpy(time, name + strlen(name)-10, 6);
        int err = 0;
        hey("Checking time %s", time);
        // Check Hr:
        sprintf(n1, "%s-hr-re-%s.sli", name1, time);
        sprintf(n2, "%s-hr-re-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        sprintf(n1, "%s-hr-im-%s.sli", name1, time);
        sprintf(n2, "%s-hr-im-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        // Check Hp:
        sprintf(n1, "%s-hp-re-%s.sli", name1, time);
        sprintf(n2, "%s-hp-re-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        sprintf(n1, "%s-hp-im-%s.sli", name1, time);
        sprintf(n2, "%s-hp-im-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        // Check Hz:
        sprintf(n1, "%s-hz-re-%s.sli", name1, time);
        sprintf(n2, "%s-hz-re-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        sprintf(n1, "%s-hz-im-%s.sli", name1, time);
        sprintf(n2, "%s-hz-im-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        if (err) ack("We found problems in %d H components at time %s.", err, time);        
        // Check Er:
        sprintf(n1, "%s-er-re-%s.sli", name1, time);
        sprintf(n2, "%s-er-re-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        sprintf(n1, "%s-er-im-%s.sli", name1, time);
        sprintf(n2, "%s-er-im-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        // Check Ep:
        sprintf(n1, "%s-ep-re-%s.sli", name1, time);
        sprintf(n2, "%s-ep-re-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        sprintf(n1, "%s-ep-im-%s.sli", name1, time);
        sprintf(n2, "%s-ep-im-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        // Check Ez:
        sprintf(n1, "%s-ez-re-%s.sli", name1, time);
        sprintf(n2, "%s-ez-re-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        sprintf(n1, "%s-ez-im-%s.sli", name1, time);
        sprintf(n2, "%s-ez-im-%s.sli", name2, time);
        err += compare_files(n1, n2, little);
        if (err) ack("We found problems in %d E components at time %s.", err, time);
      }
      free(namelist[n]);
    }
    free(namelist);
  }
}

int compare_files(const char *f1, const char *f2, double small) {
  double *s1, *s2;
  int nz=0,nr=0, err=0;
  
  read_slice(f1,s1,nr,nz);
  read_slice(f2,s2,nr,nz);

  for (int r=0;r<nr;r++) {
    for (int z=0;z<nz;z++) {
      if (fabs(s1[z+r*nz] - s2[z+r*nz]) > small) {
        printf("Err: %s and %s at %4d,%4d: %lg\n",
               f1, f2, r, z, s1[z+r*nz] - s2[z+r*nz]);
        err = 1;
      }
    }
  }
  delete[] s1;
  delete[] s2;
  return err;
}

void read_slice(const char *fn, double *&a, int &nr, int &nz) {
  FILE *f = fopen(fn, "r");
  if (!f) ack("Cannot open file %s", fn);
  int z, r, oldnz = nz, oldnr = nr;
  nz = 0;
  nr = 0;
  double v;
  while (fscanf(f," %d %d %lg", &z, &r, &v) == 3) {
    if (z >= nz) nz = z+1;
    if (r >= nr) nr = r+1;
  }
  fclose(f);
  if ((oldnz && oldnr) && (nz != oldnz || nr != oldnr))
    ack("You changed your dimensions!");
  //hey("Got %d by %d data!", nr, nz);
  a = new double[nz*nr];
  if (!a) ack("Error allocating!");
  f = fopen(fn, "r");
  if (!f) ack("Cannot open file %s", fn);
  while (fscanf(f," %d %d %lg", &z, &r, &v) == 3) {
    a[z+r*nz] = v;
  }
  fclose(f);
}

void hey(const char *ackf, ...) {
  if (verbose) {
    va_list argptr;
    va_start(argptr, ackf);
    
    vprintf(ackf, argptr);
    printf("\n");
  }
}

void ack(const char *ackf, ...) {
  va_list argptr;
  va_start(argptr, ackf);

  printf("Aaack. ");
  vprintf(ackf, argptr);
  printf("\n");
  exit(1);
}

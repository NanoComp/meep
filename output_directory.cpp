/* Copyright (C) 2003 David Roundy  
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
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>

#include "dactyl.h"

void mat::set_output_directory(const char *name) {
  outdir = name;
}

void fields::set_output_directory(const char *name) {
  outdir = name;
}

static int is_same_file(const char *a, const char *b) {
  FILE *fa = fopen(a,"r");
  FILE *fb = fopen(b,"r");
  if (!fa || !fb) return 0;
  int ca;
  do {
    ca = getc(fa);
    if (ca != getc(fb)) return 0;
  } while (ca != EOF);
  return 1;
}

static void cp(const char *a, const char *b) {
  FILE *fa = fopen(a,"r");
  FILE *fb = fopen(b,"w");
  if (!fa || !fb) return;
  int ca;
  while (1) {
    ca = getc(fa);
    if (ca == EOF) break;
    putc(ca,fb);
  }
  fclose(fa);
  fclose(fb);
}

static int is_ok_dir(const char *dirname, const char *sourcename, const char *basename) {
  const int buflen = 300;

  DIR *dir;
  if ((dir = opendir(dirname)) != NULL)
    closedir(dir);
  else {
    mkdir(dirname, 00777);
    return 1;
  }

  char drsrcn[buflen];
  snprintf(drsrcn, buflen, "%s/%s.cpp", dirname, basename);
  if (is_same_file(drsrcn, sourcename)) return 1;
  
  FILE *f;
  if ((f = fopen(drsrcn, "r")) == NULL) return 1;
  fclose(f);
  return 0;
}

const char *make_output_directory(const char *exename) {
  const int buflen = 300;
  char basename[buflen];
  const char *bnp = exename; // basename holds the actual name of the
                                  // executable (dirs removed).
  const char *t;
  for (t=exename;*t;t++) {
    if (*t == '/') bnp = t+1;
  }
  snprintf(basename, buflen, "%s", bnp);
  if (strcmp(basename + strlen(basename) - 4, ".dac") == 0) {
    basename[strlen(basename) - 4] = (char)0;
  }

  char sourcename[buflen]; // Holds the "example.cpp" filename.
  snprintf(sourcename, buflen, "%s.cpp", exename);

  char outdirname[buflen];
  snprintf(outdirname, buflen, "%s-out", basename);
  if (!is_ok_dir(outdirname, sourcename, basename)) {
    for (int i=1;i<100;i++) {
      snprintf(outdirname, buflen, "%s-out-%d", basename, i);
      if (is_ok_dir(outdirname, sourcename, basename)) break;
    }
  }
  
  char outsrcname[buflen];
  snprintf(outsrcname, buflen, "%s/%s.cpp", outdirname, basename);
  cp(sourcename, outsrcname);

  char *dirname = new char[strlen(outdirname)+1];
  snprintf(dirname, strlen(outdirname)+1, "%s", outdirname);
  return dirname;
}

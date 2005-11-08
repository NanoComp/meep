/* Copyright (C) 2005 Massachusetts Institute of Technology  
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
#include <unistd.h>

#include "meep.hpp"

namespace meep {

const char symlink_name[] = "latest_output";

void structure::set_output_directory(const char *name) {
  char buf[300];
  outdir = name;
  if (!quiet) master_printf("Using output directory %s/\n", name);
  if (readlink(symlink_name, buf, 300) > 0) {
    // Link already exists.
    unlink(symlink_name);
  }
  symlink(name, symlink_name);
  outdir = name;
}

void fields::set_output_directory(const char *name) {
  delete[] outdir;
  outdir = new char[strlen(name) + 1]; strcpy(outdir, name);
  for (int i=0;i<num_chunks;i++) chunks[i]->set_output_directory(outdir);
}

void fields_chunk::set_output_directory(const char *name) {
  outdir = name;
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

static bool is_ok_dir(const char *dirname) {
  DIR *dir;
  bool direxists;
  if (am_master()) {
    direxists = (dir = opendir(dirname)) != NULL;
    if (direxists) closedir(dir);
    else mkdir(dirname, 00777);
  }
  direxists = broadcast(0, direxists);
  return !direxists;
}

file *create_output_file(const char *dirname, const char *fname) {
  const int buflen = 300;
  char n[buflen];
  snprintf(n, buflen, "%s/%s", dirname, fname);
  file *o = everyone_open_write(n);
  if (!o) abort("Unable to create file %s!\n", n);
  return o;
}

const char *make_output_directory(const char *exename, const char *jobname) {
  const int buflen = 300;
  char basename[buflen];
  const char * const evil_suffs[] = { ".dac", ".cpp", ".cc", ".cxx", ".C" };
  char stripped_name[buflen];
  const char *bnp = exename; // stripped_name holds the actual name of the executable (dirs removed).
  const char *t;
  for (t=exename;*t;t++) {
    if (*t == '/') bnp = t+1;
  }
  
  snprintf(stripped_name, buflen, "%s", bnp);
  for (int i = 0; i < (int)(sizeof(evil_suffs) / sizeof(evil_suffs[0])); ++i) {
    int sufflen = strlen(evil_suffs[i]);
    if (strcmp(stripped_name + strlen(stripped_name) - sufflen, 
	       evil_suffs[i]) == 0
	&& strlen(stripped_name) > size_t(sufflen)) {
      stripped_name[strlen(stripped_name) - sufflen] = (char)0;
      break;
    }
  }

  char sourcename[buflen]; // Holds the "example.cpp" filename.
  snprintf(sourcename, buflen, "%s.cpp", stripped_name);

  if (jobname != NULL) {
    snprintf(basename, buflen, "%s", jobname);
  }
  else {
    snprintf(basename, buflen, "%s", stripped_name);
  }

  static char outdirname[buflen];
  snprintf(outdirname, buflen, "%s-out", basename);
  {
    int i = 0;
    while (!is_ok_dir(outdirname)) {
      if (!quiet)
	master_printf("Output directory %s already exists!\n", outdirname);
      snprintf(outdirname, buflen, "%s-out-%d", basename, i++);
    }
  }
  char outsrcname[buflen];
  snprintf(outsrcname, buflen, "%s/%s", outdirname, sourcename);
  cp(sourcename, outsrcname);

  return outdirname;
}

void trash_output_directory(const char *dirname) {
  if (am_master()) mkdir(dirname, 00777);
}

} // namespace meep

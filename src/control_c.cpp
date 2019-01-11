/* Copyright (C) 2005-2019 Massachusetts Institute of Technology
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
#include <signal.h>

#include "meep.hpp"

using namespace std;

namespace meep {

int interrupt = 0;
static int kill_time = 2;

static void handle_control_c(int i) {
  (void)i; // unused: should equal SIGINT
  interrupt++;
  if (interrupt >= kill_time) {
    abort("interrupted");
  } else if (interrupt + 1 == kill_time) {
    printf("Be patient... hit ctrl-C one more time to kill me.\n");
  } else {
    printf("Be patient... hit ctrl-C %d more times to kill me.\n", kill_time - interrupt);
  }
}

void deal_with_ctrl_c(int stop_now) {
  kill_time = stop_now;
  if (signal(SIGINT, handle_control_c) == SIG_IGN)
    signal(SIGINT, SIG_IGN); // ignore if parent process was ignoring
}

} // namespace meep

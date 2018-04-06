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

#include "meep.hpp"

using namespace std;

namespace meep {

void fields::finished_working() {
  double now = wall_time();
  if (last_wall_time >= 0)
    times_spent[working_on] += now - last_wall_time;
  last_wall_time = now;
  working_on = was_working_on[0];
  for (int i = 0; i+1 < MEEP_TIMING_STACK_SZ; ++i)
    was_working_on[i] = was_working_on[i+1];
  was_working_on[MEEP_TIMING_STACK_SZ-1] = Other;
}

void fields::am_now_working_on(time_sink s) {
  double now = wall_time();
  if (last_wall_time >= 0)
    times_spent[working_on] += now - last_wall_time;
  last_wall_time = now;
  for (int i = 0; i+1 < MEEP_TIMING_STACK_SZ; ++i)
    was_working_on[i+1] = was_working_on[i];
  was_working_on[0] = working_on;
  working_on = s;
}

double fields::time_spent_on(time_sink s) {
  return times_spent[s];
}

static const char *ts2n(time_sink s) {
  switch (s) {
  case Stepping: return "time stepping";
  case Connecting: return "connecting chunks";
  case Boundaries: return "copying borders";
  case MpiTime: return "communicating";
  case FieldOutput: return "outputting fields";
  case FourierTransforming: return "Fourier transforming";
  case Other: break;
  }
  return "everything else";
}

static void pt(double ts[], time_sink s) {
  if (ts[s]) master_printf("    %18s: %g s\n", ts2n(s), ts[s]);
}

void fields::print_times() {
  master_printf("\nField time usage:\n");
  for (int i=0;i<=Other;i++)
    pt(times_spent, (time_sink) i);
  master_printf("\n");
}

} // namespace meep

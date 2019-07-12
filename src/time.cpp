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

#include "meep.hpp"

using namespace std;

namespace meep {

void fields::finished_working() {
  double now = wall_time();
  if (last_wall_time >= 0) times_spent[working_on] += now - last_wall_time;
  last_wall_time = now;
  working_on = was_working_on[0];
  for (int i = 0; i < MEEP_TIMING_STACK_SZ - 1; ++i)
    was_working_on[i] = was_working_on[i + 1];
  was_working_on[MEEP_TIMING_STACK_SZ - 1] = Other;
}

void fields::am_now_working_on(time_sink s) {
  double now = wall_time();
  if (last_wall_time >= 0) times_spent[working_on] += now - last_wall_time;
  last_wall_time = now;
  for (int i = MEEP_TIMING_STACK_SZ - 1; i > 0; --i)
    was_working_on[i] = was_working_on[i - 1];
  was_working_on[0] = working_on;
  working_on = s;
}

double fields::time_spent_on(time_sink s) { return times_spent[s]; }

double fields::mean_time_spent_on(time_sink s) {
  int n = count_processors();
  double total_time_spent = sum_to_master(times_spent[s]);
  return total_time_spent/n;
}

static const char *ts2n(time_sink s) {
  switch (s) {
    case Stepping: return "time stepping";
    case Connecting: return "connecting chunks";
    case Boundaries: return "copying borders";
    case MpiTime: return "communicating";
    case FieldOutput: return "outputting fields";
    case FourierTransforming: return "Fourier transforming";
    case MPBTime: return "MPB";
    case GetFarfieldsTime: return "getting farfields";
    case Other: break;
  }
  return "everything else";
}

static void pt(double mean[], double stddev[], time_sink s) {
  if (mean[s] != 0) {
    if (stddev[s] != 0)
      master_printf("    %21s: %g s +/- %g s\n", ts2n(s), mean[s], stddev[s]);
    else
      master_printf("    %21s: %g s\n", ts2n(s), mean[s]);
  }
}

void fields::print_times() {
  double mean[Other + 1], square_times[Other + 1], stddev[Other + 1];
  int n = count_processors();

  for (int i = 0; i <= Other; ++i)
    square_times[i] = times_spent[i] * times_spent[i];
  sum_to_master(times_spent, mean, Other+1);
  sum_to_master(square_times, stddev, Other+1);
  for (int i = 0; i <= Other; ++i) {
    mean[i] /= n;
    stddev[i] -= n*mean[i]*mean[i];
    stddev[i] = n == 1 || stddev[i] <= 0 ? 0.0 : sqrt(stddev[i] / (n-1));
  }

  master_printf("\nField time usage:\n");
  for (int i = 0; i <= Other; i++)
    pt(mean, stddev, (time_sink)i);
  master_printf("\n");

  if (verbosity > 0) {
    master_printf("\nField time usage for all processes:\n");
    double *alltimes_tmp = new double[n * (Other+1)];
    double *alltimes = new double[n * (Other+1)];
    for (int i = 0; i <= Other; ++i) {
      for (int j = 0; j < n; ++j)
        alltimes_tmp[i*n+j] = j == my_rank() ? times_spent[i] : 0;
    }
    sum_to_master(alltimes_tmp, alltimes, n * (Other+1));
    delete[] alltimes_tmp;
    for (int i = 0; i <= Other; i++) {
      master_printf("    %21s: %g", ts2n((time_sink)i), alltimes[i*n]);
      for (int j = 1; j < n; ++j)
        master_printf(", %g", alltimes[i*n+j]);
      master_printf("\n");
    }
    master_printf("\n");
    delete[] alltimes;
  }
}

} // namespace meep

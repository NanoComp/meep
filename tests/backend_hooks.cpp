/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
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
*/

/* Reference test for the backend extension-point hook surface
 * (src/meep/backend_hooks.hpp).
 *
 * Asserts:
 *   1. With a transparent counting backend installed (every hook
 *      defers to CPU), numerical results are bit-identical to a
 *      baseline run with no backend.  This is the load-bearing
 *      claim: the hook surface itself never perturbs values.
 *   2. Hooks fire at the expected sites with sensible counts.
 */

#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
#include <meep/backend_hooks.hpp>

using namespace meep;

static double one(const vec &) { return 1.0; }

namespace {

struct hook_counts {
  int init = 0;
  int cleanup = 0;
  int step = 0;
  int sync_to_host = 0;
  int sync_from_host = 0;
  int needs_host_sync = 0;
  int read_point = 0;
};
hook_counts counts;

void test_init(fields *) { counts.init++; }
void test_cleanup(fields *) { counts.cleanup++; }
bool test_step(fields *) {
  counts.step++;
  return false; /* fall through to CPU step path */
}
void test_sync_to_host(fields *) { counts.sync_to_host++; }
void test_sync_from_host(fields *) { counts.sync_from_host++; }
bool test_needs_host_sync(const fields *) {
  counts.needs_host_sync++;
  return false; /* host always considered fresh -> sync_to_host never fires */
}
realnum test_read_point(const fields *, const fields_chunk *, component, int, ptrdiff_t) {
  counts.read_point++;
  return 0;
}

void install_transparent_backend() {
  meep_backend.init = test_init;
  meep_backend.cleanup = test_cleanup;
  meep_backend.step = test_step;
  meep_backend.sync_to_host = test_sync_to_host;
  meep_backend.sync_from_host = test_sync_from_host;
  meep_backend.needs_host_sync = test_needs_host_sync;
  meep_backend.read_point = nullptr; /* leave null so dft_ldos uses the host path */
  counts = hook_counts{};
}

void uninstall_backend() { meep_backend = backend_hooks{}; }

constexpr int n_steps = 100;

double run_sim() {
  grid_volume gv = volone(6.0, 10.0);
  structure s(gv, one);
  fields f(&s);
  /* gaussian pulse: freq=0.7, fwidth=1.0, peak t=0, cutoff at 4 sigmas */
  f.add_point_source(Ex, 0.7, 1.0, 0.0, 4.0, vec(2.5), 1.0);
  for (int i = 0; i < n_steps; i++)
    f.step();
  return f.field_energy();
}

} /* namespace */

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;

  /* 1. Baseline: no backend installed. */
  uninstall_backend();
  const double baseline_energy = run_sim();

  /* 2. Transparent backend: every hook fires but defers to CPU. */
  install_transparent_backend();
  const double hooked_energy = run_sim();
  uninstall_backend();

  /* 3. Numerical results must be bit-identical. */
  if (baseline_energy != hooked_energy) {
    master_printf("backend_hooks: FAIL baseline_energy=%.17g hooked_energy=%.17g (diff=%g)\n",
                  baseline_energy, hooked_energy, baseline_energy - hooked_energy);
    return 1;
  }

  /* 4. Hook fire counts. */
  if (counts.init != 1) {
    master_printf("backend_hooks: FAIL init=%d expected 1\n", counts.init);
    return 1;
  }
  if (counts.cleanup != 1) {
    master_printf("backend_hooks: FAIL cleanup=%d expected 1\n", counts.cleanup);
    return 1;
  }
  if (counts.step < n_steps) {
    master_printf("backend_hooks: FAIL step=%d expected >= %d\n", counts.step, n_steps);
    return 1;
  }
  /* needs_host_sync was reported as false everywhere, so sync_to_host
   * must never have fired. */
  if (counts.sync_to_host != 0) {
    master_printf("backend_hooks: FAIL sync_to_host=%d expected 0\n", counts.sync_to_host);
    return 1;
  }
  if (counts.read_point != 0) {
    master_printf("backend_hooks: FAIL read_point=%d (hook was null)\n", counts.read_point);
    return 1;
  }
  /* needs_host_sync should be queried at least once per readout site. */
  if (counts.needs_host_sync < 1) {
    master_printf("backend_hooks: FAIL needs_host_sync=%d expected >= 1\n", counts.needs_host_sync);
    return 1;
  }

  master_printf("backend_hooks: PASS  (init=%d cleanup=%d step=%d needs_host_sync=%d)\n",
                counts.init, counts.cleanup, counts.step, counts.needs_host_sync);
  return 0;
}

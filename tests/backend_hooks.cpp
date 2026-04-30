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
 * Goals:
 *   1. With no backend installed, every hook call site is a no-op and
 *      the simulation runs exactly as it would without these hooks.
 *      Verified implicitly by the rest of the test suite continuing to
 *      pass; verified explicitly by `baseline_run` here.
 *
 *   2. With a "transparent" backend installed -- one that implements
 *      every hook but defers all real work to the CPU path -- the
 *      simulation produces bit-identical output to the baseline run.
 *      This is the load-bearing assertion for backend authors: the
 *      hook surface itself does not perturb numerical results.
 *
 *   3. The hooks fire at the expected sites: `init` on construction,
 *      `cleanup` on destruction, `step` once per timestep, and the
 *      sync/read hooks are never called when their predicates say no.
 */

#include <stdio.h>
#include <stdlib.h>

#include <meep.hpp>
#include <meep/backend_hooks.hpp>

using namespace meep;
using std::complex;

static double one(const vec &) { return 1.0; }

/* ----- counting "transparent" backend ----- */

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
  return false; /* host is always considered fresh -> sync_to_host never fires */
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
  meep_backend.read_point = nullptr; /* leave read_point null so dft_ldos uses the host path */
  counts = hook_counts{};
}

void uninstall_backend() { meep_backend = backend_hooks{}; }

} /* namespace */

/* ----- the fixture: small 1-D run with a continuous source ----- */

static const int n_steps = 50;

static double run_sim_capture_field(double *out_field_at_origin) {
  grid_volume gv = vol1d(2.0, 20.0);
  structure s(gv, one);
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ez, continuous_src_time(0.3), vec(0.5));
  for (int i = 0; i < n_steps; i++)
    f.step();
  monitor_point pt;
  f.get_point(&pt, vec(1.0));
  *out_field_at_origin = real(pt.get_component(Ez));
  return f.field_energy_in_box(gv.surroundings());
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;

  /* 1. Baseline: no backend installed. */
  uninstall_backend();
  double baseline_field = 0, baseline_energy = 0;
  baseline_energy = run_sim_capture_field(&baseline_field);

  /* 2. Transparent backend installed: every hook fires but defers to CPU. */
  install_transparent_backend();
  double hooked_field = 0, hooked_energy = 0;
  hooked_energy = run_sim_capture_field(&hooked_field);
  uninstall_backend();

  /* 3. Numerical results must be bit-identical -- the hook surface
   *    cannot perturb values when the backend defers to CPU. */
  if (baseline_field != hooked_field) {
    master_printf("backend_hooks: FAIL baseline_field=%.17g hooked_field=%.17g (diff=%g)\n",
                  baseline_field, hooked_field, baseline_field - hooked_field);
    return 1;
  }
  if (baseline_energy != hooked_energy) {
    master_printf("backend_hooks: FAIL baseline_energy=%.17g hooked_energy=%.17g (diff=%g)\n",
                  baseline_energy, hooked_energy, baseline_energy - hooked_energy);
    return 1;
  }

  /* 4. Hook fire counts must match what the call sites promise. */
  if (counts.init != 1) {
    master_printf("backend_hooks: FAIL init fired %d times, expected 1\n", counts.init);
    return 1;
  }
  if (counts.cleanup != 1) {
    master_printf("backend_hooks: FAIL cleanup fired %d times, expected 1\n", counts.cleanup);
    return 1;
  }
  /* step is called once per fields::step() AND once at the start of
   * the (internal) get_field call inside the NaN guard, but only when
   * the backend chose to handle it; here we returned false so the CPU
   * path also runs.  We just check the call count matches the user's
   * step count (the NaN guard's internal calls don't go through step). */
  if (counts.step != n_steps) {
    master_printf("backend_hooks: FAIL step fired %d times, expected %d\n", counts.step, n_steps);
    return 1;
  }
  /* needs_host_sync was reported as false everywhere, so sync_to_host
   * must never have fired. */
  if (counts.sync_to_host != 0) {
    master_printf("backend_hooks: FAIL sync_to_host fired %d times, expected 0\n",
                  counts.sync_to_host);
    return 1;
  }
  if (counts.read_point != 0) {
    master_printf("backend_hooks: FAIL read_point fired %d times (hook was null)\n",
                  counts.read_point);
    return 1;
  }
  /* needs_host_sync is asked at every CPU readout site: at minimum,
   * once for the get_point at the end and once per NaN-guard read in
   * fields::step().  Just check it was queried at least n_steps times. */
  if (counts.needs_host_sync < n_steps) {
    master_printf("backend_hooks: FAIL needs_host_sync queried %d times, expected >= %d\n",
                  counts.needs_host_sync, n_steps);
    return 1;
  }

  master_printf("backend_hooks: PASS  (init=%d cleanup=%d step=%d needs_host_sync=%d)\n",
                counts.init, counts.cleanup, counts.step, counts.needs_host_sync);
  return 0;
}

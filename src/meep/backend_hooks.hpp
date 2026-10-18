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
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/* Backend extension point.
 *
 * Vanilla meep runs entirely on the CPU.  This header declares a small
 * table of function pointers that an external backend library (e.g. a
 * CUDA, ROCm, or vectorized-CPU implementation) may install at load time
 * to redirect hot paths through its own implementation.
 *
 * All entries default to null; null means "fall through to the in-tree
 * CPU implementation", and every call site is written so that a null
 * hook is a no-op.  As a result, a stock build of meep with no backend
 * loaded behaves bit-identically to one without these hooks.
 *
 * State management.  Backends store per-fields state in
 * `fields::backend_state` and per-chunk state in
 * `fields_chunk::backend_state` (both `void*`, owned by the backend).
 * Upstream meep never inspects these; backends are responsible for
 * allocating them in `init` and freeing them in `cleanup`.
 *
 * Suspension.  A backend can be suspended for a single `fields` object
 * by setting `fields::backend_suspended = true`.  The step hook will
 * not fire while suspended; the in-tree CPU step path runs instead.
 * Used by the CW solver, which iterates against host arrays directly.
 */

#ifndef MEEP_BACKEND_HOOKS_H
#define MEEP_BACKEND_HOOKS_H

#include "meep.hpp" /* fields, fields_chunk, realnum, component */

namespace meep {

struct backend_hooks {
  /* Per-sim lifecycle.  `init` is called once after a `fields` object's
   * structure is built and chunks are connected; `cleanup` is called
   * once before the `fields` object is destroyed.  Backends typically
   * allocate/free per-sim shadow state here and stash a pointer in
   * `fields::backend_state`. */
  void (*init)(fields *f);
  void (*cleanup)(fields *f);

  /* Take one FDTD timestep on behalf of the caller.  Returning `true`
   * means the backend handled the step and the in-tree CPU step path
   * should be skipped.  Returning `false` (or leaving this null) falls
   * through to the CPU step path.  Skipped entirely while
   * `fields::backend_suspended` is true. */
  bool (*step)(fields *f);

  /* Sync canonical host arrays with the backend's shadow storage.
   * `sync_to_host` brings the latest field values into the per-chunk
   * `f[c][cmp]` arrays so CPU code can read them; `sync_from_host`
   * pushes any host-side modifications back into shadow storage so the
   * next step sees them.  Backends may treat either as a no-op when
   * not active (e.g. before `init` has been called).  Called by every
   * CPU code path that reads or writes field data. */
  void (*sync_to_host)(fields *f);
  void (*sync_from_host)(fields *f);

  /* Fast single-point read.  Used on the LDOS / point-monitor hot path
   * to fetch one field value without a full sync.  Backends that do not
   * implement this should leave it null; callers fall back to a full
   * `sync_to_host` followed by a direct array read. */
  realnum (*read_point)(const fields *f, const fields_chunk *fc, component c, int cmp,
                        ptrdiff_t idx);

  /* Single predicate.  Returns true iff the host arrays do not reflect
   * the backend's current state and `sync_to_host` should be called
   * before reading them.  Returns false when no backend is active or
   * when the backend's host arrays are already in sync. */
  bool (*needs_host_sync)(const fields *f);
};

/* Single process-global hook table.  Default-initialized to all nulls.
 * The plugin populates this at library load time. */
extern backend_hooks meep_backend;

/* Convenience used at every CPU readout site.  Inlined so it compiles
 * to a single null-pointer test (and nothing else) when no backend is
 * loaded. */
inline void sync_host_if_needed(fields *f) {
  if (meep_backend.needs_host_sync && meep_backend.sync_to_host && meep_backend.needs_host_sync(f))
    meep_backend.sync_to_host(f);
}

/* Read a single cell of `fc->f[c][cmp]`.  If a backend has installed
 * `read_point`, route through it (avoiding a full host sync); otherwise
 * read the host array directly.  The caller is responsible for ensuring
 * `fc->f[c][cmp]` is non-null and (when no read_point hook is installed)
 * for calling `sync_host_if_needed` once before the loop.  Inlined to
 * compile to a single load when no backend is loaded. */
inline realnum read_field_at(const fields *f, const fields_chunk *fc, component c, int cmp,
                             ptrdiff_t idx) {
  return meep_backend.read_point ? meep_backend.read_point(f, fc, c, cmp, idx) : fc->f[c][cmp][idx];
}

} /* namespace meep */

#endif /* MEEP_BACKEND_HOOKS_H */

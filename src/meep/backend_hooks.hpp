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
 * The table is intentionally a flat C ABI so it can be written from a
 * separately-built shared library loaded via LD_PRELOAD or dlopen.  All
 * entries default to null; null means "fall through to the in-tree CPU
 * implementation", and every call site is written so that a null hook is
 * a no-op.  As a result, a stock build of meep with no backend loaded
 * behaves bit-identically to one without these hooks.
 *
 * State management.  Backends store per-fields state in
 * `fields::backend_state` and per-chunk state in
 * `fields_chunk::backend_state` (both `void*`, owned by the backend).
 * Upstream meep never inspects these; backends are responsible for
 * allocating them in `init` and freeing them in `cleanup`.
 *
 * Extending the hook table.  Adding a new function pointer at the end
 * of `backend_hooks` is forward-compatible: the global table is
 * zero-initialized so unset entries simply read as null.  Reordering or
 * removing fields breaks ABI for already-built backends.
 */

#ifndef MEEP_BACKEND_HOOKS_H
#define MEEP_BACKEND_HOOKS_H

namespace meep {

class fields;
class fields_chunk;

struct backend_hooks {
  /* Per-sim lifecycle.  `init` is called once after a `fields` object's
   * structure is built and chunks are connected; `cleanup` is called
   * once before the `fields` object is destroyed.  Backends typically
   * use these to allocate/free per-sim shadow state (e.g. device
   * buffers, communicator handles) and stash a pointer in
   * `fields::backend_state`. */
  void (*init)(fields *f);
  void (*cleanup)(fields *f);

  /* Take one FDTD timestep on behalf of the caller.  Returning `true`
   * means the backend handled the step and the in-tree CPU step path
   * should be skipped.  Returning `false` (or leaving this null) falls
   * through to the CPU step path.  Backends that cannot handle every
   * configuration may return false selectively. */
  bool (*step)(fields *f);

  /* Sync canonical host arrays with the backend's shadow storage.
   * `sync_to_host` brings the latest field values into the per-chunk
   * `f[c][cmp]` arrays so CPU code can read them; `sync_from_host`
   * pushes any host-side modifications back into shadow storage so the
   * next step sees them.  Called by every CPU code path that reads or
   * writes field data. */
  void (*sync_to_host)(fields *f);
  void (*sync_from_host)(fields *f);

  /* Predicates.  The backend owns the truth about its own state;
   * upstream code never inspects backend_state directly. */
  bool (*is_active)(const fields *f);     /* backend is steering this sim */
  bool (*host_is_stale)(const fields *f); /* host arrays need a sync_to_host */
};

/* Single process-global hook table.  Default-initialized to all nulls.
 * Backends populate this once at library load time. */
extern backend_hooks meep_backend;

/* Convenience used at every CPU readout site.  Inlined so it compiles
 * to a single null-pointer test (and nothing else) when no backend is
 * loaded. */
inline void sync_host_if_needed(fields *f) {
  if (meep_backend.host_is_stale && meep_backend.sync_to_host && meep_backend.host_is_stale(f))
    meep_backend.sync_to_host(f);
}

inline bool backend_is_active(const fields *f) {
  return meep_backend.is_active && meep_backend.is_active(f);
}

} /* namespace meep */

#endif /* MEEP_BACKEND_HOOKS_H */

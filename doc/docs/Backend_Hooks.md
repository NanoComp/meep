---
# Backend Hooks
---

`<meep/backend_hooks.hpp>` declares a small extension point that lets an external library plug into meep's hot paths at load time, without forking or patching upstream sources. The intended use is sibling backends (CUDA, ROCm, vectorized-CPU) that ship as separate shared libraries pinned to a specific meep version.

This document describes the contract. The header itself is the canonical source.

## Shape

A single process-global table of function pointers:

```cpp
namespace meep {
struct backend_hooks {
  void    (*init)           (fields *f);
  void    (*cleanup)        (fields *f);
  bool    (*step)           (fields *f);
  void    (*sync_to_host)   (fields *f);
  void    (*sync_from_host) (fields *f);
  realnum (*read_point)     (const fields *f, const fields_chunk *fc,
                             component c, int cmp, ptrdiff_t idx);
  bool    (*needs_host_sync)(const fields *f);
};
extern backend_hooks meep_backend;
}
```

All entries default to null. **Null means "fall through to the in-tree CPU implementation."** Every call site in upstream meep is a null-pointer test, so a meep build with no backend loaded behaves bit-identically to a build without these hooks at all.

A backend is installed by writing function pointers into `meep::meep_backend` at library load time (typically from a `__attribute__((constructor))`).

## Per-sim opaque state

Backends store per-simulation state in `fields::backend_state` (a `void *`) and per-chunk state in `fields_chunk::backend_state`. Upstream meep never inspects these slots; backends are responsible for allocating them in `init` and freeing them in `cleanup`.

## Lifecycle

```
fields::fields(...)             # constructs structure + chunks + connections
  └─ meep_backend.init(this)    # backend allocates per-sim/chunk state
                                # backend stashes pointers in backend_state slots
... user code calls f.step(), f.flux_in_box(), f.add_dft(), etc. ...
fields::~fields()
  └─ meep_backend.cleanup(this) # called BEFORE chunks are deleted, so the
                                # backend can read fields_chunk::backend_state
  └─ chunk teardown
```

`init` is also called from the `fields` copy constructor, so a backend that supports `fields(const fields&)` must be prepared to attach to a freshly-copied object.

## Hook contracts

### `step`

```cpp
bool step(fields *f);
```

Take one full FDTD timestep on behalf of the caller. Returning `true` means the backend handled the step and the in-tree CPU step path is skipped. Returning `false` (or leaving the hook null) falls through to the CPU step path -- useful for backends that handle most configurations but not all.

The step hook is **bypassed entirely** when `fields::backend_suspended == true`. The CW solver sets this flag for the duration of its CG iterations, since it operates against the host arrays directly.

### `sync_to_host` / `sync_from_host`

```cpp
void sync_to_host(fields *f);
void sync_from_host(fields *f);
```

Sync the canonical host arrays (`fc->f[c][cmp]`) with the backend's shadow storage.

- `sync_to_host`: bring the latest field values into the host arrays so CPU code can read them. Called at every CPU readout site (DFT, monitors, integration, dump, CW solver).
- `sync_from_host`: push host-side modifications back into shadow storage. Called after `fields::load`, after `fields::solve_cw`.

Backends may treat either as a no-op when not active. Both hooks are also called via the convenience helper `sync_host_if_needed(f)`, which checks `needs_host_sync` first.

### `needs_host_sync`

```cpp
bool needs_host_sync(const fields *f);
```

Returns `true` iff the host arrays do not reflect the backend's current state and `sync_to_host` should be called before reading them. Returns `false` when no backend is active or when the backend's host arrays are already in sync.

### `read_point`

```cpp
realnum read_point(const fields *f, const fields_chunk *fc,
                   component c, int cmp, ptrdiff_t idx);
```

Fast single-cell read. Used on the LDOS / point-monitor hot path to fetch one field value without triggering a full `sync_to_host`. Backends that don't implement this should leave it null; callers fall back to a sync followed by a direct array read.

### `init` / `cleanup`

```cpp
void init(fields *f);
void cleanup(fields *f);
```

Per-sim setup and teardown. Typical body: allocate device buffers and per-chunk shadow storage, stash pointers in `f->backend_state` and each `f->chunks[i]->backend_state`. The matching `cleanup` releases everything.

## Minimal example

A "transparent" backend that counts hook invocations and defers all work to the CPU. Suitable as a starting skeleton.

```cpp
#include <meep.hpp>
#include <meep/backend_hooks.hpp>

namespace my_backend {

static void on_init(meep::fields *)   { /* allocate device state, stash in f->backend_state */ }
static void on_cleanup(meep::fields *) { /* free device state */ }

static bool on_step(meep::fields *) {
  // run the FDTD step on the device
  // return true to skip the CPU path; return false to defer to CPU
  return false;
}

static void on_sync_to_host(meep::fields *)   { /* device -> host arrays */ }
static void on_sync_from_host(meep::fields *) { /* host arrays -> device */ }
static bool on_needs_host_sync(const meep::fields *) { return false; }

__attribute__((constructor))
static void install() {
  meep::meep_backend.init            = on_init;
  meep::meep_backend.cleanup         = on_cleanup;
  meep::meep_backend.step            = on_step;
  meep::meep_backend.sync_to_host    = on_sync_to_host;
  meep::meep_backend.sync_from_host  = on_sync_from_host;
  meep::meep_backend.needs_host_sync = on_needs_host_sync;
}

}  // namespace my_backend
```

Build it as a shared library, then load it before running meep:

```sh
LD_PRELOAD=libmy_backend.so python -c "import meep; ..."
```

See `tests/backend_hooks.cpp` in the meep source for a complete working example used as a CI guard.

## ABI notes

There is no formal ABI versioning on the hook table. Backends should be pinned to a specific meep version (typically by submodule or distro package). Adding a new function pointer at the end of `backend_hooks` is forward-compatible with already-built backends because the global is zero-initialized; reordering or removing fields breaks ABI.

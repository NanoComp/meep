# Scalable Padé Extrapolation for DFT Monitors

Status: Proposed<br>
Issue: [NanoComp/meep#3217](https://github.com/NanoComp/meep/issues/3217)<br>
Related implementation: [NanoComp/meep#2440](https://github.com/NanoComp/meep/pull/2440)

## Summary

Add optional, bounded-memory Padé extrapolation to Meep's existing DFT monitors. Each MPI-local `dft_chunk` will retain a fixed number of recent weighted field samples while continuing the ordinary DFT accumulation. When a result is requested, Meep will estimate only the uncomputed future tail and add it to the accumulated DFT.

The initial implementation supports only the frequencies configured on the monitor. Arbitrary-frequency evaluation, pure-Padé operation without configured DFT frequencies, and LDOS support are follow-up projects.

A follow-up phase uses the same corrected DFT infrastructure to accelerate adjoint optimization. It replaces neither the adjoint equations nor the material-gradient calculation; it reduces the number of forward and adjoint FDTD timesteps by stopping after the extrapolated frequency-domain quantities have converged.

## Goals

- Parallelize naturally using the existing `dft_chunk` domain decomposition.
- Bound additional storage independently of the number of timesteps.
- Preserve current behavior and allocation when extrapolation is disabled.
- Make corrected values available to all `dft_chunk`-backed consumers, including fields, flux, force, energy, near-to-far, eigenmode coefficients, and adjoint gradients.
- Expose the magnitude of the correction and a full-history versus half-history error estimate.
- Fail safely when a Padé fit is singular, ill-conditioned, or close to a pole.
- Reduce forward and adjoint timesteps for high-Q optimization problems without changing the default adjoint workflow.
- Bound Padé fitting and memory overhead for multi-objective and MPI adjoint calculations.

## Non-goals for the initial implementation

- Replacing or changing the existing Python `PadeDFT` class.
- Evaluating monitors at frequencies not configured during construction.
- Supporting a monitor with an empty frequency list as a pure-Padé transform.
- Supporting `Ldos` or `dft_ldos`, which use separate `Fdft` and `Jdft` accumulators rather than `dft_chunk`.
- Providing a rigorous error bound for nonlinear derived quantities such as flux or force.
- Differentiating through the finite-time Padé fitting algorithm, pole rejection, or fallback branches.
- Enabling automatic Padé stopping for continuous sources or LDOS objectives in the first adjoint integration.

## Public API

Add a keyword-only `pade_samples` argument to each chunk-backed monitor constructor:

```python
sim.add_dft_fields(..., pade_samples=0)
sim.add_flux(..., pade_samples=0)
sim.add_mode_monitor(..., pade_samples=0)
sim.add_force(..., pade_samples=0)
sim.add_energy(..., pade_samples=0)
sim.add_near2far(..., pade_samples=0)
```

Semantics:

- `pade_samples=0` disables extrapolation and preserves existing behavior.
- `pade_samples=k>0` retains at most the latest `k` samples for every local DFT point and component.
- Samples use the monitor's existing decimation schedule, so their spacing is `dt * decimation_factor`.
- The Padé numerator and denominator orders are selected internally from the available sample count.
- The accepted values are `0` or integers greater than or equal to four. Boolean values are rejected.
- Before enough samples are available, getters return the ordinary accumulated DFT and report that extrapolation is not ready.

The existing Python `PadeDFT` remains available and unchanged. Its full-history storage and SciPy implementation make it unsuitable as the backend for this feature, but it remains useful for compatibility and comparison.

## Internal architecture

### Ownership

Store the additional state directly on `dft_chunk`. A chunk already owns:

- one component and its spatial DFT points;
- the configured frequencies and accumulated DFT values;
- component weights, symmetry phases, and interpolation metadata;
- decimation and update scheduling;
- local MPI ownership and persistent-adjoint lifetime.

This placement keeps sampling, storage, fitting, and evaluation local to the rank that owns the field chunk. Tail histories must never be gathered across MPI ranks.

### New chunk state

Conceptually add:

```cpp
size_t pade_samples;
size_t pade_count;
size_t pade_head;
double pade_last_time;
double pade_sample_period;
complex<realnum> *pade_history;

mutable complex<realnum> *effective_dft;
mutable uint64_t state_generation;
mutable uint64_t cached_generation;
mutable pade_status cached_status;
```

The history layout is point-major, with `N * pade_samples` entries. The raw `dft` array remains unchanged. `effective_dft` is a lazy cache containing the raw accumulator plus the extrapolated future tail.

No history or cache allocation is performed when `pade_samples == 0`.

### Sampling

During `dft_chunk::update_dft`, store the coefficient immediately before multiplication by the absolute-time Fourier phase. If the current update is

```text
D_i += a_k exp(i omega_i t_k),
```

the ring stores `a_k`. This preserves the existing spatial interpolation, integration weight, symmetry/Bloch phase, component weight, and quadrature scaling exactly once while keeping the history independent of frequency.

The ring is updated only on actual DFT sampling steps. The ring head is advanced outside the parallel point loop so each thread writes independent spatial entries without contending on shared indexing state.

Electric and magnetic chunks retain their actual sample times; the existing half-timestep staggering must not be removed.

## Numerical method

For `L` available samples, use a near-diagonal Padé approximant with

```text
q = floor(L / 2)
p = L - q - 1.
```

Given retained coefficients `a_0, ..., a_(L-1)`, solve

```text
sum_(j=1)^q Q_j a_(k-j) = -a_k,
```

for `k = p+1, ..., p+q`, with `Q_0 = 1`, and recover the numerator coefficients from the matched power series.

Use a small in-tree partial-pivoted complex Gaussian solver with `complex<double>` scratch arithmetic, including in single-precision builds. Typical systems are about 10 by 10 for the suggested 20-sample history, so this avoids making LAPACK a mandatory dependency. A later SVD or rank-revealing implementation may improve noisy-data robustness.

### Future-only correction

The ordinary DFT already includes every retained sample. The Padé contribution must therefore begin strictly after the latest real sample.

Let the retained series start at time `t0`, and define

```text
z = exp(i omega dt_sample).
```

If `R(z) = P(z) / Q(z)` matches the retained power series, the correction is

```text
exp(i omega t0) * (R(z) - sum_(j=0)^(L-1) a_j z^j).
```

This subtraction is required to prevent double counting. The half-history estimate uses the most recent `floor(L/2)` samples and its own shifted `t0`, but must predict from the same next-sample boundary.

### Failure behavior

Reject a fit when it contains non-finite data, has a pivot below a scale-relative threshold, or evaluates with a non-finite or near-zero denominator. For a rejected point/frequency fit:

- use zero correction, leaving the raw accumulated DFT unchanged;
- record an invalid-fit status and count;
- report an infinite error contribution;
- do not abort the simulation or return NaNs.

The thresholds must be deterministic and covered separately in double- and single-precision tests.

## Corrected-value access

Add a const-facing chunk method such as:

```cpp
const complex<realnum> *dft_values() const;
```

It returns raw `dft` immediately when extrapolation is disabled or still warming up. Otherwise it lazily prepares `effective_dft`, caches the result and status for the current generation, and returns the cache until another sampled update or mutation invalidates it.

Prepare the cache once per chunk before entering downstream parallel loops, particularly the near-to-far OpenMP frequency loop. Lazy initialization must not occur concurrently from multiple worker threads.

All reads of stored DFT values must use the accessor. This includes:

- generic DFT array and HDF5 output;
- flux and complex flux;
- electric and magnetic energy;
- Maxwell stress and force;
- near-to-far propagation;
- eigenmode overlaps and coefficients;
- material-grid adjoint gradients;
- normalization snapshots.

Keep `dft_norm` on the raw accumulator in the initial implementation so `stop_when_dft_decayed` retains its current behavior.

## Error reporting

Add a uniform monitor method without changing existing result types:

```python
monitor.get_pade_error()
```

The returned value should include:

```python
PadeError(
    relative_correction,
    relative_estimate,
    samples,
    ready,
    invalid_fits,
)
```

The two arrays contain one entry per configured frequency. For every component/list representation, compute a relative spatial norm, then report the maximum ratio across those representations:

```text
relative_correction = max_c ||tail_n,c|| / ||accumulated_c + tail_n,c||
relative_estimate = max_c ||tail_n,c - tail_n/2,c|| / ||accumulated_c + tail_n,c||.
```

This is a conservative convergence diagnostic for the underlying DFT fields, not an error bar on nonlinear derived observables.

The C++ layer returns status data. Python may emit at most one `RuntimeWarning` for an unchanged monitor generation. Repeated queries must not repeat the warning. A sampled update or mutation begins a new generation.

A Padé-specific stopping condition is required for the adjoint-acceleration follow-up described below. It must remain opt-in and must not silently replace or alter `stop_when_dft_decayed`.

## Adjoint acceleration

### Opportunity and current limitation

`OptimizationProblem` performs one forward simulation followed by one adjoint simulation for each objective function. Both the forward and adjoint runs currently terminate with `stop_when_dft_decayed`, which observes changes in the raw accumulated DFT norm. `MeepJaxWrapper` uses the same raw stopping condition for its forward and reverse simulations.

Corrected getters alone therefore improve a manually shortened result but do not reduce runtime automatically. Adjoint acceleration requires an explicit Padé-aware convergence controller.

The speedup comes from fewer FDTD timesteps. Padé does not change the asymptotic cost of the material-gradient loop.

### Quantities that must be corrected together

An accelerated gradient is consistent only when Padé is enabled for the complete adjoint dataflow:

1. Forward objective monitors, whose corrected values determine the objective and `dJ/dq` used to construct adjoint sources.
2. The three persistent forward design-region DFT monitors for every design region.
3. The three persistent adjoint design-region DFT monitors created for each adjoint solve.
4. Every direct forward and adjoint DFT read in `material_grids_addgradient`.

Correcting only the objective or only one design-field leg can produce a plausible objective with a systematically truncated gradient.

Prepare corrected pointers once per forward/adjoint chunk before the native material-gradient loops. Do not fit inside the spatial inner loop.

### Driver API

Add the same two opt-in controls to `OptimizationProblem` and `MeepJaxWrapper`:

```python
pade_samples: int = 0
pade_tolerance: Optional[float] = None
```

- `pade_samples=0` preserves the current adjoint registration, allocation, and stopping behavior.
- `pade_samples>0` enables corrected fixed-frequency objective and design-region monitors while retaining the current raw-decay stop when `pade_tolerance` is `None`.
- Setting `pade_tolerance` enables automatic Padé convergence stopping and requires `pade_samples >= 4`.
- Keep `decay_by` or `dft_threshold` as a raw-decay fallback; do not reinterpret either existing parameter.
- Values such as `pade_samples=20` and `pade_tolerance=1e-4` are benchmark candidates, not recommended defaults until validated.

Thread `pade_samples` through every chunk-backed `ObjectiveQuantity.register_monitors` implementation and through `install_design_region_monitors` for both forward and adjoint setup. Reject automatic acceleration when an LDOS objective is present.

### Padé convergence controller

Use one explicit combined stopping predicate. Do not pass multiple callable conditions through the current `until_after_sources` list path: its wrappers close over the loop index, so multiple callables may invoke the final condition rather than their corresponding condition.

Automatic convergence initially supports finite-duration sources only. A check is eligible only after every active monitor has collected a complete `pade_samples` history after `last_source_time`. Continuous or otherwise infinite-duration sources must be rejected instead of waiting indefinitely.

At each eligible check, require every active monitor and frequency to be ready and valid, and compute:

```text
fit_delta = ||F_full - F_half|| /
            max(||F_full||, absolute_floor)

drift = ||F_full(current) - F_full(previous)|| /
        max(||F_full(current)||, ||F_full(previous)||, absolute_floor)
```

Both quantities must be no greater than `pade_tolerance`. `F_half` uses the most recent half-history and predicts from the same future boundary as `F_full`. A negligible zero-valued history is a valid zero tail rather than a failed singular fit.

Exact corrected-vector drift requires a previous corrected snapshot, adding a second `N * Nfreq` cache while automatic convergence is active. Evaluate no more frequently than every `ceil(pade_samples / 2)` new DFT samples so repeated small solves do not erase the timestep savings. Require two passing checks separated by this interval.

The combined predicate records which condition fired:

- `pade_converged`: all Padé checks passed;
- `raw_decay_fallback`: the existing raw DFT condition passed first;
- `maximum_time`: neither convergence condition passed before the cap.

For automatic adjoint optimization, `maximum_time` is a failure and raises rather than silently returning a known-unconverged gradient. Fits that are ready, valid, and based on a complete post-source window may still be exposed diagnostically; all other points fall back to raw accumulated values.

### Nonlinear objective confirmation

Field-level Padé diagnostics do not bound error in arbitrary nonlinear user objectives. Ratios near zero, phase objectives, logarithms, and sharp nonlinearities can amplify small monitor errors in both the objective and `dJ/dq`.

`OptimizationProblem` therefore uses a two-level forward gate:

1. Run the inexpensive field-level full/half and drift checks until every objective and forward-design monitor becomes a convergence candidate.
2. Evaluate the corrected objective arguments, every objective value, and every `dJ/dq`; require finite values and save this first high-level confirmation.
3. Continue for at least `ceil(pade_samples / 2)` new samples and require field-level candidate convergence again.
4. Evaluate the corrected arguments, objectives, and derivatives a second time. Stop only if their absolute-plus-relative drift also passes.

If either confirmation fails, discard the candidate and resume field checks. This bounds expensive eigenmode, near-to-far, user-objective, and Autograd work to convergence candidates instead of every sampling interval.

Automatic Padé stopping requires `OptimizationProblem` objective callables to be pure and deterministic. Side-effecting, time-varying, or externally stateful objective functions are unsupported. Record the number and wall time of high-level confirmations.

The JAX wrapper cannot inspect the eventual outer objective or upstream cotangent during its forward invocation. Its initial Padé stop is therefore monitor-only, explicitly experimental, and does not imply a universal objective or gradient error bound.

### Adjoint-source normalization

`ObjectiveQuantity._adj_src_scale` currently computes the transform of an internally synthesized Gaussian source over `T = sim.meep_time()`. A shorter forward run can therefore change the adjoint-source normalization even when the monitored fields have converged.

Refactor this normalization to integrate over the synthesized source's known complete finite support, or evaluate it analytically, so it is independent of the forward stopping time. If that refactor is deferred, automatic stopping must extend the minimum forward runtime through the complete synthesized-source support.

### Gradient meaning

The accelerated adjoint estimates the derivative of the converged, infinite-time physical frequency-domain objective. It is not the exact derivative of the finite-time Padé-corrected objective because the existing adjoint does not differentiate Padé coefficients, pole rejection, or fallback decisions.

Primary correctness comparisons must therefore use:

- a long-run Padé-disabled adjoint gradient; and
- central finite differences of long-run Padé-disabled objective values.

Finite differences of short Padé-corrected objectives are useful diagnostics but are not the correctness contract.

### Multi-objective memory lifetime

The current `OptimizationProblem` retains the forward design monitors plus every objective's adjoint design monitors until `calculate_gradient`. With `Q` objectives, this can keep `1 + Q` persistent design-monitor sets alive.

For Padé-accelerated adjoints, calculate and store each objective's gradient immediately after its adjoint run, then remove that adjoint monitor set before starting the next objective. This bounds live design storage to one forward and one adjoint set instead of allowing history and corrected caches to scale with `Q`.

### Required run diagnostics

Record one compact result for the forward leg and each adjoint leg before its chunks are removed:

- exit reason;
- simulated time and timestep count;
- retained sample count;
- maximum full-versus-half and corrected-drift ratios;
- readiness and invalid-fit count or invalid-signal ratio;
- Padé fitting/checking wall time;
- for `OptimizationProblem` forward runs, high-level confirmation count and wall time;
- for JAX, an explicit marker that convergence was monitor-only.

## Mutation and persistence semantics

### Continued runs and time discontinuities

Ordinary consecutive `Simulation.run` calls continue both raw accumulation and tail history. Queries do not consume or materialize the extrapolation.

If the next sampled timestamp is not approximately `pade_last_time + pade_sample_period`, clear the history and cache while retaining the raw accumulator. This prevents fitting across `restart_fields`, which resets simulation time while retaining DFT data, or any other discontinuity in sampling.

### Scaling

`scale_dfts` scales only the raw accumulated DFT and clears the Padé history/cache. It must not scale retained history and later mix it with unscaled samples, and it must not materialize predicted future samples before the simulation continues.

### Subtraction

Treat direct subtraction as snapshot algebra:

```text
lhs.raw = effective(lhs) - effective(rhs).
```

Then clear the left-hand history/cache. If accumulation continues, it builds a fresh tail.

### Normalization snapshots

Monitor-level save/load operations and Python `get_*_data`/`load_*_data` APIs serialize a materialized corrected snapshot using the existing fixed-size `N * Nfreq` representation. Importing a snapshot populates raw DFT values and clears history. Histories from separate simulations must never be combined.

### Full simulation checkpoints

Full `fields::dump`/`fields::load` checkpoints have different semantics. They must serialize and restore:

- the raw DFT accumulator;
- `pade_samples`, ring contents, head, and valid count;
- sample period and last sample time;
- cache-invalid state.

The current normalization and checkpoint paths share low-level HDF5 helpers, so the implementation must add separate modes or entry points. If resumable Padé checkpoints are deferred, dump/load must explicitly reject enabled Padé monitors rather than silently saving a materialized prediction that would be double-counted after restart.

## Memory and performance

Additional persistent storage per local chunk for corrected reads is

```text
O(N * pade_samples),
```

independent of elapsed timesteps. A lazy corrected cache may add `O(N * Nfreq)` transient or persistent storage. Automatic Padé convergence additionally retains the previous corrected vector for the temporal-drift check, bringing convergence-cache storage to `2 * N * Nfreq`.

For one three-component design-monitor set in double precision, the approximate additional storage is

```text
3 * Nlocal * (pade_samples + 2 * Nfreq) * 16 bytes,
```

plus bounded fitting scratch. Forward objective monitors also contribute history and caches until adjoint sources are constructed. Streaming each objective's gradient immediately after its adjoint solve keeps only one forward and one adjoint design-monitor set live.

Sampling adds one complex write per monitored spatial point. Fitting is lazy and rank-local, with approximate work

```text
O(N q^3 + N Nfreq q),
```

where `q` is roughly half the retained sample count.

Convergence checks are separated by at least half a history. Padé fitting, drift checks, and high-level objective confirmations must be timed separately; they should remain a minority of accelerated wall time.

Thread `pade_samples` into the pre-run DFT metadata and fragment memory estimator so chunk balancing sees the additional history and cache cost before selecting a layout.

## Implementation sequence

### 1. Freeze contracts and register tests

- Adopt `pade_samples` consistently across Python, SWIG, and C++.
- Finalize correction, failure, mutation, checkpoint, and normalization semantics.
- Add a new C++ test executable to both `check_PROGRAMS` and `TESTS`.

### 2. Implement the numerical kernel

- Add the small complex solver and Padé coefficient construction.
- Add future-tail and full-versus-half-history evaluation.
- Test exact damped exponential sequences and all failure states.

### 3. Add chunk storage and memory accounting

- Add the ring, timing metadata, generation, cache, and status to `dft_chunk`.
- Capture samples inside the existing DFT update loop.
- Include history/cache cost in fragment statistics and chunk balancing.
- Test wraparound, fixed capacity, cache invalidation, decimation, and time discontinuities.

### 4. Complete fixed-frequency integration

- Thread `pade_samples` through Python, SWIG, and C++ constructors.
- Convert every direct raw DFT consumer to the corrected-value accessor.
- Split normalization and checkpoint serialization behavior.
- Implement and test scaling, subtraction, persistence, and restart semantics.

This stage constitutes the minimum viable implementation of issue #3217.

### 5. Add the error API

- Expose readiness, sample count, invalid-fit count, relative correction, and the full-versus-half-history estimate.
- Reduce only compact summaries across MPI.
- Add warning suppression by monitor generation.
- Add corrected-vector drift support and the optional previous-snapshot cache required by automatic convergence.

### 6. Add opt-in adjoint acceleration

- Refactor adjoint-source normalization so it is independent of the shortened forward runtime.
- Thread `pade_samples` through objective monitors and forward/adjoint design-region monitors.
- Replace raw reads in `material_grids_addgradient` with prepared corrected pointers.
- Add one combined finite-source stopping predicate with sparse full/half and drift checks.
- Add the two-level `OptimizationProblem` objective/`dJ` confirmation gate.
- Record per-leg convergence diagnostics and enforce strict maximum-time failure.
- Calculate each objective's gradient immediately after its adjoint run to bound peak memory.
- Add experimental monitor-only JAX integration.

### 7. Deferred extensions

- Design arbitrary-frequency evaluation, including interpolation, ordering, output shape, lifetime, and out-of-band behavior.
- Fix the current zero-frequency-count assumptions before adding pure-Padé operation.
- Design LDOS support separately.

## Test plan

### Numerical unit tests

- One damped complex exponential with an exact geometric tail.
- Two or three modes with distinct complex poles and varied amplitudes.
- Absolute time origin, decimated spacing, and the first predicted `K+1` sample.
- Real and complex field sequences and odd/even history lengths.
- Zero, constant, rank-deficient, near-singular, near-pole, and non-finite inputs.
- Full-history versus half-history error output.

### Chunk and parallel tests

- Ring capacity remains `N * pade_samples` after arbitrarily many timesteps.
- Aggregate capacity across ranks equals the expected globally owned entries without replicated histories.
- Corrected arrays and status agree between one and multiple MPI ranks.
- Electric and magnetic components preserve their half-step phase difference.
- Repeated queries reuse the cache; sampled updates invalidate it; decimation-skipped steps do not.
- Cache preparation is safe before OpenMP-enabled downstream loops.

### Python and consumer tests

- Omitted `pade_samples` and explicit zero are bit-for-bit compatible with current output.
- A short corrected resonator run is closer to a long-run DFT reference than the same uncorrected short run.
- `get_dft_array` and HDF5 output agree.
- Flux, force, energy, near-to-far, eigenmode coefficients, and one persistent adjoint-gradient case use corrected values consistently.
- Error values and statuses are finite and useful on controlled cases without requiring monotonic improvement at every checkpoint.

### Lifecycle tests

- Query, continue, and query again without double counting.
- `restart_fields` clears history at the time discontinuity while retaining raw accumulation.
- Scale, continue, and compare with a raw-history reference.
- Subtract corrected snapshots, clear history, and continue with a fresh tail.
- Test normalization round trips separately from exact full-checkpoint continuation.

### Adjoint correctness tests

- Compare corrected forward objective quantities, objective values, and `dJ/dq` with a long-run Padé-disabled reference.
- Compare corrected forward and adjoint design-region fields before forming the gradient.
- Compare per-frequency and optimizer-facing gradients with a long-run Padé-disabled adjoint.
- Compare directional derivatives with central finite differences of long-run Padé-disabled objectives.
- Treat finite differences of short Padé-corrected objectives as diagnostics only.
- Cover single- and multi-frequency objectives, multiple objective functions, real and complex fields, normalization subtraction, and one persistent-adjoint case.
- Verify that the forward two-level gate performs exactly two high-level confirmations per successful candidate sequence and rejects unstable objective or `dJ` values.
- Verify that each adjoint leg stops and reports diagnostics independently.
- Verify raw-decay fallback, strict maximum-time failure, finite-source rejection, LDOS rejection, and experimental JAX monitor-only behavior.
- Verify source normalization is independent of the shortened forward stop time, including pulse-bandwidth variation.
- Verify streamed multi-objective gradient calculation matches the current retained-all-adjoints result.

### CI coverage

Register the numerical and integration tests so they run in the existing serial-double, MPI-double, and MPI-single workflows. The default-zero path must add no allocation and remain bit-for-bit compatible with the existing test suite.

Keep one small DFT-field adjoint case and one compact multi-objective/eigenmode case in per-commit testing. Put broader cylindrical near-to-far, damping, frequency, and finite-difference sweeps in the scheduled adjoint suite.

### Performance evaluation

Use the waveguide-crossing example as a single-objective benchmark and the six-frequency, two-objective mode converter as the representative multi-run benchmark. Test candidate histories of 20 and 40 samples and candidate tolerances of `1e-3` and `1e-4` at one and two MPI ranks.

Record forward and per-adjoint timesteps, total and FDTD-only wall time, Padé fit/check time, confirmation time, objective and gradient error, cosine similarity, directional-derivative error, exit reasons, invalid fits, and per-rank peak RSS. Use one warm-up and five measured repetitions.

The following are provisional rollout targets for these selected benchmarks, not universal guarantees or merge-blocking mathematical claims:

- at least 40% fewer total forward-plus-adjoint FDTD timesteps;
- at least 1.25x median end-to-end speedup for the single-objective case and 1.5x for the multi-objective mode converter;
- Padé fitting and convergence checks below 20% of accelerated wall time;
- double-precision objective error near `1e-4` and gradient/directional-derivative error near `1e-3`;
- single-precision objective error near `2e-3` and gradient/directional-derivative error near `2e-2`;
- observed RSS increase within the predicted bounded history/cache/scratch allocation plus 10%; and
- no non-negligible invalid-fit signal on an exit labeled `pade_converged`.

Do not document a recommended sample-count/tolerance preset until measurements across representative workloads satisfy or deliberately revise these targets.

## Acceptance criteria

1. `pade_samples=0` preserves current behavior and storage.
2. Additional history remains bounded and distributed across MPI-owned chunks.
3. Exact synthetic tails verify sign, phase, decimation, and the future-only boundary.
4. Invalid fits fall back to raw values with deterministic, observable status.
5. One-rank and multi-rank corrected results agree within precision-appropriate tolerances.
6. A short corrected FDTD run improves against a long-run reference.
7. Every advertised chunk-backed consumer uses the corrected-value path.
8. Scaling, subtraction, restart, normalization, and checkpoint semantics pass dedicated continuation tests.
9. Memory estimation includes Padé history before chunk balancing.
10. Adjoint acceleration corrects objective, forward-design, and adjoint-design DFTs together and validates gradients against long-run adjoints and finite differences.
11. Automatic adjoint stopping uses complete post-source histories, full/half agreement, corrected-vector drift, and required per-leg diagnostics.
12. Multi-objective adjoint gradients are consumed incrementally so peak design-monitor storage does not scale with the number of objectives.
13. The default adjoint path remains bit-for-bit compatible and allocates no Padé state.
14. Arbitrary-frequency, pure-Padé, and LDOS work do not block the fixed-frequency merge.

## Principal risks

- A missed direct access to `dft_chunk::dft` would produce inconsistent corrected and uncorrected outputs.
- Incorrect time origin or failure to subtract retained terms would double-count samples.
- Ill-conditioned fits or poles could dominate quadratic quantities without failure-closed evaluation.
- Treating normalization snapshots as resumable checkpoints would double-count predicted future samples after restart.
- Scaling or subtracting histories with incompatible accumulation epochs would create invalid fits.
- Failing to include history in memory estimation could produce poor chunk layouts despite bounded storage.
- Correcting only the objective, forward design fields, or adjoint design fields would yield an inconsistent gradient.
- Full-versus-half fit agreement can share systematic bias; automatic stopping also requires temporal drift across separated corrected snapshots.
- Field-level convergence does not bound arbitrary nonlinear objective or `dJ/dq` error. `OptimizationProblem` therefore requires bounded high-level confirmations, while JAX remains monitor-only and experimental.
- The accelerated adjoint estimates the converged physical gradient; it does not differentiate the nonlinear Padé fitting and rejection logic.
- Shortening the forward run without decoupling `_adj_src_scale` from `sim.meep_time()` would change adjoint-source normalization.
- Retaining every objective's adjoint history and caches would make peak memory scale with the objective count.
- Continuous sources cannot provide the required post-source history window and must be rejected for automatic Padé stopping.
- Calling expensive fitting or objective confirmation too frequently could consume the intended FDTD speedup.

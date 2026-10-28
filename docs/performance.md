# Workflow performance

Measured on 2026-10-28 using `examples/simple-boat`: the generated 24 ft by
9 ft, 20-degree, 2,500 kg boat, observed sea with seed 42, 100 s dwell,
10 s rests and 5 s ramps. Both revisions compute the full 130 s case with
13,001 motion samples, 1,626 surface/video frames and 5,911 exterior/boundary
field points. Video remains 800 by 800 at 12.5 fps.

## Results

| Command | Baseline | Optimized |
|---|---:|---:|
| Boat and sea generation | Not separately timed | 0.33 s |
| Solve | 619.69 s | 313.70 s |
| Render | 198.27 s | 72.54 s |
| Sum | 817.96 s | 386.57 s |

The measured workflow is **2.12 times faster**, reducing time by 52.7%:
13 min 38 s to 6 min 27 s. The optimized total includes generation.

These are profiled wall times from one complete run per configuration on the
same Windows 11 machine with 32 logical CPUs, Python 3.12.10, NumPy 2.5.4,
SciPy 1.18.1, Matplotlib 3.11.2, Capytaine 3.0.0 and MATLAB R2025b. The shared
MATLAB session was already running. Python cProfile and the MATLAB profiler
were enabled for both solve measurements. Baseline solve and render were
profiled in separate processes and summed; the optimized commands ran in
sequence. Timing claims use wall clocks and stage logs, rather than sums of
inclusive function-profile times. These are case measurements, not CI gates
or speed guarantees for other hardware and inputs.

The baseline is main commit `3f73d30`. The final optimized run used a warm
Green-function table cache. Its first cache miss took 23.4 s in an earlier
run; the final warm lookup was below the stage log's 0.1 s display resolution.
The former implementation regenerated this same table on every solve.

Every named numeric array in `boat.npz`, `sea.npz` and `result.npz` matches the
baseline exactly. The plot and encoded video are byte-identical. The MP4
SHA-256 is `07347a6a2b4c15f82149d9a3e7700cfc2564850e10673fc3770b6b0443f0ae74`.
All 456 pytest tests passed, including the new reuse/invalidation and native
field-matrix equivalence checks. Test temporary directories were removed and
their removal verified.

## Implemented changes

- **Reuse immutable source geometry.** Capytaine's `build_S_matrix` merges
  and copies its source mesh for every call, losing cached face geometry.
  Radiation reconstruction now merges once per spatial chunk and calls the
  same native Green-function evaluator on that mesh. The single-layer matrix,
  source coefficients, quadrature and complex multiplication are unchanged.
- **Control threads around batches.** The previous `n_threads=1` argument
  caused Capytaine to rediscover thread pools for every scalar solve. The
  existing one-thread policy now surrounds the batch, with caller settings
  restored afterward. Coefficient preparation dropped from 15.4 s to 0.8 s.
  Together with source-geometry reuse, radiation reconstruction dropped from
  381.1 s to 127.1 s.
- **Reuse identical numerical tables.** The CLI supplies
  `~/.cache/topwater/capytaine-3.0.0` to Capytaine's parameter-keyed tabulation
  cache. This caches numerical integration tables, not boat coefficients or
  simulated motion. Capytaine documents this facility in its
  [Delhommeau API](https://capytaine.org/stable/developer_manual/api/capytaine.green_functions.delhommeau.html).
- **Cache static rendering work.** An identical hole-free XY projection and
  connectivity reuse their validated triangulation/trifinder. Changed XY,
  changed connectivity and interior holes retain the original validation
  path. Heights and physical contact still update every frame. The static
  HUD also reuses its processed image on the fixed canvas. This extends the
  existing [Agg blitting approach](https://matplotlib.org/stable/users/explain/animations/blitting.html).
- **Reuse connectivity proofs when loading results.** Identical local face
  indices and vertex counts reuse range, duplicate-face and vertex-coverage
  checks. Every frame still receives its own geometric degeneracy check.

No timestep, frequency bank, field grid, finite-band policy, damping model,
motion integration, native audit tolerance, archive schema or rendering
resolution was changed.

## Reproduce and inspect

From the checkout, with the usual shared MATLAB session running:

```bash
uv run --locked --group wecsim python scripts/profile_workflow.py --stage all --matlab-profile
uv run --locked --group wecsim python scripts/profile_workflow.py --stage render
```

The harness prints and retains a unique OS-temporary directory. `--output`
can select another directory outside the checkout. Inputs are linked/copied
into a private case and all outputs publish there. Reports include wall times,
input hashes, Python source/resource hashes, dependency versions, `.prof`
files, command logs and optional `matlab_profile.json`. The published example
is the benchmark input; profiling does not replace its outputs.

Evidence directories under this machine's OS temporary directory:

- Baseline solve: `topwater-profile-sq8ts1jv`.
- Baseline render: `topwater-profile-86li2r22`.
- Intermediate cold-cache run: `topwater-profile-j8w456cw`.
- Final run: `topwater-profile-n11d631o`, including `comparison.json` with
  exact array-member and artifact hashes.

## Further opportunities

1. **Native MATLAB execution is now a leading cost.** Its stage remained
   149.4 s before and after. MATLAB profiling attributes about 104 s to
   simulation and 26.4 s to the separate native excitation reference loop.
   Investigate caching constant FullDir force factors and compiled simulation
   reuse while retaining the applied-force and integration audits. Accelerator
   modes trade compilation cost against execution speed; their benefit needs
   measuring on this actual model, as described by
   [MathWorks](https://www.mathworks.com/help/simulink/ug/what-is-acceleration.html).
2. **Parallel field kernels are promising if roundoff differences are acceptable.**
   A bounded experiment on the actual 33-panel wetted mesh, 128 exterior
   points and 128 frequencies measured median Green-function evaluation times
   of 0.130 s with one OpenMP thread and 0.047 s with four (three repetitions).
   The maximum matrix difference was 2.86e-17, or 1.60e-16 of peak magnitude.
   This is a matrix microbenchmark, not a demonstrated end-to-end speedup or
   full surface/video equivalence result. The shipped default remains one thread.
3. **Rendering now spends time on real geometry and drawing.** The final
   run spent 26.5 s intersecting physical hull/water triangles and roughly
   24 s drawing the Gouraud water mesh. Spatial indexing or compiled contact
   kernels should be evaluated against the existing tangency, hole and
   overtopping fixtures before changing this path.
4. **Incident-spectrum synthesis costs about 11 s.** Coherent aggregation
   of equal headings could reduce repeated exponentials, but would change
   floating-point summation order. Measure field differences before adopting it.

Existing mesh-resolution, irregular-frequency and linear-model applicability
warnings still apply. Performance measurements and exact software parity do
not establish real-boat physical accuracy.

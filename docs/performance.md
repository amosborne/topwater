# Workflow performance

## Full-history integration and pixel-defined water grid

Native integration now uses all available velocity history, without a memory
cutoff option. Water spacing is derived from ten pixels between samples in the
solver-owned view. The solver saves the intended video dimensions in the result;
rendering reads them from the archive.

Complete unprofiled simple-boat solves including result publication:

| Policy | GPU median, three runs | CPU, one run |
|---|---:|---:|
| Exact packet excitation, 10 s memory, 20 cm grid | 9.50 s | 14.59 s |
| Exact packet excitation, full history, 20 cm grid | 16.98 s | 22.25 s |
| Full history, ten pixels between water samples | 18.46 s | 25.60 s |

The controlled 20 cm comparison isolates approximately 7.48 s of additional GPU
workflow time for full history. The final policy produced a 13.799 m square view,
a target spacing of 0.17248 m, 7,396 water vertices per frame and 1,301 saved
frames. CPU and GPU motion matched bit-for-bit. A complete 800 by 800 video
rendered and verified in 62.20 s. Reports and source/input provenance are in
`tests/data/native_preparation.json` and `tests/data/full_history_pixel_grid.json`.

Native excitation now evaluates the actual packet frequencies and directions
coherently, without a frequency/heading interpolation table. Its sources are
reused in reconstruction. Native preparation omits reference-only diffraction
cases and NetCDF/MAT/hydrostatics exports. The optional WEC-Sim reference path
continues to produce its exchange files.

## Simplified CLI and 10 cm surface grid

The current CLI uses `sea`, a named `--dwell` with a 100 s default, and a single
`--device` selector accepting automatic modes or a GPU name. Seeds are recorded
only in `sea.log`. Rendering has no elapsed-time budget.

A complete unprofiled simple-boat solve with `--device gfx1201`, the default
10 cm water grid and 1,301 saved water frames took 19.10 s including publication.
Each water frame contained 21,609 points. Motion matched the saved regression
result to 1e-10. Input hashes, executed source identities and dependency versions
are in `tests/data/simplified_cli_solve.json`. Earlier settings benchmarks below
used the previous 18.75 cm default spacing.

## Configurable solve and output defaults

The calculation controls now live in `native/settings.py` and are exposed through
`SolveOptions` and matching CLI arguments. New solves default to 10 fps, retain
an exact off-cadence endpoint in numeric results, and save their view bounds.
The renderer uses those bounds and omits only the off-cadence movie endpoint.
The following complete, unprofiled simple-boat runs include publication:

| Calculation | Complete solve |
|---|---:|
| Defaults, CPU | 18.15 s |
| Defaults, GPU | 10.92 s |
| 20 ms timestep, GPU, 130.02 s endpoint | 16.14 s |
| 5 ms timestep, 5 mm packet budget, 0.25 m water spacing, GPU | 16.20 s |

CPU/GPU default motion and topology were identical; maximum water-coordinate
difference was 1.38e-15 m. All 549 tests passed, including cold/warm/backoff cases,
option forwarding, legacy archives, and real video encoding. Wheel and rebuilt
source installations passed CPU-only smoke checks without MATLAB or PyOpenCL.

The tighter custom solve preserved the 130.02 s endpoint and rendered 1,301 frames
at 10 fps, ending at 130.0 s, in 42.50 s including verification/publication. This
sea's 5 mm packet budget could not be certified on a 20 ms packet grid; selecting
5 ms succeeded. Packet tolerance governs construction, rather than a post-solve
motion check. Raw evidence is in `tests/data/solve_options_benchmark.json`.

Earlier timings below use the previous 12.5 fps output and surface-view policy;
they document successive optimizations rather than identical current settings.

## Shared CPU/GPU reconstruction pipeline

CPU and GPU now use one coordinator for source caching, lazy transfer handling,
batching, allocation backoff, FFT preparation, field synthesis, and assembly.
Only numerical kernels and memory operations differ. Common code is directly
under `native/`; mirrored `cpu/` and `gpu/` packages contain backend-specific code.

Complete unprofiled simple-boat solves, including result publication and excluding
rendering, measured after each change:

| Successive implementation | Complete solve |
|---|---:|
| Original CPU path | 156.40 s |
| Shared coordinator, cache and Fourier preparation; stock CPU Green function | 152.57 s |
| Potential-only CPU kernels, geometry/buffer reuse, fused source and incident kernels; first run | 19.45 s |
| Final CPU, median of three warmed-process runs | 17.74 s |
| Final GPU, one run | 12.15 s |

The final CPU runs took 17.86, 17.71, and 17.74 s with eight Numba workers.
Motion, clocks, topology and offsets were exactly unchanged. Maximum saved
surface-coordinate differences were below 1.4e-15 m. All 528 tests passed,
including full CPU/GPU cold-cache, warm-cache and forced-allocation-backoff cases.
These are measurements on the current example and hardware, not speed guarantees.

Both wheel and rebuilt-source installations passed distribution verification.
A complete automatic-CPU solve from an isolated base wheel with neither MATLAB
nor PyOpenCL succeeded in 28.28 s; it includes separate first-use startup and
compilation state and Numba 0.68.1 instead of the lockfile's 0.68.0. Its surface
also matched to roundoff. All reports, actual installed source hashes, dependency
versions and output comparisons are in `tests/data/shared_reconstruction.json`.

## GPU spatial batch sweep

Three shuffled complete, unprofiled simple-boat solves per batch size measured
the following medians on the RX 9070 XT. Timings include result publication and
exclude rendering. Two additional warm-up solves were excluded.

| Surface points per batch | Complete solve median |
|---|---:|
| 8 | 14.12 s |
| 16 | 13.05 s |
| 32 | 12.60 s |
| 64 | 12.23 s |
| 128 | 12.23 s |
| 256 | 12.36 s |
| 512 | 12.33 s |
| 1024 | 12.35 s |

The GPU now starts at 64 points: it is the smallest tested batch within 1% of
the fastest median, with half the main FFT workspace of 128. The initial 5%
screen also included 32, but that batch was measurably slower. CPU's default
remains 128. This is a starting policy measured on this case/device, not a
universal optimum or a fixed memory allowance.

Workspaces now use immediate allocations and are reused. If the complete
radiation workspace cannot allocate, partial buffers are released and the batch
halves until it fits; a one-point failure stops the solve. Three further complete
64-point solves with this implementation had a median of 12.22 s. Every result
archive member matched the original sweep exactly. The integration test also
forces allocation failures at 64 and 32 and verifies the 16-point result.
Evidence is retained in `tests/data/gpu_batch_sweep.json`.

Reproduce with `uv run python scripts/benchmark_gpu_batches.py`. The harness
retains reports, logs, plots and result-member hashes outside the checkout;
large result archives are removed after parity verification.

Capytaine mesh-resolution and irregular-frequency checks now run once before
Green-function tabulation against the complete requested bank. They cover cached
sources too; scalar BEM solves skip the completed checks. Complete unprofiled
simple-boat solve times before/after were 12.31 s and 12.22 s, with exactly
identical result archive members. This is a small difference in single-run
measurements, not evidence of a material speedup. Provenance is retained in
`tests/data/bem_check_preflight.json`.

Automatic linear-model applicability assessments have been removed: solve no
longer recomputes exact calm-water restoring at every motion sample, and optional
surface analysis no longer evaluates a wet-set departure threshold. Users assess
physical applicability from the model assumptions and final results.
Complete unprofiled simple-boat solves before and after removal took 14.24 s and
12.31 s, including publication and excluding rendering. Every result archive
member was exactly unchanged. These are single-run measurements; provenance is
retained in `tests/data/linear_check_removal.json`.

## Successive GPU reconstruction improvements

Each row below is a complete unprofiled solve of the unchanged simple-boat case,
including result archive/plot publication and excluding rendering. Runs were
performed successively in the same locked environment on the RX 9070 XT.

| Successive change | Complete solve |
|---|---:|
| Fresh production GPU baseline | 23.29 s |
| Persistent source-bank cache, first population | 20.96 s |
| Persistent source-bank cache, separate-process hit | 18.95 s |
| Shared velocity FFT and interpolation preparation | 17.35 s |
| GPU interpolation and arbitrary-length FFT synthesis | 14.98 s |
| Resident transfer data through GPU synthesis | 13.89 s |

These are single-run measurements, not speed guarantees. The first source-cache
population and baseline vary with initialization; only the hit row demonstrates
reuse. Motion, clocks, topology, and offsets match the saved CPU reference exactly
for every stage. Maximum surface-vertex differences remain below 1.4e-15 m.
Full reports, input/source hashes, and run directories are retained in
`tests/data/successive_gpu_optimizations.json`.

The final wheel was also installed into an isolated environment with the GPU
extra and no MATLAB. Its first complete solve took 24.94 s; a second process
took 13.96 s. This captures the additional first-use startup cost separately
from the successive locked-environment measurements. Both outputs matched the
reference to the same roundoff margin. The wheel contains both OpenCL kernels
and their license.

All Fourier operations retain the original native transform length, complex
phase, interpolation rule, band, and output cadence. GPU Bluestein convolution
padding is an implementation of the same DFT, not future-history padding. FFT
tests cover odd/even lengths, the native 13,001-sample case, phase, virtual zero
transfer rows, and the historical padded policy. No extra GPU FFT package or
system SDK is needed.

## Initial production GPU solve (before the successive improvements)

Native solve now uses an explicit per-solve GPU field backend when the optional
GPU package and a compatible driver are available. A complete unprofiled solve
of the unchanged simple-boat case took **20.93 s** in the locked development
environment. This includes preparation, motion, physical diagnostics, complete
water-surface reconstruction, and publication of `result.npz` and `result.png`.
It excludes video rendering and dependency installation. A first solve from an
isolated GPU-extra wheel installation took **26.23 s**. The isolated environment
uses Numba 0.68.1 versus the lockfile's 0.68.0, and has a separate installation
and startup state; this is a distinct measurement, not a repeat of the locked run.

Paired complete solves from isolated installations of the wheel measured:

| Installation / device | Complete solve |
|---|---:|
| Base package, no PyOpenCL or MATLAB, automatic CPU fallback | 162.72 s |
| GPU extra, no MATLAB, automatic GPU selection | 26.23 s |

That is **6.20x faster** for the complete solve. Both runs used Python 3.12.10,
NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, and Matplotlib 3.11.2. The tested GPU
was the RX 9070 XT with OpenCL driver 3683.0. These are unprofiled single-run
measurements. Inputs, final source hashes, parity results, and profile locations
are retained in `tests/data/production_gpu_solve.json`. All 494 regression and
hardware tests passed. CPU and GPU result motion is identical; surface differences
are at roundoff level. The CPU-path source was verified unchanged between its
isolated package run and final verification.

Two measured opportunities were implemented after initial GPU integration:

- The first production GPU solve took **38.92 s under cProfile**. Its profile
  exposed repeated Fourier exponentials across independent wave-packet carriers
  and CPU multiplication of large GPU potential matrices.
- A bounded per-case packet-basis cache reuses only identical Fourier vectors.
  Every carrier still selects its own cutoff and computes the same diagnostics.
  Kernel matrices and source multiplication now stay on GPU in bounded batches,
  transferring the resulting fields. The resulting full solve took **28.93 s
  under cProfile**, followed by the unprofiled **20.93 s** run above.

Motion matches the published CPU native result exactly. Surface topology and
offsets are identical; maximum saved-vertex differences are around 1e-15 m.
No mesh, frequency bank, time step, damping, warning policy, physical diagnostic,
archive format, or rendering behavior changed. CPU reconstruction remains
available through `--device cpu`; automatic fallback occurs only at initialization.

### Opportunities identified before the successive changes

That production profile pointed to FFT/interpolation and radiation source
preparation as the next useful targets. Interpolation currently repeats the same
frequency search across points/modes, and the velocity FFT is recomputed per
spatial chunk. A shared interpolation map and one per-motion FFT could reduce
this work. GPU FFT synthesis could follow if transfer and precision costs justify
it. Radiation source solutions depend on the mean hull and frequencies, so a
versioned persistent source-bank cache could benefit repeated controller studies;
cache identity and memory/storage cost needed measurement first. These changes
are now implemented and measured in the successive table above.

Rigid-body integration remains about four seconds without profiling. Its
compilation was intentionally deferred at the user's request. Physical-warning
and diagnostic triage is a separate effort; those policies were retained.
Rendering remains outside this solve-only optimization.

For installation, backend selection, reproduction, and supported hardware, see
[GPU reconstruction](gpu.md). `--no-profile` on the profiling harness records
whole-command wall time without instrumentation; omitting it collects cProfile.

## Potential-only GPU experiment

An opt-in OpenCL prototype alongside unmodified Capytaine evaluated the complete
1,920-frequency field workload in **1.12 s including transfers**, versus
**126.28 s** for the current one-thread CPU evaluator. It uses float64 and agrees
to floating-point roundoff. Expanded GPU-assisted surface reconstruction took
**13.46 s**, with identical topology and maximum saved vertex difference
**8.88e-16 m** against the published CPU result. Corrected potential-only CPU and
GPU radiation-only totals are 34.01 s and 33.24 s. Earlier 52-53 s prototype
totals included repeated warning-rendering overhead. The default workflow is unchanged.
See [prototype scope, reproduction, and measurement limits](gpu-prototype.md).

## Native Python solver

The default solver now integrates the complete 130 s example in about **4.05 s**
without profiling. The freshly captured WEC-Sim preparation/integration/audit
stage took 155.4 s. These are different stage boundaries: the WEC-Sim stage also
includes model preparation, a separately evaluated excitation reference, and
audits. They do not imply a 38-fold speedup of the complete workflow.

Maximum native/reference differences over all 13,001 motion samples were
**1.995 mm** CG distance and **0.03369 degrees** orientation, within the 1 cm /
0.1 degree limits. Excitation differed by at most 7.05e-10 in its dimensional
force/moment components; the radiation kernel differed by at most 5.60e-10.
The reference and comparison are retained under `tests/data/`; see
[native solver details](native.md). The revised integrator removes the residual
added-mass delay and uses causal stage-local radiation quadrature.

The native workflow still computes the existing full-complex water surface and
renders every frame. A complete cProfile run from an isolated wheel installation
with no MATLAB Engine installed measured:

| Command | Native | Earlier optimized WEC-Sim |
|---|---:|---:|
| Solve | 175.46 s | 313.70 s |
| Render | 75.08 s | 72.54 s |
| Solve + render | 250.54 s | 386.24 s |

This is a measured 35.1% reduction for solve plus render. Native integration
took 8.4 s under cProfile, including its kernel preparation. The isolated run
used Python 3.12.10, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, and Matplotlib
3.11.2. First-process boat generation also incurred import/JIT setup and took
5.45 s; sea generation took 0.01 s. Timing remains hardware/cache dependent.
The regenerated example contains the verified native outputs.

Profile without MATLAB:

```bash
uv run --locked python scripts/profile_workflow.py --stage all
```

## Earlier WEC-Sim workflow optimization

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
uv run --locked --group wecsim python scripts/profile_workflow.py --stage all --solver=wecsim --matlab-profile
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

Existing mesh-resolution and irregular-frequency warnings still apply. Performance measurements and exact software parity do
not establish real-boat physical accuracy.

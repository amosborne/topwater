# Potential-only GPU prototype

This page records the experiments preceding production integration alongside
unmodified Capytaine 3.0.0. The normal workflow now supports explicit per-solve
GPU selection as described in [GPU reconstruction](gpu.md). Its CPU and GPU
backends now share one production reconstruction pipeline. The benchmark script
retains the paired evaluator experiment; reconstruction modes invoke that shared
production pipeline. The historical radiation-only substitutions described
below have been removed. The standard package does not depend on PyOpenCL unless
the GPU extra is installed.

Run from the checkout on a machine with a suitable OpenCL GPU driver:

```bash
uv run --locked --with pyopencl python scripts/prototype_gpu_surface.py
uv run --locked --with pyopencl python scripts/prototype_gpu_surface.py --points 5911 --frequencies 1920 --reconstruct
```

The prototype was tested with PyOpenCL 2026.1.4 and the AMD RX 9070 XT
(`gfx1201`) on Windows. `--device` accepts a device-name substring; the script
requires an explicitly matching GPU with double precision. Device/runtime/build
failure is an error, with no silent fallback. It installs no system driver/SDK.
Outputs and comparisons are retained in the printed OS temporary directory;
the published example is never overwritten. `--threads` selects the CPU baseline's
thread count (default one). Build/startup time and GPU time including transfers
are reported separately from device-kernel time.

## Implementation

`topwater.solve.gpu.OpenCLPotential` consumes Capytaine's existing source mesh
and Green-function tables. It returns complex single-layer influence matrices
for batches of wavenumbers and surface targets. Capytaine's source strengths,
body BEM solves, and default field evaluator are unchanged.

The translated OpenCL code computes only the scalar potential. It retains the
reference's near-panel Rankine integral, distant-panel approximation, quadratic
table interpolation, asymptotic wave term, surface-panel self term, and runtime
normalization constants. The latter are read from the installed Fortran library
rather than replaced with different rounded constants. Geometry and tables remain
resident; frequency-independent Rankine terms are evaluated once per spatial chunk.
No fast-math option or reduced precision is enabled. Upstream attribution and the
Apache-2.0 license accompany the derived kernel.

The supported domain is deliberately narrow: float64, infinite water depth,
targets exactly at computational z=0, submerged source quadrature points,
`low_freq` singularities, and `scaled_nemoh3` tabulation. This matches the current
surface workflow. It does not provide gradient matrices, finite-depth evaluation,
body-force coefficients, or a general replacement for Capytaine.

The reconstruction experiment temporarily substitutes field routines inside its
own benchmark process, restoring them in `finally`.
It applies GPU influence matrices to the same six Capytaine source columns and
uses the existing FFT/reconstruction code. It then compares every saved surface
vertex and requires identical connectivity and offsets against the published CPU
result. `--gpu-fields` also evaluates diffraction potentials and the physical
incident-wave sum on GPU; body/source solves, FFTs, and packing remain on CPU. This temporary substitution
is benchmark instrumentation, not a production plugin API or thread-safe adapter.

## Validation

### Total reconstruction comparison

| Reconstruction backend | Total wall time |
|---|---:|
| Stock CPU evaluator, eight threads | 56.44 s |
| Potential-only CPU evaluator, eight threads | 34.01 s |
| Potential-only GPU evaluator, radiation | 33.24 s |
| GPU radiation, diffraction, and incident synthesis | 13.46 s |

These are single full-case runs, with identical surface topology and only
roundoff-level vertex differences. The compiled CPU experiment uses the same
potential-only formulas and frequency-independent geometry reuse as the GPU
experiment. Compilation is warmed before reconstruction timing. At the complete
reconstruction level, the two potential-only implementations are close; the
large evaluator-only speedup below is not an end-to-end speedup.

All rows use production warning aggregation. Earlier prototype totals around
52-53 s included accidental repeated warning rendering and have been retired.
The expanded GPU run's maximum saved vertex difference was 8.88e-16 m, with
identical topology. These are unprofiled single-run total reconstruction times.

```bash
uv run --locked --with pyopencl python scripts/prototype_gpu_surface.py --reconstruct --reconstruct-only --gpu-fields
```

Reproduce the potential-only CPU run without a GPU runtime:

```bash
uv run --locked python scripts/prototype_gpu_surface.py --threads 8 --reconstruct-cpu --potential-only-cpu --reconstruct-only
```

### Evaluator experiment

The complete comparison used 1,920 frequencies, 5,911 field points plus 23 special
waterline/panel points, and all 33 source panels:

| Measurement | Result |
|---|---:|
| CPU evaluator, one thread | 126.285 s |
| GPU evaluation including transfers | 1.121 s |
| Evaluator speedup | 112.65x |
| Maximum potential difference | 2.93e-15 |
| Relative RMS potential difference | 3.51e-15 |

All surface connectivity and offsets matched the saved CPU result exactly.
The original CPU reconstruction was about 150 s under cProfile. The total
reconstruction table above uses unprofiled runs. The evaluator timings are paired unprofiled
measurements in one process. They compare the existing general CPU evaluator
with a GPU implementation that also skips derivatives and reuses Rankine terms;
the improvement cannot be attributed to GPU hardware alone.

The report is retained in `tests/data/simple_boat_gpu_comparison.json`, including
device/driver identity, source hashes, counts, initialization, and kernel timing.
The full reconstruction still uses CPU source solves, NumPy source multiplication,
FFTs, and packing. The expanded experiment moves diffraction potentials and
incident synthesis to GPU as described above. No faster video-rendering claim
is made. These are one-run prototype measurements on the stated hardware.

Opt-in GPU tests cover near-panel, vertex, edge, free-surface self, and far-field
evaluation; resident-point reuse and invalidation; and invalid inputs:

```bash
TOPWATER_TEST_GPU=1 uv run --with pyopencl python -m pytest tests/test_gpu_potential.py
```

Normal pytest skips these hardware tests. Ordinary source and native-solver tests
continue to require neither a GPU nor PyOpenCL. Performance numbers describe the
tested device and workload, not a guarantee for other OpenCL implementations.

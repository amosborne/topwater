# GPU reconstruction and video drawing

The native solver supports GPU reconstruction with unmodified Capytaine 3.0.0.
GPU acceleration covers radiated and scattered water-surface fields, incoming-wave
spectra and radiation FFT synthesis. Boat integration, body BEM solves on cache
misses and archive publication remain on CPU. Video drawing uses hardware
WebGPU when available, with native Skia CPU rasterization otherwise.

## Install and run

The normal installation works on CPU. To enable GPU discovery, install the
optional extra:

```bash
uv sync --locked --extra gpu
uv run --extra gpu topwater solve examples/simple-boat-native --dwell 100 --device auto
uv run --extra gpu topwater render examples/simple-boat-native
```

For an installed distribution, use `pip install "topwater[gpu]"`. The extra adds
PyOpenCL and wgpu-py. A working OpenCL GPU driver with double-precision support must already
be installed. Kernels compile through that driver; no separate compiler or HIP
SDK is required. The tested configuration is Windows, AMD Radeon RX 9070 XT
(`gfx1201`), driver 3683.0, PyOpenCL 2026.1.4. Other OpenCL devices need their own
numerical/performance validation before being described as tested combinations.

| Selection | Behavior |
|---|---|
| `--device auto` (default) | Select an available float64 GPU, otherwise use CPU during initialization. |
| `--device cpu` | Run the shared reconstruction pipeline with CPU kernels. |
| `--device gpu` | Require GPU initialization; report an error if unavailable. |
| `--device "GPU name"` | Select a GPU by name. An unavailable or ambiguous explicit selection is an error; available names appear in the error. |

The selected backend and driver appear in console output and `result.log`.
Automatic selection prefers dedicated GPUs over shared-memory integrated GPUs.
GPU failures after initialization fail the solve and preserve previous published
outputs; they are not silently retried on CPU. Device settings apply to the native
solver; `--solver=wecsim` uses the same reconstruction pipeline with CPU kernels.

The Python API accepts the same choices:

```python
from topwater.solve import simulate

result = simulate(boat, sea, 100., device='auto', influence=model)
```

## Implementation and compatibility

Per-solve field resources are passed explicitly into reconstruction. Production
uses no monkeypatches or global backend assignment. GPU resources are released
after the solve, and CPU imports require neither PyOpenCL nor MATLAB.

`native/fields.py` owns source caching, source-to-field coupling, batching,
allocation backoff, and cleanup for both devices. `native/fourier.py` owns the
shared motion-spectrum cache, interpolation map, and transform-size policy.
Backend-specific construction and numerical operations are mirrored under
`native/cpu/{backend,potential,fourier}.py` and
`native/gpu/{backend,potential,fourier}.py`; OpenCL allocation helpers and kernels
also live under `gpu/`. Common code stays directly under `native/`.
The CPU uses NumPy/SciPy FFTs and potential-only Numba kernels; the GPU uses
OpenCL kernels. Capytaine still owns all body BEM solves and hydrodynamic
coefficients. Its installation is unmodified. Topwater maintains the specialized
CPU and GPU potential translations, which consume the same Capytaine tables.

The Capytaine adapter and GPU kernels live under `src/topwater/solve/`.
The adapter uses Capytaine's existing source strengths and tabulation; the scalar
potential kernels retain its near-panel, surface-self, interpolation, and
asymptotic formulas. The derived code includes upstream Apache-2.0 attribution.
Compatibility is pinned to Capytaine 3.0.0. The supported numerical domain is the
current workflow: infinite depth, zero forward speed, z=0 computational surface,
float64, low-frequency singularity treatment, and scaled Nemoh3 tables.

Geometry and tables remain resident. Frequency-independent Rankine terms are
reused. Potential matrices and source multiplication stay on GPU in bounded
batches; only the resulting fields return to CPU. There is no mesh refinement,
frequency reduction, timestep change, damping change, or runtime CPU/GPU accuracy
comparison. Reference comparisons run in tests. Capytaine mesh-resolution and
irregular-frequency checks run once at the start against the full coefficient,
wave-sideband, and surface-radiation frequency bank, before tabulation or BEM
solves. The same warnings are produced on source-cache hits. Per-problem
wavelength checks are disabled afterward; malformed-input and matrix-NaN checks
remain enabled. Linear-model applicability is left to users inspecting results.

The production CPU backend reuses geometry and host buffers, fuses potential
evaluation with source multiplication, and computes signed incident waves in a
compiled kernel. It uses up to eight Numba workers, respecting a lower current
Numba thread limit and restoring the caller's thread state afterward. First use
may compile kernels; subsequent processes reuse Numba's on-disk compilation cache.

## Reproduce full solve measurements

Radiation source banks are reused between solves under
`~/.cache/topwater/radiation-sources-v1`. `TOPWATER_SOURCE_CACHE_DIR` can redirect
this cache. Keys include numerical geometry, mode motion fields, frequencies,
solver settings, physical parameters, and library versions. Entries contain exact
numeric arrays (about 6 MB for the example), published by atomic replacement.
Invalid entries rebuild; cache-write failures leave computation available.
Changing waves or a controller can reuse sources when the mean hull is unchanged.

Each reconstruction shares one velocity spectrum and interpolation map. The GPU
performs double-precision Fourier synthesis at the original, potentially odd,
native FFT length. An OpenCL radix-2/Bluestein implementation handles arbitrary
lengths without a separate SDK or FFT dependency. Its convolution padding is
internal to the DFT algorithm; it does not pad the physical history, alter the
frequency support, or add future motion. Transfer matrices stay resident and only
the requested surface-time samples return to CPU. First-use kernel compilation
and source-cache misses can make the first solve slower than subsequent solves.

The measured GPU starting batch is 64 surface points. Before radiation field
evaluation, the solver reserves reusable FFT, transfer, source, and potential
matrix buffers using PyOpenCL's immediate allocator. For an actual transform
length `N`, Bluestein uses `M = next_power_of_two(2*N - 1)`; its two complex128
work arrays require `32*M*max(batch_points, 6)` bytes, also accommodating the
six-channel velocity transform. Other buffers are sized from the
case's frequency count, mesh panels, and requested output samples. A failed
allocation releases partial workspaces and halves the batch, stopping with an
error if one point cannot fit. Other computation failures propagate normally.
The workspaces are released after radiation, before later field stages.

```bash
uv run --extra gpu python scripts/profile_workflow.py --stage solve --device gpu --no-profile
uv run --extra gpu python scripts/profile_workflow.py --stage solve --device gpu
```

The first command records complete wall time without profiler overhead. The
second also writes a cProfile file and source hashes, including GPU kernels.
Both operate on a private copy of the simple-boat case in OS temporary storage.
See [performance](performance.md) for measured totals and remaining opportunities.

The opt-in assembled regression runs complete simple-boat CPU and GPU solves with an
empty source cache and again with a cache hit, verifies cache/device consistency,
and publishes/reloads the result and plot. Like-for-like historical motion-model
comparisons with the WEC-Sim reference live separately in `tests/test_native_solver.py`.
The assembled regression also forces allocation failures
and verifies automatic reduction to 16 points on both devices. Enable it with
`TOPWATER_TEST_INTEGRATION=1` when running `tests/test_solve_integration.py`;
also set `TOPWATER_TEST_GPU=1` to include the GPU cases.
Hardware tests select a GPU automatically; `TOPWATER_GPU_DEVICE` selects a
device-name substring for the potential and Fourier component tests.

## Video drawing

The render workflow shares its scene, camera, physical contact, palette, HUD and
encoding code under `topwater/render/`. Mirrored `cpu/backend.py` and
`gpu/backend.py` implementations consume that same scene. CPU uses Skia's native
rasterizer; GPU uses wgpu-py's offscreen hardware drawing. No window or display
canvas is created. GPU emulation adapters are deliberately excluded.

The GPU extra enables automatic physical-GPU discovery for `topwater render`.
Without it, the base installation renders on CPU without a graphics driver.
The render log names the selected implementation and device. Windows hardware
drawing is tested on AMD Radeon RX 9070 XT through Vulkan; wgpu also supports
Direct3D 12 on Windows, Metal on macOS and Vulkan on Linux when compatible
hardware/drivers are installed; see [wgpu platform backends](https://wgpu-py.readthedocs.io/en/stable/backends.html).
Those other combinations have not been validated here. Linux installations need
the Fontconfig and OpenGL shared libraries required by the Skia wheel, even when
only its CPU rasterizer is used; no GPU context is created on that path. See
[Skia installation](https://skia-python.github.io/skia-python/install.html).

The Python API accepts `render_result(result, path, device='cpu')` or
`device='gpu'`; default `device='auto'` uses hardware GPU when initialization
succeeds, otherwise native CPU. Explicit GPU requests fail if no physical adapter
is usable. Errors after drawing begins fail the render and preserve the previous
video. The CLI adds no rendering-device option.

Both backends preserve saved geometry, view bounds, resolution and cadence.
Their antialiasing differs at triangle/line edges: GPU uses four-sample coverage;
Skia uses native analytic coverage. Small pixel differences are expected, while
the physical scene and meter-scaled HUD remain the same. Complete benchmarks
and visual comparison evidence are in [render performance](render-performance.md).

# Current workflow accounting

Code follows the CLI: `topwater/boat`, `sea`, `solve`, and `render`.
Shared numeric inputs, archive contracts and constants stay at the package root;
shared physical geometry lives in `geometry/`. There are no compatibility facades.
Solve has one workflow with native/WEC-Sim motion backends and CPU/GPU field
backends. Render has one scene/encoding workflow with native CPU and hardware GPU
rasterizers. The solver saves the view, cadence and video size in `result.npz`.

## CLI inputs

| Step | Input | Default | Purpose |
|---|---|---|---|
| All | DIRECTORY | Current directory | Input/output location |
| All | -h / --help | — | Show help |
| boat | --length | Required | Length, metres |
| boat | --beam | Required | Beam, metres |
| boat | --deadrise | Required | Transom deadrise, degrees |
| boat | --mass | Required | Operating mass, kilograms |
| sea | --record | Required | Observed sea-state ID |
| sea | --seed | Random; recorded only in sea.log | Reproducible realization |
| solve | --dwell | 100 s | Full-strength wave interval |
| solve | --rest | 10 s per side | Quiet intervals |
| solve | --ramp | 5 s per side | Wave ramps |
| solve | --time-step | 0.01 s | Packet and motion sampling |
| solve | --packet-tolerance | 0.01 m | Packet-envelope construction budget |
| solve | --solver | native | Python or WEC-Sim motion backend |
| solve | --device | auto | auto, cpu, gpu, or a GPU name; native field backend |

Render has only the common directory/help inputs. There are no FPS, render-stride,
surface-spacing, radiation-memory, GPU-substring or render-budget CLI options.

## Internal numerical and output policy

| Item | Current value or rule | Owner |
|---|---|---|
| Hull/deck mass | 0.1 L_ft^3 + 70 kg | boat/generation.py |
| Payload distribution | 50/50 at (0.05L,0,0.03L) and (0.50L,0,0.20L) | boat/generation.py |
| Hull shape/integration | Eight inherited profile stations; 41 structural integration stations | boat/generation.py |
| Sea realization | Eight components per observed band; 720-angle MEM grid; uint32 PCG64 seed | sea/states.py |
| Video size | 800 x 800 pixels; saved in the result | solve/settings.py |
| Water sample pitch | At most 10 output pixels; metres derived from view span | solve/settings.py, solve/surface.py |
| Saved cadence | Every 10 accepted steps; also archive the exact endpoint | solve/settings.py, solve/sampling.py |
| Video endpoint | Omit an off-cadence final endpoint | render/geometry.py |
| View padding | One beam per side of the hull-trajectory square | solve/settings.py |
| Native radiation memory | Full case history | solve/native/dynamics.py |
| Motion coefficient spacing | 0.0159155 Hz | solve/settings.py |
| Motion coefficient coverage | Starts at 0.0159155 Hz; through at least 1.90986 Hz, extended to packet coverage | solve/settings.py |
| Kernel frequency integration | Analytical integral of piecewise-linear B/omega; no quadrature count | solve/native/radiation.py |
| Water radiation frequencies | Every velocity FFT bin: spacing 1/(N dt) Hz through the highest bin | solve/fourier.py |
| CPU/GPU starting batch | 128/64 points; halve and retry allocation failures | solve/fields.py, solve/cpu, solve/gpu |
| CPU potential workers | Up to eight available Numba threads | solve/cpu/potential.py |
| Packet basis cache | At most 64 MiB | solve/packet_basis.py |
| GPU matrix workspace | Minimum of 256 MiB, device allocation limit, and device memory / 8 | solve/gpu/potential.py |
| Legacy archive framing | 1.45 times trajectory span when saved bounds are absent | render/frames.py |
| Video drawing | Native Skia CPU or hardware WebGPU; indexed geometry | render/cpu, render/gpu |
| Raster output | RGBA8; GPU four-sample coverage; CPU native antialiasing | render/cpu, render/gpu |
| CPU HUD storage | Cached native premultiplied BGRA; one final RGBA conversion | render/cpu/backend.py |
| Water colour/gauge range | Symmetric maximum saved elevation, rounded up to 0.1 m | render/workflow.py |
| Presentation | Cached HUD/graticule; fixed style scaled from an 800-pixel design | render/hud.py |
| Video encoding | H.264, 3 Mbit/s, one thread, even pixel dimensions | render/workflow.py |
| Packet sampling | At least 64 samples and 16 samples per carrier period/ramp | solve/packets.py |
| Startup quietness scale | 0.001 of force/velocity scales; elevation uses packet tolerance | solve/native/excitation.py, solve/preparation.py |
| Incident certification budget | 480 s | solve/preparation.py |
| Physical defaults | Gravity 9.81 m/s^2; seawater density 1025 kg/m^3 | constants.py |

The old incident-frequency/beam eligibility screening and repeated sideband
warning have been removed. Structural input checks, sampling requirements, and
the once-at-start Capytaine mesh/irregular-frequency checks remain.

The motion coefficient grid remains a chosen numerical policy, independent of
reconstruction. Analytical kernel integration is exact for its chosen
interpolant, not an exact continuum hydrodynamic solution. Reconstruction uses
all raw, unpadded velocity FFT bins; endpoint behavior remains unchanged.
The optional WEC-Sim reference retains its historical ten-second memory and
coefficient-table/export policy. The hull-shape and mass assumptions are detailed
in [the boat model](model.md#standard-boat-generator).

## Prepared and cached items

| Item | When | Policy |
|---|---|---|
| Mean flotation, wetted hull, mass/restoring matrices | Before integration | Frozen linear boat model |
| Packet Fourier coefficients | Before integration | Construct to envelope budget; lazily reuse harmonic basis |
| Radiation damping and infinite-frequency added mass | Before integration | Fixed motion coefficient grid |
| Radiation kernel and RK-stage history weights | Before integration | Analytical kernel over full duration; causal time quadrature |
| Coherent packet excitation | Before integration | One combined diffraction solve per active packet frequency |
| Packet diffraction sources | During excitation preparation | Reused for water reconstruction |
| Green tables and mesh geometry | Initialization / batch use | Reuse unchanged geometry |
| Water-radiation BEM source bank | Reconstruction | Exact disk cache or missing solves; share across spatial batches |
| Velocity FFT and bin identity | Reconstruction | Once per waveform/settings; share across batches |
| Water grid | After motion | Derive from solved view and video pixels |
| MATLAB exchange files/reference-only cases | WEC-Sim only | Absent from native preparation |
| Video frame adapters | Render startup | Shared saved-result or explicit-mesh scene contract |
| Indexed mesh, HUD and drawing buffers | Render startup / reuse | Common scene, device-specific rasterization |

## Complete workflow measurements

| Configuration | Complete solve |
|---|---:|
| GPU, cold source cache | 19.29 s |
| GPU, warm source cache (two-run median) | 13.69 s |
| CPU, warm source cache | 30.28 s |

These include result publication and exclude video rendering. They use full
history, analytical kernel integration and exact FFT-bin water transfers on the
130 s simple-boat case. CPU/GPU motion and topology agree exactly; water vertices
differ by at most 1.24e-14 m. The record has 13,001 samples and 6,501 real FFT bins,
spaced 0.007691716 Hz through 49.996154 Hz; these follow from N and dt.

Complete rendering, encoding and full decode verification measured **53.02 s
native CPU / 29.03 s hardware GPU**, excluding 1.46 s archive loading from both.
GPU emulation is not used for CPU rendering. Matplotlib remains for static plots
and some geometry/font utilities, not video drawing.

Changing only kernel quadrature to analytical integration changes motion by at
most 5.27e-7 m / 1.93e-5 degrees. Independent adaptive quadrature tests check zero,
tiny and long lags. The changed default model differs from historical WEC-Sim
by **0.09070 m / 0.4353 degrees**, outside the earlier 1 cm / 0.1 degree target:
exact packet excitation and full history differ from the reference's interpolated
excitation and ten-second memory. Like-for-like historical-model tests retain
the original tolerance; this is separate from the changed default model.

`examples/simple-boat-native/` and `examples/simple-boat-wecsim/` have identical
boat/sea inputs and independent results, plots, videos and logs.

Evidence: [solve benchmark](../tests/data/analytical_radiation_fft_bins.json),
[render benchmark](../tests/data/native_render_backends.json),
[render performance](render-performance.md).

## Final verification

The final reviewed tree passed **583 tests**, with GPU and full integration cases
enabled. The unique pytest temporary base was removed and checked. Both wheel
and source distributions passed isolated installation/resource checks, including
native CPU drawing without MATLAB or GPU packages. A fresh CLI solve in a private
copy took 13.58 s and reproduced every published native result member exactly.

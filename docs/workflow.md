# Workflow

A case is a directory containing `boat.npz` and `sea.npz`. All archives use the
[numeric file formats](formats.md); no generator replay is required to solve.
Commands default to the current directory when DIRECTORY is omitted.

## Generate a boat

```powershell
uv run topwater boat DIRECTORY --length 7.3152 --beam 2.7432 --deadrise 20 --mass 2500
```

All four flags are required: hull length and maximum beam in metres, transom
deadrise in degrees, and total operating mass in kilograms. Length excludes
outboards and brackets. This example is 24 ft by 9 ft. Generation replaces
`boat.npz` and `boat.log`; it needs no MATLAB or Capytaine. A failed generation
preserves the previous boat archive. Existing sea, result and video files are
left unchanged and may belong to the previous boat until solve/render are run.

The surface is the convex envelope of Marine Gyro's eight warped-V control
stations, with its upper closure omitted from physical panels. Mass consists
of distributed hull/deck structure and two equal payload lumps. See the
[standard boat model](model.md#standard-boat-generator) for formulas, coordinates
and limitations. Geometry, mass, inertia and stable dry-opening flotation are
checked before publication. The original nominal family is 20–26 ft with
18–25 degree transom deadrise; generation reports extrapolation outside it.

## Generate a sea

```powershell
uv run topwater sea DIRECTORY --record RECORD_ID --seed 42
```

This replaces `sea.npz` and `sea.log`. Omitting the seed draws one and records it only in `sea.log`.
Use `topwater.sea.states.catalog()` to inspect the packaged observed records.
Save the generated arrays for durable replay across generator changes.
Generation retains observed band energy and directional information; phases are
synthetic, so it reproduces a spectrum rather than the measured elevation trace.
The original carrier arrays are not rewritten by solve.

## Solve motion and water surface

```powershell
uv run topwater solve DIRECTORY --dwell 100 --rest 10 --ramp 5
```

`--dwell` is the positive full-strength interval in seconds (default 100). Total time is
`2*rest + 2*ramp + dwell`. Ramp must be positive, rest nonnegative, and all three
durations must align to `--time-step` (default 0.01 s). One water frame is saved every 10 motion steps.
At the default 10 ms timestep this produces 10 fps. The result also saves the exact
endpoint; rendering omits it if it falls between regular video frames. Total
duration does not have to be a multiple of the video interval.

The everyday calculation controls are CLI arguments. Developer controls remain
fields of
[`SolveOptions`](native.md#calculation-options). For example:

```sh
uv run topwater solve DIRECTORY --dwell 100.02 --time-step .005 --packet-tolerance .005
```

| Option | Default | Effect |
|---|---|---|
| `--time-step` | 0.01 s | Packet sampling and motion integration step |
| `--packet-tolerance` | 0.01 m | Packet envelope fitting budget, allocated by carrier amplitude |

Frame rate is derived as `1 / (10 * time-step)`. Render stride, one-beam view
padding and numerical developer controls are configured through `SolveOptions`
in the Python API. Native integration uses the entire available velocity history
without a memory cutoff. The optional WEC-Sim reference backend retains its
historical fixed timestep, ten-second memory and quadrature. All selected controls are recorded in `result.log`.

Water samples are spaced at most ten video pixels apart, using the internal
output dimensions (800 by 800 pixels) and the solved world view. Physical spacing
is `10 * min(view_span_xy_m / video_size_px)`, and the regular grid fits that
spacing to its domain. For the simple boat this is approximately 17 cm.
There is no surface-spacing or radiation-memory CLI option. Native radiation
uses the full history since simulation start, including each provisional RK stage.

The solver saves `video_size_px`, view bounds and cadence in `result.npz`.
Rendering reads those values, so it uses the resolution the solver planned for.
The 10 ms motion step and 1 cm packet budget preserve the established calculation
controls. Radiation discretization choices remain developer controls, not
whole-solution accuracy guarantees. Native excitation uses the packet frequencies
and directions directly.

Sea phases refer directly to the dwell midpoint. Each carrier's target is
multiplied by a quintic rest/ramp/dwell/ramp/rest envelope. Fourier sidebands use
their own dispersion and retain their carrier's direction. See [the model](model.md).

The default backend is `native`, implemented in Python without MATLAB/Simulink.
Use `--solver=wecsim` for the optional MATLAB reference backend. Both produce the
same archive format. The [native API](native.md) supports external Python models.
Native reconstruction selects its device with `--device auto|cpu|gpu|"GPU name"`; `auto`
uses GPU when available and CPU otherwise. Install the optional GPU extra and
driver as described in [GPU setup](gpu.md). Use `--device "GPU name"` to select a GPU directly by name.

Solve writes `result.npz`, `result.png` and `result.log`. The result contains the
computed motion, physical hull and total incident/diffraction/radiation surface.
The plot combines six motion histories with the ORIGINAL sea's spectrum; packet
sidebands and reconstructed water do not supply that spectrum.

The polar chart uses 15 degree by 0.025 Hz arrival-direction/frequency bins and
`sqrt(0.5*sum(amplitude**2))` RMS elevation. Default bow heading points up; a dotted
line marks the saved initial bow. Its linear cividis colorbar includes zero bins.
Carrier/spectrum metadata is not stored in the result archive or video.

The output archive and plot are prepared before publication. Handled publication
failures attempt to restore the previous pair. Abrupt process termination between
replacements can leave an inconsistent pair; a log without successful END means
the invocation is incomplete. Input boat/sea files are not changed. An existing
video is left unchanged and may belong to an earlier result.

## Render a saved result

```powershell
uv run topwater render DIRECTORY
```

This reads only `result.npz` and replaces `video.mp4` and `video.log`; it does not
run a solver or require MATLAB. The saved result's hull is authoritative, even
if `boat.npz` was edited later. Rendering computes the posed hull's waterline.
The video uses the solver's saved square view, metre graticules, telemetry
and a water-elevation gauge. It has no sea-spectrum inset. Legacy archives without
view bounds retain the original origin-centred view.

## Logs and the published example

Each command replaces its corresponding log when invoked. Logs include stage
messages, backend diagnostics, warnings, errors, timing and publication status.
Read warnings alongside the [model assumptions](model.md).

Each of `examples/simple-boat-native/` and `examples/simple-boat-wecsim/`
contains nine files: `boat.npz`, `boat.log`, `sea.npz`,
`sea.log`, `result.npz`, `result.png`, `result.log`, `video.mp4` and `video.log`.
Both cases have 13,001 motion samples over 130 s. The native example saves
1,301 surface/video frames at 10 fps and records its 800 by 800 output resolution
in the archive. The preserved WEC-Sim reference has 1,626 frames at 12.5 fps.
The two directories share identical inputs and keep solver outputs separate.

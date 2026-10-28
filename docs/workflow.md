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
uv run topwater generate-sea DIRECTORY --record RECORD_ID --seed 42
```

This replaces `sea.npz` and `sea.log`. Omitting the seed draws and prints one.
Use `topwater.sea_states.catalog()` to inspect the packaged observed records.
Save the generated arrays for durable replay across generator changes.
Generation retains observed band energy and directional information; phases are
synthetic, so it reproduces a spectrum rather than the measured elevation trace.
The original carrier arrays are not rewritten by solve.

## Solve motion and water surface

```powershell
uv run topwater solve DIRECTORY DWELL --rest 10 --ramp 5
```

DWELL is the positive full-strength interval in seconds. Total time is
`2*rest + 2*ramp + dwell`. Ramp must be positive, rest nonnegative, and all three
durations must align to 0.01 s. Motion/packet sampling is 0.01 s; surface cadence
targets 0.08 s and adjusts slightly when needed to include both endpoints.
The published case uses exactly 0.08 s. The end time is fixed; the solver does
not extend it to force settling.

Sea phases refer directly to the dwell midpoint. Each carrier's target is
multiplied by a quintic rest/ramp/dwell/ramp/rest envelope. Fourier sidebands use
their own dispersion and retain their carrier's direction. See [the model](model.md).

WEC-Sim is the sole backend. `--solver=wecsim` remains accepted for compatibility;
the unavailable `native` selection has been removed.

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
The video uses a fixed origin-centred square view, metre graticules, telemetry
and a water-elevation gauge. It has no sea-spectrum inset.

## Logs and the published example

Each command replaces its corresponding log when invoked. Logs include stage
messages, backend diagnostics, warnings, errors, timing and publication status.
Read warnings alongside the [model assumptions](model.md).

The simple-boat example contains exactly nine files: `boat.npz`, `boat.log`, `sea.npz`,
`sea.log`, `result.npz`, `result.png`, `result.log`, `video.mp4` and `video.log`.
Its 130 s case has 13,001 motion samples and 1,626 surface/video frames.
The published video is 800 by 800 pixels at 12.5 fps, with 130.08 s encoded duration.

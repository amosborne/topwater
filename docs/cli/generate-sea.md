# topwater generate-sea

Generate a replayable directional wave realization from one packaged observed
spectrum. This creates incoming wave components, not a boat response or a
measured buoy elevation trace: the observed spectrum constrains energy and
directional moments, while individual phases are synthetic.

```powershell
uv run topwater generate-sea [DIRECTORY] --record RECORD_ID [--seed N]
uv run topwater generate-sea examples/simple-boat --record ndbc-46086-20231223T0040Z --seed 42
```

## Inputs and outputs

DIRECTORY defaults to the current directory and is created when needed.
`--record` selects a packaged observation. `--seed` is an optional unsigned
32-bit integer; omission draws a seed, printed with the source ID on success.
No boat file or network retrieval is required. This command has no dwell/ramp
arguments: it preserves original carriers as the reusable input. Packet timing
belongs to `solve`, which constructs its duration-dependent Fourier components
without rewriting `sea.npz`.

The files written are DIRECTORY/sea.npz and DIRECTORY/sea.log. The NPZ contains explicit amplitudes,
angular frequencies, propagation directions and original carrier phases; see
[the NPZ contract](../npz-format.md). Successful generation atomically replaces
an old sea.npz. Failures preserve it. Existing result/video files are untouched;
the user manages their relationship to changed inputs.

For solving, `phases_rad` specifies each target carrier phase **directly at the
dwell midpoint**, at the solver's wave reference point. The physical convention
is `cos(k dot x - omega*(t-midpoint) + phase)`, with propagation directions
unchanged. A zero phase means a target carrier crest there. The finite Fourier
approximation can depart from the
target; the original phase numbers are not shifted or normalized in this archive.
An internal Fourier sign conversion preserves this physical convention and does
not move its midpoint reference.
For example, `topwater solve DIRECTORY 100 --rest 10 --ramp 5 --solver=wecsim`
requests100s dwell and130s total, with midpoint65s. Reusing the same sea for
another dwell keeps the carrier arrays unchanged and changes packet construction
inside solve. See [solve](solve.md) and [packet search](../wave-packet-search.md).

`sea.log` records the latest invocation, selected record/seed, stages, errors and
publication status. It is replaced even when generation fails; existing data
remain protected. See [the common log format](logging.md).

## Assumptions and hard-coded choices

- Eight carriers per observed frequency band; PCG64 random generator.
- Stratified within-band frequencies and positive-MEM directional reconstruction
  on a 0.5-degree nautical-FROM grid, converted to east/CCW propagation radians.
- Each measured band's variance is retained; no amplitude rescaling, response
  fitting, pruning or boat-specific support filtering.
- Fixed deep-water conventions and gravity9.81m/s2. Timed quintic/Fourier
  packet construction belongs to solve, not this original carrier file.
- The original nine-observation catalog is preserved; the accepted stronger
  observation46086/2023-12-23T00:40Z is additionally packaged for this example.

## Fatal errors and detection

Argument parsing rejects missing record IDs or malformed seed arguments. Before
writing, generation rejects unknown observations, out-of-range seeds, malformed
observed bands/directional moments and invalid component arrays. Publication
fails on filesystem errors. Fatal errors return nonzero and identify the stage.

## Warnings and detection limits

Generation does not assess a particular boat's linear range, immersion,
flooding, wave breaking, numerical resolution or physical packet startup.
Those checks belong to solve. Source coverage is limited to measured frequency
bands; this command cannot infer unmeasured spectral tails or observed phases.
It emits no scenario-specific physical warning without a boat simulation.

# Observed sea states

`topwater.sea_states` supplies nine actual Southern California observations from
NOAA NDBC station 46086 (San Clemente Basin), with a canonical baseline. Each
record keeps its complete measured frequency-band energies and paired first two
directional harmonics. The catalog has no runtime dependency on marine-gyro.

The baseline is **2025-07-05 20:40 UTC**: measured `Hm0 = 1.3923 m`,
`Tm02 = 6.9405 s`. Its peak spectral-density band's centre corresponds to
14.8148 s; this is a band descriptor, not the period of a single wave.

The 42,575 complete 2021–2025 observed hours define Hm0/Tm02 marginal quantiles.
The baseline is the actual joint observation nearest both medians, using squared
distances normalized by each quantity's observed interquartile range; earlier
UTC timestamps break ties. Eight additional examples target the remaining
combinations of 10th/50th/90th percentiles. Already selected hours are excluded.
The nine cases are coverage examples, not nine equally probable sea states.

```python
from topwater import sea_states

state = sea_states.observed()  # canonical SoCal observation
waves = sea_states.realize(state, beam_m=9*.3048, gravity=9.81, seed=42)
saved = sea_states.dumps_realization(waves)
replayed = sea_states.loads_realization(saved, beam_m=9*.3048)
```

For the checkout's canonical boat, pass `tests.reference_boat.BEAM` and
`tests.reference_boat.GRAVITY` explicitly. These are wave inputs; combining them
with a mean-water equilibrium does not initialize a consistent nonlinear wave
field or establish a safe floating pose.

Omit `seed` for a fresh realization. The actual generated uint32 seed is retained
in `waves.seed`. A supplied seed repeats frequency, directional and phase draws.
Save the explicit arrays for durable replay across generator/version changes.

Each band amplitude is a cosine **crest amplitude in metres**, with variance
`a²/2`. Eight stratified frequencies per band share that exact variance. A
positive maximum-entropy reconstruction of four measured directional moments
supplies spreading; finite directional samples approximate that reconstruction.
Any uniform mixing used to regularize the moments is recorded per band.
Directions in `WaveComponents` are propagation angles counterclockwise from
world east; phases use `cos(k·x - omega*t + phase)`.

Phases are synthetic independent uniform draws from `[0, 2*pi)`. This recreates
the observed **spectrum**, not the recorded buoy elevation history. Finite-run
RMS varies with interference and phases; it is never normalized to force a
particular RMS or wave height.

The catalog covers **0.03–0.495 Hz**; the documented 0.020-Hz noise band is
excluded. Unmeasured tails remain unknown. Its deep-water generator rejects
any positive-energy observed support above 0.6 Hz or with wavelength shorter
than twice the supplied boat beam; equality is allowed. It neither crops nor
renormalizes an incompatible record. A 9 ft beam admits these observed bands;
an 11 ft beam can fail. Passing this incident-band screen does not establish
validity against breaking, slamming, flooding or loss of contact, and does not
qualify a full nonlinear simulation or its runtime.
Finite depth shortens wavelength at the same frequency; a finite-depth solver
must recheck its own dispersion relation. This generator does not certify it.

The packaged JSON includes NOAA URLs, frozen source hashes, exported-library
hash, observed UTC identities and selection targets. All selected energies and
paired directional moments were checked against their frozen native NOAA rows.
Directional moments retain the source export's int16/32767 precision.

To rebuild from the verified current marine-gyro export, without network access:

```powershell
uv run --locked python scripts/import_socal_seas.py `
  --library C:/Users/alex/marine-gyro/artifacts/sea-direction-check/sea_states.json `
  --source-manifest C:/Users/alex/marine-gyro/references/sea-states/config.json `
  --output src/topwater/data/socal_sea_states.json `
  --native-control tests/data/socal_baseline_native.json
```

The builder reads marine-gyro files and writes only its explicit output paths.
The stale `artifacts/report-build` library lacks directional data and is rejected.

# Existing sea-state generator: integration review and neutral export

Bounded assignment from coordinator, 2026-10-24: inspect the existing generator
for the user's higher-energy 5-10 degree roll case and dominant-direction box
alignment. The initial read-only review performed no generator call, test or
solver job. Catalog descriptors below are arithmetic on the frozen JSON, not
predicted boat response. A later authorized export is recorded separately below.
Topwater solver implementation remains held; its source was not modified.

## API and energy

`src/topwater/sea_states.py` owns `catalog()`, `observed(case_id=None)`,
`realize(case=None, *, beam_m, seed=None, components_per_band=8, gravity=9.81)`,
`dumps_realization(waves)`, `loads_realization(text, *, beam_m=None)` and
`validate_envelope(waves, *, beam_m, gravity=9.81)`.

The generator already supports higher-energy actual observations. Select a
different whole catalog record, then call `realize(record_id, ...)`; there is
no amplitude ceiling, Hs multiplier or energy normalization. Changing seed or
components-per-band changes the realization, not its measured band energies.
The API has no independent Hs/Tp energy-selection knob. A custom
`ObservedSeaState` is accepted, but labeling fabricated/scaled data as observed
would be incorrect.

Frozen source observations (Hm0 = 4 sqrt(sum(a_band^2/2)); Tm02 uses flat spectral
density inside each measured band):

| record ID | Hm0 (m) | Tm02 (s) | peak density-band period (s) |
|---|---:|---:|---:|
| ndbc-46086-20250705T2040Z (default) | 1.39226434 | 6.94052537 | 14.81481481 |
| ndbc-46086-20220424T0740Z | 2.27842051 | 6.94524674 | 8.51063830 |
| ndbc-46086-20250307T0540Z | 2.27613708 | 5.49922144 | 7.27272727 |
| ndbc-46086-20250105T1140Z | 2.27033037 | 9.29644097 | 12.90322581 |

The three high-energy records target the observed Hm0 90th percentile paired
with Tm02 50th/10th/90th percentile targets. These are joint observed records,
not freely combined independent parameters or population weights. They are not
certified to produce 5-10 degree roll, linear response, or safe boat operation.

## What must survive neutral export

- Entire measured spectrum: 48 bands covering 0.03-0.495 Hz. Per-band cosine
  amplitudes are divided by sqrt(component count), preserving each band's
  variance exactly. No active band may be cropped or renormalized.
- Each band's four measured nautical FROM moments remain paired with its
  energy. Positive MEM reconstruction uses a 720-bearing (0.5 degree) grid;
  uniform-mixing regularization scales are recorded. Finite directional samples
  approximate the reconstructed PDF, not uniquely observed wave directions.
- PCG64 uint32 seed controls stratified within-band frequencies, directional
  quantile samples and independent uniform phases. `seed=None` generates and
  records a new seed. Save explicit arrays as well as seed/provenance.
- Deep-water dispersion, gravity and actual hull beam. The incident support
  rule is f <= min(0.6, sqrt(g/(4*pi*beam))). For the 9 ft beam, the complete
  catalog support is admissible; this is not a nonlinear validity certificate.

The smallest existing neutral format is `dumps_realization`: four explicit
equal-length arrays `amplitudes` [m], `angular_frequencies` [rad/s], `directions`
[rad], `phases` [rad], plus seed and string provenance. Replaying checks exact
measured variance closure and current support rules. Preserve that immutable
JSON; consumer-specific spectral tables are derived artifacts with a proved
units/phase/directional mapping.

## Dominant direction and box orientation

Source moments a1,b1,a2,b2 describe nautical FROM bearings beta clockwise from
true north. Generator output is propagation theta = (3*pi/2-beta) mod 2*pi,
counterclockwise from world east, with cos(k dot x - omega*t + phase).
`WaveComponents` and `InitialState` document this explicitly.

There is **no existing dominant-direction field or API**. The existing
`tests/reference_case.py` deliberately starts the bow east; it does not align
the boat to its generated sea. A case adapter must name the dominant statistic.
One defensible source-based statistic is the energy-weighted first harmonic:
E_b=a_b^2/2, A=sum(E_b*a1_b), B=sum(E_b*b1_b),
theta_mean_prop=atan2(-A,-B). Its resultant sqrt(A^2+B^2)/sum(E_b) must be
reported so broad/bimodal ambiguity is visible. This is a mean propagation
direction, not necessarily the mode of the integrated MEM distribution.

The coordinator subsequently fixed the convention: default box forward yaw psi
points toward the arriving source, realized mean propagation + pi. Explicit
relative heading offsets are internal selection alternatives, not a replacement
for the generated directional spread. Both measured and actually realized first
moments/resultants are stored. Do not silently replace the spread with a beam sea.

Recommended frame adapter: retain theta_global in immutable generator JSON,
store psi explicitly, and use theta_boat = wrap(theta_global-psi) in a
boat-aligned WEC-Sim frame. This keeps the original full hull/inertia and makes
native roll about the initial box longitudinal axis. Store Rz(psi), both frame
definitions, and consistent output/render transforms. Rotating only the drawn
box while leaving physical waves/inertia unchanged is invalid. Rotating the
hull in geographic world axes instead also requires tensor/DOF transformations
and a careful definition of longitudinal roll at nonzero mean yaw.

Current native data at only 0/45/90 degrees does not cover a faithful MEM
realization's complete relative angular spread. New producer coverage or an
independently reviewed symmetry/interpolation policy is needed; do not discard
opposite-direction tail components.

## Startup and next bounded work

`realize` creates a stationary random-phase harmonic field. It does not create
an upstream packet or consistent calm-start nonlinear fluid state. Documentation
explicitly warns that combining these waves with calm-water equilibrium is not
an initialized simulation. A startup adapter/initialization policy therefore
needs review; silently applying a global ramp, rephasing into the old cos4
packet, or normalizing a finite time trace would change the requested sea.

Next proposal: select genuine higher-energy observation(s), freeze record/seed/
component count and explicit arrays, define dominant statistic and forward
convention, then derive faithful consumer inputs and new native heading
coverage. A measured response/nonlinear check is still required for the requested
roll target; generator energy alone cannot qualify it. No generator source
change is presently needed merely to obtain the catalog's higher-energy seas.

## Source evidence

- `src/topwater/sea_states.py`: dataclass and descriptors lines23-79; catalog
  lines83-114; support rules117-141; MEM144-176; realization179-227;
  serialization230 onward.
- `src/topwater/inputs.py`: WaveComponents157 onward; InitialState197 onward.
- `docs/sea-states.md`: observed records, conventions and initialization limits.
- `src/topwater/data/socal_sea_states.json`: selected joint targets and full
  observed data; no observations were generated or altered in this review.
- `tests/test_sea_states.py`: read-only inspection of native NOAA, energy,
  harmonic, seed/replay and support checks; tests were not run.

Suggested current-state update: record this generator-only read-only review and
the source/API facts while retaining the held status of all Topwater solver
implementation, tests and numerical campaigns.

## Authorized neutral-export proposal (not executed)

The coordinator subsequently authorized preparing a neutral data-export producer,
without modifying the generator or running it. Files:
`reference/sea_states/export_generator_candidates.py` and `generator-candidates.json`.
Three high-energy catalog records above, seed42, eight components/band; heading
offsets0/45/90degrees are internal cheap candidate metadata. ONE final canonical
realization/heading must later be selected; no nine-case simulation campaign is
proposed or authorized by this export.

The export stores unchanged global realization JSON plus derived boat-frame wave
arrays, rigid yaw/translation transform and frame-phase closure diagnostics. The
default heading is realized energy-weighted mean propagation + pi. A pure yaw
about the known mean CG gives x_global=R*x_boat+t, t=CG-R*CG; body phase is
phase_global+k_global dot t. Current canonical CGxy=0 makes that phase correction
zero. Exact explicit-array replay, each observed band's variance, whole support,
source invariance and independent frame elevation equivalence are checked.
Versioned source/data/config/physical input snapshots, Python/NumPy/SciPy versions,
Python executable hash and actual NumPy core/random binary hashes are recorded.
The output explicitly remains a stationary realization, with startup and dynamic
qualification unresolved.

Proposed command from repository root:

```powershell
.venv/Scripts/python.exe -u reference/sea_states/export_generator_candidates.py --config reference/sea_states/generator-candidates.json --output artifacts/reference/generated-sea-candidates-20261024
```

Use a fresh output directory. Expected: three generator realizations, each384
components, source-band/replay/frame checks pass; zero BEM/dynamics/Topwater
solver calls. One compute thread, expected under5s, budget60s. Coordinator must
allocate the sole execution owner and lease before launch. AST and config JSON
parse checks passed; producer has not been executed.

## Executed neutral export

Coordinator allocated the single command above after all previous artifact jobs
were verified terminal. Proposal SHA256
`2f61b0683ddecca1bcc0af663f750c3a9447de58ad517343fa85a67971a56231`.
Producer SHA256
`1d9f7cca702d3d8fa43a138c7827acadb609e2486c7e982f2fd0c9219cce390e`;
config SHA256
`d369bc2ee34d984721fc9650ca0d59ca84dd311951c76add572e5ff63efaf690`.
Actual generator source hash
`cad9b4a13aee0b08281f77a7379d8adc286b98c7855cffaa96d4c80453c26a4f`;
catalog hash
`138e59a360b30ec4028f77e8d1f1c629911d05c7b37d5542c2d59b51d2628ae9`.

The one export completed in0.1564761s; each record has384 components. All observed
band variances, exact explicit-array replay, entire incident support, global/body
phase-coordinate closure, frozen source/dependency invariance and output hashes
pass. Maximum per-band variance differences are6.94e-18/6.94e-18/1.39e-17m2.
No BEM, dynamics, Topwater solver, pytest or phase/seed optimization ran.

| observation suffix | realized mean propagation (east/CCW deg) | resultant | default bow toward arrival (deg) |
|---|---:|---:|---:|
| 20220424T0740Z | 345.96313610 | 0.75906727 | 165.96313610 |
| 20250307T0540Z | 335.23017729 | 0.69098338 | 155.23017729 |
| 20250105T1140Z | 343.21754067 | 0.79015173 | 163.21754067 |

The actual global realization hashes, respectively:

- `53e48556657cdedceb16a4005315630800b0e089c5f5126ed3e15fab7ce736bb`
- `16eee33cb7a3c73a0b1651e30cb71f09fe00fa6405fa71364ebb052aee7bd431`
- `8594ce2bb7a4b08e76f9df4d30fd02f6e6dcdf5ed3c8fafae99b47d710a39c58`

Authoritative output:
`artifacts/reference/generated-sea-candidates-20261024/manifest.json`, per-record
`realization-global.json`, `summary.json` and heading-offset data. Exact executed
proposal/source archives are under `generated-sea-candidates-executed-source`;
runtime snapshots also sit under the output's `source/` directory.

Terminal evidence: exact PID1960 `WaitForExit` returned and `HasExited=true`;
subsequent PID/Python process inventories were empty. Producer launches no child
processes. Synchronous shell tool chunkabb526 exited0 with complete-true output
and empty stderr. Windows Start-Process followed by Refresh exposed a null
Python ExitCode; this is preserved, not relabeled as directly observed0. The
producer's complete manifest is written only after every check and source hash
invariance assertion. Terminal/postcheck evidence is in
`generated-sea-candidates-terminal.json` and `generated-sea-candidates-postcheck.json`.
Ledger idle was written while holding the atomic claim; removal of the lease
was the final operation. No numerical job remains allocated to this owner.

These are internal candidate inputs, not three new canonical simulations.
No actual roll, startup compatibility, wetted/righting validity or eventual
single-canonical selection is established by this generation job.

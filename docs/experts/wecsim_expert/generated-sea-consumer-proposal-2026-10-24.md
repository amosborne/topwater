# REF-002 faithful generated-sea consumer proposal

Status: source/data inspection only, 2026-10-24. No new hydrodynamic or motion
job launched. This is a proposal requiring coordinator adoption and independent
review before implementation/physical reliance. REF-001 inputs/results unchanged.
The requested deliverable is ONE selected case, not the internal candidate bank.

## Source identities and known limits

- Maintained WEC-Sim checkout inspected at
  `f962103a348d63101ce669c2611f7d7c109755c3` (selected v7.1.0 based fork).
- `artifacts/reference/generated-sea-candidates-20261024/manifest.json` records
  three observed records, seed42, 384 explicit components each, complete MEM
  spreading, Hm0 approximately2.27m, and no amplitude multiplier. Full measured
  band energy is retained. Finite realization is not the measured buoy time trace.
- Every bank frequency lies in .0302418613 to .4948102898Hz. Full positive-energy
  direction support matters; REF-001 0/45/90deg excitation is inadequate.
- Previous nonlinear-adapter/weak-nonlinear comparison is superseded. This plan
  uses fixed-mean linear excitation/restoring/radiation, with an A PRIORI
  engineering operating-envelope assumption from Capytaine owner's geometry
  curves. It does not establish nonlinear dynamic accuracy or guarantee the
  requested5–10deg roll is achievable inside that envelope.

## Proposed hydrodynamic interface

Prefer direct native excitation at every actual component pair
`(omega_j, wrap(theta_global_j-yaw_initial))`. Export dimensional complex
`H_j = diffraction_force + Froude_Krylov_force`, six named modes about the same
mean CG, and an explicit component-index-to-solve map. Keep zero-amplitude array
entries in the bank; algebraically zero entries need not require a BEM solve.
Any deduplication is exact pair equality only, with provenance back to all entries.
No directional binning, narrowband replacement, spectrum fitting or interpolation
is needed for excitation. There are at most384 excitation solves per selected
heading/record, plus independent radiation work; the producer owns exact cost.
An irregular pair table is a supplemental native producer artifact, rather than
pretending it is a rectangular BEMIO frequency-by-heading grid.

Use official BEMIO for the heading-independent A/B/radiation/Ainf/H5 path, and
exact native MAT metadata for full M, CG/CB, volume and K as in REF-001. The
existing radiation grid0.1:0.1:12rad/s and native Ainf may be reusable for the
identical hull/mean frame; its finite bandwidth/memory remains a documented
approximation, not newly proven convergence. Proposed first dynamics keeps
the previously exercised10s CIC memory anddt.02 pending explicit coordinator
decision. No automatic longer-kernel/frequency refinement campaign is implied.

Capytaine owner agreed the pair-table approach and fixed native-frame physics
are appropriate in direct correspondence on2026-10-24; no new producer job or
coverage contract has yet been allocated or frozen.

## Exact forcing and case-only application

For generator convention `cos(k.x-omega*t+phi)` and native Capytaine
`exp(-i*omega*t)`, prepare all six histories as

`F(t) = real(sum_j a_j H_j exp(i*phi_local_j-i*omega_j*t))`.

Use the exported heading file's exact local phase, including its translation
correction `phi_local = phi_global + k_global.translation`. Do not regenerate
phases or apply BEMIO's imaginary sign conversion a second time. Independently
check a unit real/imaginary coefficient and a single-component analytic history
before relying on the whole-bank total. Direct harmonic forcing preserves the
entire stationary bank and avoids finite excitation-IRF truncation. This is
native BEM excitation, not a fitted/external motion trajectory or Topwater load.

`bodyClass.hydroForce` has private SetAccess in this fork. Do not overwrite it
from a case script. Proposed case-only configuration, after official cold-start
preprocessing and native body construction: use proper Simulink API to point
the ONE selected imported-excitation FromWorkspace block at a named case-local
`[time,F]` table. Keep the native parent sum, six force channels, native body,
mass adjustment, radiation and restoring intact. Retain one `elevationImport`
wave object for the already-reviewed single-branch dispatch and incident data;
clearly declare its scalar elevation as visualization/initialization input, not
the source of the direct total. Archive the old/new expression and assert exactly
one active positive path. This is an explicit case forcing adapter, not a generic
fork correction. Its source, phase and actual applied-force test need independent
review before dynamics. If the coordinator requires exclusively stock excitation
preprocessing, choose the slower full-grid/384-import alternative explicitly;
do not silently switch approaches or use unreviewed direction bin normalization.

Reproduce the independent full native signed excitation bus against the prepared
table on EVERY actual output time, all six columns. Also preserve the native
preprocessing result separately; do not label that different finite-IRF result as
the applied direct harmonic input. Native mass/H5/tensor checks remain mandatory.

## Heading, longitudinal roll and global output

Run dynamics in the unchanged mean boat axes (x forward, y port, z up), with
native mean-oriented geometry and inertia unchanged. The generated heading
file maps global sea to that frame. Native Roll is then rotation about the boat's
initial longitudinal x axis, rather than global east. Passive yaw stays off;
there is no direction-following coefficient update.

For global output/rendering, use the heading file's rigid transform once:
`p_global = R_heading p_native + translation`,
`R_global = R_heading R_native`, and transform world vectors with `R_heading`.
The generator retains its original global sea; render that exact field in global
coordinates, or use the validated local field in native coordinates. Label both
representations. Do not yaw the Solid/inertia again in the physics model or
conflate a viewport rotation with a changed wave heading.

## Proposed initialization, duration and candidate selection

Propose ONE final120s record at native t=0:.02:120, ode4, rampTime=0,
zero perturbation displacement/velocity and zero radiation memory at t=0.
The stationary sea is already present at t=0: force/elevation can be nonzero.
This is an explicitly forced-from-rest linear problem, NOT calm startup or a
steady-state fluid/body initialization. Preserve the entire transient; no prerun,
hidden force ramp, phase shift, clipping, or discarded initial trajectory. Initial
acceleration follows the stock native mass-feedback discretization documented
for REF-001. Future comparisons must use these same declared conditions.

Before any first dynamic candidate, adopt the geometry envelope and freeze a
prospective finite selection policy. Suggested cheap screen after producer
coverage: linear frequency-response/analytic forcing diagnostics across the
three immutable records and0/45/90deg heading options, reporting all six response
estimates, restoring conditioning and zero-frequency/free-mode caveats. A
frequency-domain estimate does not initialize the eventual time-domain run or
replace actual native dynamics. Do not infer bounded surge/sway/yaw from a
stationary RAO where unrestrained zero-frequency modes make it inappropriate.

Choose a first actual record/heading only after geometry/coverage evidence,
targeting5–10deg longitudinal roll while reporting other DOFs in metres/degrees
and relative to boat dimensions. "Appreciable" needs the coordinator's explicit
threshold or documented assessment; do not secretly optimize phase or force.
Preserve failed candidates, and deliver only the final frozen case. Report if
no unchanged bank meets the motion target inside the assumed operating envelope.
Assess immersion/righting relative to local water elevation and slope; global
wave-following heave alone does not prove invalid immersion. A large incident
sea on this small hull may fail the geometry screen: no result is promised.

## Coordinator feedback and revised prospective controls

The120s/10s-memory suggestion above is NOT adopted. The coordinator accepts
direct-pair H and an explicit case-only force table as a direction, conditional
on native normalization/phase/frame and all-six applied-bus checks. Initial
relative immersion must be screened before selecting a pose.

Concrete revised proposal, still unrun/unapproved:

- First record60s, ode4dt.02, CICdt.02, complete3001sample trajectory retained.
  Sixty seconds covers multiple typical6–9s periods and almost two periods at
  the lowest generated frequency; it is a finite realization, not statistical
  convergence or a steady-state claim. No automatic extension/discard interval.
- Radiation omega=.05:.05:12rad/s, exact native Ainf, causal memory40s and
  kernel spacing.02s. This is240frequencies, six radiation modes each before
  any exact reuse. Frequency recurrence2pi/domega=125.66s exceeds twice40s;
  that avoids the obvious recurrence collision, but does not prove cutoff/tail
  accuracy. Producer must agree cost/metadata and independent reviewer must
  assess the actual kernel and reconstruction diagnostics before adoption.
  Exact existing even-grid radiation values may be reused only with identical
  physical mesh/frame/solver/dependency identities; no assumption of reuse.
- Report kernel tail, and reconstruct finite-kernel radiation impedance at all
  actual incident frequencies against direct native A/B. Also report first/last
  frequency endpoint behavior. These diagnose the proposed model; they are not
  an invented coefficient tolerance or automatic refinement authorization.
- Direct component excitation uses each actual pair, with no lower-frequency
  extrapolation; the bank minimumomega is about.190rad/s. Radiation's separate
  low/high frequency integration band is not a physical incident-spectrum crop.

Proposed canonical initialization is a **held-body release** at generated t=0,
with all six initial velocities zero and radiation memory zero. Prior to release
the body is idealized as held motionless at its selected pose while the unmodified
stationary incident realization exists. Thus zero radiation memory has an explicit
model history. Holding forces are absent from the simulated t>=0 model; no wave
or force ramp is imposed. This is a declared initial-value experiment, not a
previously free-floating stationary state or calm startup. The stationary linear
diffraction/excitation field is represented by the native H coefficients; no
initialized nonlinear fluid field is claimed.

Determine the initial position/orientation prospectively for EACH unchanged
record/heading as follows, before any dynamic selection:

1. Evaluate exact full-bank eta, both horizontal gradients, eta_t and gradient_t
   at the meanCG at t=0, and eta on the full hull footprint. Record extrema and
   departure from the tangent plane. No phase/time-origin fitting is allowed.
2. Use the local tangent plane through etaCG as a geometric reference. Obtain
   the mass-balanced free-heave/trim pose with zero heel RELATIVE to that plane,
   using the already-reviewed geometry equilibrium operation in a rigidly
   transformed frame. This is a proposed bounded geometry evaluation for Cap
   owner, not a new nonlinear trajectory model. Transform the resulting body
   pose back into the unchanged native mean axes; initial surge/sway/yaw
   perturbations remain zero. Cross-check exact gravity direction and full
   world restoring torque in the tilted plane; do not blindly rotate calm-water
   equilibrium and assert world force balance. A tilted local plane is not a
   static ocean equilibrium under vertical gravity: report the required initial
   holding force/torque rather than calling the released pose globally balanced.
   If no suitable geometrically mass-balanced pose exists, or curved
   actual eta gives opening immersion/invalid wet set, reject the initial case.
3. Evaluate that pose against BOTH the exact generated curved surface at t=0
   and the adopted relative-immersion/righting envelope. A tangent-plane result
   alone is insufficient if curvature over the hull is material. A geometry-only
   diagnostic can clip against the sampled surface with explicitly bounded
   interpolation error, or conservatively bracket the surface by planes using
   the measured residual. No unsupported finite-wave force model is inferred.
4. Compute the6x6 linear harmonic response diagnostic using native H and
   `Z(omega)=K-omega^2*(M+A(omega))-i*omega*B(omega)` for exp(-iomega*t).
   Preserve full coupling/nonreciprocity; do not add viscous damping. Compare
   predicted heave/heel/trim RELATIVE to exact local eta/slopes to the adopted
   geometry range, and report all six finite-time response estimates. Treat
   near-singular Z explicitly rather than using an unreported pseudoinverse.

The release pose is not selected to maximize roll or cancel a velocity; it follows
the same geometric rule for every candidate. Initial free-mode momentum is zero.
Surge/sway/yaw have no hydrostatic restoring, so release under a finite harmonic
force generally introduces a constant-velocity homogeneous component and drift.
Steady-state RAO amplitudes alone cannot rank finite-time admissibility or predict
drift. Coordinator scope clarification forbids another custom time-domain
integrator/state-space model. The cheap screen therefore ranks ONLY harmonic
response and explicitly leaves startup/free-mode drift unresolved. Select ONE
candidate, then use the actual native WEC-Sim trajectory to assess the complete
declared initial-value problem; preserve/report rejection and seek a bounded
next data-derived candidate if needed. No parallel suite of native trajectories.
Report accumulated translational/angular momentum and integrated drift over
the full60s, with Coriolis/rotational terms distinguished when interpreting the
native nonlinear rigid-body coordinates. Do not zero detrend, recenter, constrain
or damp those modes. If fixed-mean excitation or geometry becomes inapplicable
because of that drift/yaw, reject/report the candidate; changing to passive-yaw
or a translating-body model would require a separate scoped decision.

This initial-state rule and conservative geometry check are proposed, not executed.
If the held-release experiment is unsuitable for the future comparison, the
alternative is a fully initialized harmonic particular solution INCLUDING native
radiation history and mass-feedback state. That requires a separately reviewed
state-initialization implementation; setting only RAO q/v with zero memory is
not a stationary initialization. No such extension is included here.

### Joint paired-export detail

The selected official `readCAPYTAINE.m` requires `wave_direction`,
`Froude_Krylov_force` and `diffraction_force`, and actually reads both force
variables. Therefore provide real heading0 FK/diffraction on the240 radiation
frequencies in the canonical BEMIO netCDF. Do not fabricate zeros or make a
radiation-only file appear to contain physical excitation. That heading0 H is
not the generated-bank forcing.

Proposed separate per-record paired MAT schema, agreed direction with producer:
`component_index`1:384; column vectors `omega_rad_s`, `heading_native_rad`,
`amplitude_m`, `phase_local_rad`; complex384x6 `fk_complex`,
`diffraction_complex`, `excitation_complex`; `dof_names`, `cg`, `rho`, `g`,
`globalbank_sha256`, `headingfile_sha256`, `yaw_global_rad`, and
`global_from_native_R`/`global_from_native_t`. Excitation is dimensional per
unit surface elevation, N/m or Nm/m, with native exp(-iomega*t). Preserve every
component index and exact input array, including zero amplitudes. No additional
BEMIO imaginary conversion is applied to these direct native complex arrays.
Capytaine owner controls exporter implementation and final schema/version/hash.
Current proposal covers only three DEFAULT offset0 incoming-aligned banks;
internal45/90deg alternatives are excluded from the new producer job.

## Execution boundary and stopping point

Current stopping point is this proposal and joint producer discussion. No
implementation, numerical screening, BEM/IRF processing, motion run, library
change or Topwater solver work is authorized by this note. Next coordinator
decision: adopt/modify direct-pair forcing, initialization/duration and selection
policy after geometry results; then allocate producer coverage and a bounded
source/regression implementation with independent reviewer and serial lease.

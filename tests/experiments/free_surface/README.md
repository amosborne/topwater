# Free-surface solver validation experiments

`hos_graph.graph_field` is a small CPU reference for the periodic deep-water
nonlinear-graph Dirichlet-to-velocity component. It is not a boat solver or an
open-water near/far coupling implementation. For flat-plane harmonic traces
`f_m`, it uses `f_1=psi` and
`f_m=-sum(eta**j * |D|**j f_(m-j)/j!)`, then Taylor-extends Cartesian velocity
to the current graph through the chosen total HOS order. This follows the
surface Taylor construction in [Shi, Zhu and Xu (2023)](https://doi.org/10.1016/j.joes.2023.07.003).
Fourier modes describe space; this is an instantaneous elliptic component for
time-domain evolution, not a frequency-domain body solver.

Run `uv run --locked --python 3.12 pytest tests/test_hos_graph.py -q -s`.
The manufactured 96 by 80 m periodic patch has two oblique exact harmonic
potentials and a curved graph with nominal slope scale 0.12. The independent
checks compare all three velocity components, upward unnormalised normal flux,
and boundary trace obtained by direct exponential Fourier continuation (not
the cancelling Taylor recursion). Order, grid and product-padding controls
must be read together: padding cannot repair underresolved input and an
unresolved unpadded calculation can look better through alias cancellation.
The chosen patch scale/slope exercises moderate ocean-wave geometry; it is
not the observed sea realization, a force comparison, or trajectory evidence.
Higher HOS order is not guaranteed to converge for an arbitrary steep graph.

## Time-domain solver and validation

Topwater's intended solver architecture is time-domain: advance the nonlinear
fluid and all six body degrees of freedom. Capytaine belongs in the validation
tools. Its frequency-dependent radiation coefficients must be interpreted as
time-dependent reference forces, including instantaneous added mass and causal
fluid memory. This reference must not replace Topwater's nonlinear fluid solve.
The three existing frequencies alone do not establish a broadband transient
reference; frequency sampling, tails and memory truncation require checks.
`tests.capytaine_time_reference` now supplies a test-only finite-band kernel,
causal six-mode force convolution and coefficient back-check. Four analytic
coupled-memory tests pass, including startup, causality and time refinement.
No suitable broadband canoe dataset has yet been supplied to this adapter;
the analytic checks do not establish Capytaine correlation.

The P4 time candidate is `p4_canoe_release`. Its 1.6 s finite-pose,
six-velocity release was stopped after about 46 minutes for approximately
1.65 ms of inferred accepted motion. The Windows interruption did not save
accepted trajectory arrays. This runtime is unusable; another long integration
must wait for a measured reduction in stage cost and an explanation of the
small timesteps. The completed P4 trajectory remains the earlier 1 ms smoke case. The
prescribed-heave time runner and matched closed-tank linear control are ready
but not yet exercised. The closed tank still needs a validated time-domain
outgoing boundary before open-water comparison.

### Measured P4 stage cost

`profile_p4_canoe_stage` executes one unchanged `prepare/evaluate` pair at the
saved flat-surface, finite-pose, six-velocity canoe state. On this machine,
with one BLAS thread, the cold cProfile diagnostic measured 19.67 s of geometry
startup, then 63.53 s of preparation and 52.77 s of nonlinear evaluation.
Both the physical state and complete canonical RHS matched the saved values
bit for bit. The recurring 116.29 s divides into:

| Work | Wall seconds |
| --- | ---: |
| Current fluid stiffness assembly | 23.60 |
| Current fluid and surface sparse factorizations | 37.46 |
| Remaining preparation | 2.46 |
| Full nonlinear force/rate evaluation | 52.77 |

Cell integrals account for 49.60 s within force evaluation. Cell metrics are
computed twice per cell, once during assembly and once during force evaluation;
their 49.54 s cumulative time overlaps those rows and must not be added again.
The separate startup geometry factors cost 12.50 s, included in startup only.
Raw evidence is `artifacts/p4-calm-rhs-profile-20261022.{json,npz,prof}`;
these are local artifacts. This is a profiled single-state cost diagnostic,
including instrumentation overhead, not qualified trajectory throughput.

The complete simulation must cost at most 20 wall seconds per simulated second;
real time is preferred. With three fresh RHS evaluations per accepted RK23
step, even zero rejections would require `3 * RHS_seconds / step_seconds <= 20`.
This candidate needs substantial reductions in recurring work as well as an
independent timestep diagnosis. Neither a factorization-only speedup nor a
tolerance change alone is justified by this measurement.

The first retained optimization consolidates the physical metrics into batched
small matrix products and evaluates kinetic energy/shape force without unused
volume/gravity terms. On the same complete saved-state RHS, preparation fell
to 50.02 s and evaluation to 11.02 s: 61.04 s recurring total, about 1.91 times
faster. Stiffness assembly took 9.38 s and kinetic cell integration 7.81 s.
The physical-state maximum difference was 1.25e-16; the complete RHS relative
L2 difference was 1.51e-13. Quadrature and geometry guards are unchanged.
Fifteen focused physical cell, shape-work and geometry-trial checks also pass.
The result remains far from the runtime ceiling; these are instrumented
single-state measurements, not trajectory throughput.

The exact saved 79,261-row/5,501,437-entry current fluid matrix also permits
solve experiments without more geometry or assembly. MMD ordering preserved
all six response fields (relative difference 1.96e-14) and all 36 added-mass
entries (maximum absolute difference 3.44e-12), but factorization still took
33.68 s. This modest ordering result does not solve repeated setup cost.
`p4_coupled_solve` instead eliminates rigid velocity into the positive rank-six
action `K + rho*B*M^-1*B.T`, avoiding construction of six auxiliary fields.
The saved calm case and all six body directions agree with the existing block
elimination to roundoff: maximum velocity difference 3.22e-14, field relative
difference 2.27e-14, fluid scaled residual below 7.50e-14. Each needed six PCG
iterations and 0.68--1.34 s with an exact-current fluid preconditioner; its
28.61 s setup was measured separately. This verifies current-matrix algebra,
not reuse after motion, nonzero-surface-potential mapping or coupled throughput.
Those remain required before adopting this solve path.

The next bounded probe assembles one actual nonlinear Euler-predicted state at
delta=1e-4 s, including nonzero surface potential (L2 norm 2.15e-5). Its fluid
matrix differs from the saved one by 0.111% in Frobenius norm. Holding only the
old solver setup fixed, a coupled Woodbury preconditioner takes four iterations
and 1.09 s to solve the **current** equations; field relative difference from
fresh direct elimination is 9.11e-15 and maximum body velocity difference is
2.78e-17. Current assembly still costs 10.36 s. Setup took roughly 24 s and is
separate from the reused solve. Evidence is
`artifacts/p4-calm-predictor-operator-reuse.{json,npz}`. This proves one nonzero-psi
mapping and lagged-setup comparison, not an operator-refresh policy or trajectory.
Approximation/caching of the physical operator remains a separate experiment.

`p4_metric_action.MetricAction` is the next exact-current candidate: retain six
weighted contravariant metric coefficients per cell/quadrature point, then
apply reference derivatives, metrics and their transpose without assembling
element/global stiffness matrices. On all 7,076 actual nonlinear canoe cells,
FP64 actions take 42--50 ms warm including host/device transfers; the maximum
recorded whole-batch sample was 57.6 ms. Actual-field, smooth and contact-supported
actions agree with independently formed CPU stiffness actions within 2.5e-13
relative. Physical curved-cell work, reciprocity and invalid-Jacobian checks pass.
The current-geometry cache occupies 175.32 MB and setup took 2.95 s, including
0.99 s of metric calculation. Evidence is `artifacts/metric-action-full.json`.
This cache currently describes one geometry; it is not a refresh controller,
coupled iterative solve or full nonlinear RHS timing. Force evaluation and
geometric updates still need integration and measurement.

The frequency-only experiments below are supporting diagnostics, not a second
Topwater solver architecture. The long quarter-PML experiment was stopped after
saving complete six-mode matrices at 2 and 4 rad/s so the actual nonlinear time
simulation takes priority. The 6 rad/s matrix and full-family gate were not run.

## Earlier P1 nonlinear validation candidate

The earlier P1 candidate is `variational_canoe`, on `experiment/canonical-canoe`.
It advances all six body coordinates and the nonlinear moving free surface
using forces differentiated from the same discrete action as the fluid solve:

```text
uv run --locked --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.variational_canoe --duration 1 --initial-scale 10
uv run --locked --group correlation --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.linear_canoe_radiation --n 24 --exterior-scale 2 --omega 2 4 6 --capytaine
uv run --locked --group correlation --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.canoe_radiation --variational --mode heave --periods 6 --capytaine
```

The older `moving_canoe` pressure formulation below remains a comparison
control: its actual canoe energy-work imbalance did not improve with spatial
refinement. Agreement with its trajectories is therefore insufficient.
The variational candidate instead matches independent fixed-momentum energy
differences, preserves constant-potential gauge invariance, and balances an
independently evaluated physical energy derivative. Its complete nonlinear
one-second n=12/wet2/initial-scale10 trajectory completes with 201 accepted
steps and 605 RHS evaluations. All six body velocities evolve. Maximum energy
change is 1.55e-5 J, or 0.00455% of the 0.341 J initial disturbance energy;
water volume changes by at most 2.70e-9 m3. Minimum rim clearance is 0.199 m.
The same full case with maximum timestep halved to 0.0025 s and tolerance
tightened to 1e-5 completes with 401 accepted steps and 1205 RHS evaluations.
Energy drift decreases by 7.97 times to 1.95e-6 J (0.000571% of the disturbance),
and maximum water-volume drift decreases to 3.37e-10 m3. At the identical final
time, position components differ by less than 4e-9 m, world rotation components
by less than 4e-8 rad, and world linear/angular velocities by less than
7e-8 m/s and 5e-7 rad/s. Reference-domain endpoint eta/psi L2 differences are
6.41e-5 and 2.62e-5 relative. The only intervening implementation change on
this free-evolution path extracts the unchanged operator assembly into a helper.
Two time levels show consistency; they do not establish convergence order or
spatial/open-water accuracy. Diagnostic coarse/fine wall times were 2478/6966 s
during concurrent validation work. This remains a major performance limitation;
its cost must be reduced before further expensive trajectory validation.

`linear_canoe_radiation` solves the exact rest linearization of that action for
all six radiation modes. It retains the full contact-node hull flux loads,
consistent weak sponge projection, phase, reciprocity and dissipation checks.
An actual nonlinear-stage amplitude-refinement test checks the fluid and all
six body acceleration first variations. This provides a complete matrix diagnostic against Capytaine; it
does not establish nonlinear force accuracy, outgoing-boundary convergence,
long-run stability or full-case performance. Current coarse tank comparisons
do not pass every coefficient gate. Optimization now takes priority so that
nonlinear six-DOF and Capytaine validation become affordable.

`canoe_radiation --variational` integrates surface heights and their conjugate
momenta while prescribing one body motion. Its controller reaction comes from
the same nonlinear action and is converted to six world-axis hydrodynamic
forces/moments. The far sponge uses consistent weak projection and includes
the geometry-dependent momentum transformation. Instantaneous acceleration,
nonlinear directional refinement and undamped gauge checks pass. Full-period
force fits, amplitude/time refinement and all-six-mode correlation remain
required. Reports label the formulation so the older pressure controls cannot
be combined with variational results in a coefficient family.

`linear_canoe_transient <radiation-report.json>` compares a saved variational
trajectory to the identical rest-linearized tank, including zero initial fluid
state, quintic motion startup, all six reaction components and the same weak
sponge. Its independent forced-mode check retains the natural startup response;
a steady harmonic response would miss it. Surface differences use the constant
rest consistent mass and trapezoidal time weights. This separates transient
integration discrepancies from finite-domain errors in the Capytaine comparison.
`release_refinement <coarse.json> <fine.json>` compares completed releases at
their identical endpoint, including world rotation and reference-domain surface
L2 errors. Two levels measure a difference and do not establish temporal order.

The variational heave case now completes six periods (9.425 s), with 474
accepted steps and 1425 RHS evaluations at amplitude 1e-4 m. Compared with its
identical zero-start linear tank transient, heave-load time RMS differs by
0.0271% of the dynamic load; mass/time L2 elevation and potential differences
are 0.0333% and 0.0203%. Identical-window fits give added mass/damping
325.629 kg / 1206.766 kg/s in the nonlinear run versus
325.644 kg / 1206.770 kg/s in the linear transient. The same tank's steady
harmonic values are 330.660 kg / 1243.511 kg/s, while matched Capytaine gives
352.732 kg / 1200.774 kg/s. Thus transient and finite-domain errors remain;
the close damping fit alone does not qualify the method. Successive-period
force changes are still 2.62% and 0.629%. Amplitude/timestep refinement and
completed all-six-mode correlation remain required.

The current P1 tetrahedral exterior also fails the declared wave-phase gate:

```text
uv run --locked python -m tests.experiments.free_surface.exterior_dispersion
```

This command deliberately exits 1 while the gate fails. One periodic column
uses the actual exterior stiffness and consistent surface mass, checked against
an independently derived stencil and the finite-depth gravity-wave dispersion
relation. At n=48 and omega=6 rad/s, axial frequency errors are 10.7% and 6.32%,
predicting 3.37 and 1.99 rad phase errors after five periods. These are spatial
errors without a hull, sponge or timestep. Separate horizontal/vertical controls
show that concentrating only the vertical layers cannot resolve the failure.
The nonlinear body/work checks above remain useful, but do not qualify this
exterior propagation method. Higher-order cell candidates are being tested
before changing the coupled solver; no accuracy gate is relaxed.

`p4_pocket` tests a contained higher-order candidate against an independent
oblique harmonic field on the same flat box. At n=6 and n=12, whole-surface
flux errors for full P4 FEM decrease from 0.349% to 0.0224%. Coupling a P1 BEM
pocket to that exterior gives 25.3% and 9.65%, despite shared potential and
weak-flux balance residuals below 3e-14. The mixed error includes local P1
surface approximation and interface trace restriction; it is not a measured
wave-reflection coefficient. These static tests do not select a new nonlinear
solver or establish canoe, outgoing-wave or performance accuracy.
The additional P1-interface control retains P4 elements on both sides of
the same P1 trace restriction: flux errors are 55.5% and 19.4%. Thus accurate
individual element spaces and algebraic interface continuity are insufficient.

`p4_canoe` uses an optional temporary `meshpy` installation to tetrahedralize
the actual current canonical canoe's closed local fluid pocket, preserving
every marked boundary triangle. The bounded depth-one-metre test gives 3716
tetrahedra and 43057 P4 nodes. Its tetrahedral/boundary volumes agree within
6e-15 relative error; an independent harmonic-quartic weak-form check has
2.83e-14 relative residual. All six inertia diagonal differences from P1 BEM
on the identical physical domain lie between -5.56% and -0.556%. This is a
single-mesh comparison, with unresolved waterline singularity convergence;
neither calculation is an exact reference. Its side walls and bed differ
from Capytaine's unbounded depth-four-metre reference, so this result is not
a Capytaine correlation or a nonlinear solver selection.

`p4_whole_canoe` joins that local volume to the original depth-four-metre
exterior, sharing all P4 trace modes without a P1 restriction. Its 7076
tetrahedra have 82429 active P4 nodes; three unused input vertices are excluded.
The independent quartic-patch, linear-system and energy/work residuals are
below 3e-13. Inertia diagonal differences from the identical finite-box mixed
solver remain between -5.58% and -0.564%. No outgoing-wave or spatial
convergence follows from this static result.

The contained `p4_geometry` calculation checks prerequisites for actual
nonlinear P4 motion. Curved physical harmonic flow, independent boundary
volume/gravity integrals, fixed-potential shape work, translation and gauge
invariance pass. A rational kinetic-integral example checks quadrature
refinement for both energy and its geometry gradient. Its current curved
volume assembly recovers a known physical harmonic field by solving interior
unknowns on a six-tetrahedron curved domain. Positive Jacobians at quadrature
points do not certify positivity between those points.

`p4_canoe_contact` enriches the actual canonical hull's boundary chart. All
surface heights are independent; new waterline nodes slide within the unchanged
physical planar facets. Current-pose finite differences verify all six body
Jacobians, while off-grid checks verify facet tangency and a shared curved
waterline trace. The sampled chart guard stops excessive contact modes.
`p4_surface` integrates the physical surface one-form from those geometry
Jacobians. Extending all motion through a distance-based P4 graph inverted
thin cells even under 5 mm heave. A reference P1 FEM extension keeps the
corners valid on the same mesh and pose; `p4_volume_map` now lifts that corner
motion and extends only the genuine high-order boundary residual. The exact
transpose supplies geometry-force pullback. `p4_whole_geometry` passes the
unchanged finite-pose volume guards, closed boundary/cell moment comparisons
and six-motion derivative checks. Its fixed maps move coordinates; they do
not freeze the current fluid operator.

`p4_canoe_field` and `p4_variational_canoe` compose a current curved operator
with the full boundary velocity action. One actual stage with all six body
velocities and independently moving high contact modes has relative weak
normal-flow residual 1.12e-13 and fluid work-identity error 4.44e-16 J.
Independent fixed-potential/rate action differences reduce the tested force
direction error from 4.67e-8 to 1.15e-8 when halving the difference step.
The tracked `p4_canoe_stage` driver reproduces this finite-pose, all-six-velocity
seed, including independent high contact modes and scalar-action differences.
From a fresh checkout, generate a new seed and release it with:

```text
uv run --locked --with meshpy==2026.1.1 python -m tests.experiments.free_surface.p4_canoe_stage --output artifacts/p4-current-stage.json
uv run --locked --with meshpy==2026.1.1 python -m tests.experiments.free_surface.p4_canoe_release artifacts/p4-current-stage.npz --duration .001 --max-step .001 --output artifacts/p4-current-release.json
```

These commands temporarily add the optional MeshPy preprocessor and build the
actual mesh and current operators; they are substantial validation runs.
The stage report records package versions. The stage NPZ includes explicit reference coordinates, element
connectivity and ordered surface master IDs. Readers verify these arrays exactly;
legacy frozen seeds retain their documented deterministic-ordering assumption.
New outputs leave the frozen `p4-nonlinear-canoe-stage.*` evidence unchanged.
For new independent checks, pass `p4_stage_audit
artifacts/p4-current-stage.npz` and an explicit `--output`; run
`p4_rest_geometry_checks --output <new-rest-report.json>` before
`p4_first_variation --stage artifacts/p4-current-stage.npz --rest-checks
<new-rest-report.json> --physical-audit <new-work-audit.json> --output
<new-first-variation.json>`. These remain experimental diagnostics.

`p4_stage_audit` then rebuilds five current fluid systems to check independent
physical energy and volume along that saved full RHS. Direct velocity quadrature
uses order 10, versus order 8 in the primal, and closed-boundary gravity/volume
uses order 12. The independently computed initial excitation agrees with the
Hamiltonian within 4.61e-15 J. Halving the directional difference step reduces
energy-work error from 4.68e-7 to 1.17e-7 W and volume-rate error from
7.25e-10 to 1.81e-10 m3/s, approximately fourfold in each case. Both satisfy
the declared instantaneous diagnostic tolerances. The saved stage lacks master
IDs; rebuilding assumes deterministic ordering and verifies an exact fresh
physical-state match. This remains one stage, not a P4 trajectory, general
rim/failure coverage, spatial accuracy or Capytaine qualification.

`p4_body_action_audit` additionally checks six separate body directions at the
same finite-pose, nonflat-surface, all-six-velocity state. One current stage is
compared with 24 geometry-only scalar-action evaluations, holding nodal
potential and canonical rates fixed. Fine-step force errors are 4.40e-11,
4.29e-10 and 1.37e-4 N; angular errors are 2.52e-8, 9.79e-6 and 2.46e-10 Nm.
Vertical force, roll and pitch errors decrease about fourfold when the step is
halved. The small yaw discrepancy does not decrease; its cause is not established
by this check. Maximum fine relative error is 2.08e-6. These are generalized
action forces with explicit world angular tangents, not raw pressure torques.
This validates separate components of the discrete action derivative at the
primal quadrature; spatial accuracy and finite-motion evolution remain separate.

`p4_canoe_release` now completes a 1 ms integration smoke check from that actual
state: three accepted adaptive RK23 steps and 13 fresh nonlinear RHS calls,
including a rejected initial trial. All six velocities evolve. Maximum accepted
Hamiltonian change is 1.21e-9 J and fluid-volume change is 1.14e-12 m3.
The conservative minimum rim-clearance bound remains above 0.203 m; it uses
all P4 Bernstein height coefficients, rather than nodal maxima. A physical
polynomial test catches an interior peak missed by every nodal sample.
This very short run establishes execution through the coupled integrator,
not temporal accuracy, representative motion or qualified throughput.

For the longer physical release, `p4_canoe_stage --calm` selects the same finite
pose and six velocities with initially flat elevation and zero surface-potential
trace. The interior flow remains compatible with the moving hull; the water is
not globally motionless. The full nonlinear RHS generates the surface motion.
This separate seed defaults to `artifacts/p4-calm-stage.json/.npz` and is consumed
by the existing `p4_canoe_release` runner. The original nonflat/high-contact audit
and its saved results remain separate. The calm initial stage has completed:
the weak-normal relative residual is 1.30e-13 and independent scalar-action
difference errors fall from 4.67e-8 to 1.15e-8 when the step is halved.
The attempted 1.6 s release stopped before its first accepted step because
the initial 0.01 s proposal produced an invalid sampled cell Jacobian.

`GeometryRK23` now rejects that specific numerical trial failure and retries
with a smaller step from the unchanged accepted state. Physical contact/rim
failures and other errors still propagate. Analytic ODE controls check recovery,
backoff after an ordinary accuracy rejection, physical-error propagation and
bitwise agreement with ordinary RK23 on valid stages. Together with the cell
and offline comparison checks, 20 focused tests pass. The longer release retry
and its tighter same-case repeat remain pending. `compare_p4_releases` will
retain six-axis motion, raw energy/volume drift and exact common-endpoint
differences, requiring identical initial stage, masters and runtime sources.
If a later step fails, the release CLI preserves complete accepted snapshots
in separate `-failed.json/.npz` files and rethrows the original error. Five
focused checks cover this failure path; failed trial states are excluded.

`p4_rest_operator_bridge` compares the actual nonlinear and frequency helpers
on the same saved coarse canoe mesh without a fluid solve. Two zero-stretch
frequencies isolate stiffness and surface mass separately; their relative
differences are 4.51e-15 and 2.78e-15. The combined six-column body-load
difference is 1.61e-15. An independent harmonic-quartic boundary-flux control
has relative weak-form errors 2.35e-14/6.36e-14. This verifies the affine-rest
assembly connection only. Different fluid resolutions and the missing causal
time-domain outgoing boundary still prevent combining the closed-tank nonlinear
and open-frequency results into an end-to-end open-water qualification.

`p4_rest_geometry_checks` separately checks the actual canoe without assembling
another fluid operator. At rest, Ceta agrees with the surface mass within
2.09e-17, Cv is zero, and all six normal-flux loads agree with independent
affine integration within 1.46e-14 relatively. The full six-coordinate gravity
Hessian approaches independently clipped-hull hydrostatics: halving the
difference step reduces its largest error fourfold, to 0.00169 N m/rad in
pitch against 5478.914 N m/rad. Selected far-surface and high-contact gravity
actions agree with rho*g*M to roundoff; their tiny errors do not show an
asymptotic convergence order.

`p4_first_variation` now checks one mixed physical disturbance with all six
body components, surface heights and potentials. Four signed current-geometry
nonlinear stages approach the independently written undamped rest equations;
the canonical momentum prediction retains fluid/body coupling and the
independent hydrostatic stiffness. Halving amplitude reduces all resolved rate
and momentum errors approximately fourfold. At the smaller amplitude, linear
and angular body momentum-rate vector errors are 4.51e-8 and 8.83e-8 relative;
surface height-rate rest-mass L2 error is 5.99e-7 relative. Translation-rate
errors are already at roundoff. This is one mixed direction, not six
independent columns, prescribed-motion reaction or Capytaine qualification.

`p4_radiation` completes the n=12/extent1/sponge4 rest-candidate family at
omega=2,4,6 rad/s against the same immersed-geometry/CG/depth4 Capytaine reference.
The full-family gates fail: maximum added-mass/damping tolerance ratios are
2.31/16.78. Raising sponge quadrature from 8 to 12 changes matrix responses
by less than 3e-6 relatively and does not resolve that discrepancy. This uses
an analytic nonnegative physical sponge, differing from the old P1 nodal
interpolant, so the comparison does not isolate potential order alone.
The largest exposed LU arrays occupy 3.55 GB, with 7.21 GB cumulative process
peak working set; these static resource diagnostics are not simulation timings.

`sponge_strip` isolates the absorber in the full finite-depth elliptic problem,
including oblique waves, consistent P4 surface projection and both original
sponge equations. Modal ports retain propagating and evanescent depth modes.
Across omega=2,4,6 and incidence 0/45 degrees, fine-grid transparent controls
have maximum reflected amplitude 1.24e-8 and transmission phase error 2.19e-6
rad; the exact reflecting-wall complex-amplitude error is at most 4.40e-6.
Independent oscillatory integration checks the port projection. Four focused
tests pass, including physical transmission, reflection phase and geometry.

Two widths and two strengths were fixed before measuring the absorber, with
two spatial levels for every frequency/angle. At omega=2 rad/s (wavelength
14.48 m), fine-grid reflected wave amplitudes are:

| Absorber width (m) | Strength (1/s) | Normal incidence | 45 degrees |
| ---: | ---: | ---: | ---: |
| 1.594 | 2 | 0.314 | 0.335 |
| 1.594 | 4 | 0.543 | 0.612 |
| 7.878 | 2 | 0.206 | 0.400 |
| 7.878 | 4 | 0.327 | 0.500 |

The short, strength-four layer is the current P4 tank's transverse setting.
Across the full fixed family, spatial refinement changes complex reflection
by at most 1.13e-4, while independently integrated physical power balances
within 6.11e-15 of incoming power. Selected port truncation/integration controls
change complex reflection by 2.89e-10 / 1.69e-14. Thus large reflection is
resolved in these body-free cases: dissipation alone did not make the sponge
an accurate open boundary. This does not quantify canoe-load error, corner
effects or damping of the hull's near field; no outgoing method is qualified.

`pml_strip` tests a separate matched-coordinate candidate against an analytic
full finite-depth wave field, including both propagating and evanescent port
modes. It transforms the volume metric and surface measure together; it does
not increase the old physical sponge strength. At fixed width 7.878 m and
stretch strength 8/s, the omega=2,4,6 and 0/45-degree family has maximum fine-grid
complex-reflection error 1.03e-8. The largest reflected amplitude, 6.33e-5,
matches the predicted finite-wall return. Maximum quiet-region physical
velocity L2 error is 0.0446%, and independently integrated transformed-form
power balances within 1.40e-13. Two spatial levels and selected rational-metric
quadrature controls are retained; three focused tests pass. These are
frequency-domain propagating-wave controls only.

`pml_evanescent_strip` also checks the lowest evanescent depth mode against
its independent full analytic field. Across the same frequency family, two
transverse wavenumbers, two widths and two spatial levels, maximum fine-grid
complex-reflection error is 6.72e-9; spurious propagating amplitude is below
2.51e-9. Three focused tests pass. The near-field return decays with physical
layer length, not imaginary stretch strength: at zero transverse wavenumber,
the 1.594 m wall returns 12.8%/22.9%/26.1% amplitude at omega=2/4/6, while
7.878 m gives 0.00390%/0.0694%/0.131%. Strength-zero/four/eight controls retain
those amplitudes and change phase. An evanescent incident mode carries zero
power, so this diagnostic uses field/reflection errors without a power ratio.
Whole-canoe radiation, domain convergence and a causal nonlinear time-domain
boundary remain required before adopting this candidate.

`pml_corner` extends the coordinate coefficients to signed x/y layers and
their intersections. Known finite-depth plane waves satisfy independently
integrated volume, surface and conormal flux balances for all 35 degree-four
test polynomials in the physical pocket, one layer and all four corners at
omega=2,4,6. Two focused tests pass, including a deliberately incorrect corner
Jacobian that fails the balance. This checks the 3D equations, without a mesh
solve; discrete propagation and full six-mode canoe correlation remain open.
Compact-layer refinement and speed work are deferred while those correctness
checks and finite nonlinear motion are completed.

`canoe_wet_reference` retains the original 704-panel Capytaine hash, CG, depth
and frequency-family checks. The pinned geometry regenerates exactly without
TetGen or a fluid solve. Different numerical hull subdivisions additionally
require coverage of every original clipped planar facet: containment, winding,
area and explicit pairwise overlap checks. Two focused tests pass on the
canonical canoe, including a missing panel, reversed panel and an equal-area
hole/duplicate counterexample. Identical physical coverage does not establish
the accuracy of the field computed on a coarser mesh.

`quarter_canoe_control --reference <pinned-report.json>` now executes a coarse
canonical canoe on one quarter mesh and its exact four-copy full mesh. The
1049/3761-node comparison independently assembles all six full body loads and
the complex wave operators, including reflected volume/surface quadrature.
Full versus quarter raw six-mode response differs by 1.69e-14 relative; the
full residual is 1.32e-13 and four parity matrix/load restrictions agree within
2.1e-16. The mirrored hull covers the pinned reference's 36 wet physical parent
facets without overlap. `--mesh-only` runs the geometry and coverage checks
without a fluid factor. These deliberately coarse cells do not resolve the
outgoing-wave family. This establishes the discrete symmetry reduction only;
the wave-resolved mesh and Capytaine coefficients remain unqualified.

`quarter_canoe_preflight --envelope 16 --reference <pinned-report.json>`
now builds the actual uncompressed candidate: 110435 P4 nodes, 9482 tetrahedra
and 6485 surface nodes. It retains the tested depth levels, pocket/absorber
planes and maximum spacing envelope. The 44 quarter hull triangles cover the
same physical facets as the pinned full reference after reflection. The
conceptual full mesh has 435913 nodes; this preflight does not allocate it.
No fluid operator or factor is built. Mesh counts do not bound LU fill or
establish wave/Capytaine accuracy; a memory-bounded frequency solve is next.

`quarter_canoe_frequency` runs the first actual omega=2 heave comparison under
a queried 9 GiB process commit limit. The 110435-row scalar matrix has 7837923
nonzeros and uses 157 MB of sparse arrays. Default COLAMD factorization failed
with `MemoryError`; it produced no coefficient. The unchanged physical case is
solved successfully with `MMD_AT_PLUS_A` ordering. At omega=2 the heave added
mass is 584.239864 kg against Capytaine's 589.851139 kg (0.9513% difference),
and damping is 443.391742 against 446.360957 kg/s (0.6652%). The relative linear
residual is 5.75e-13; measured peak process commit is 9629724672 bytes, close
to the 9 GiB limit. This worker saves the exact scalar
input before factorization and rejects reuse after changes to equations, mesh,
reference, loads or data. Seven focused integrity tests pass. This is one
coefficient pair rather than the required six-mode frequency family; spatial
and outgoing-boundary stability remain unqualified.

`p4_canoe_radiation` and `p4_linear_canoe_transient` now provide the nonlinear
prescribed-heave startup and identical-rest closed-tank comparison drivers.
They retain all six load rows and the zero-start natural response. These
drivers have not yet been numerically exercised. Prescribed heave complements
the separate free six-DOF release; it does not establish a six-axis trajectory
or an outgoing-water Capytaine comparison.

The archived Capytaine reference itself passes the existing whole-family
stability gate when refined from 704 to 2816 panels at omega=2,4,6. The maximum
added-mass/damping tolerance ratios are 0.08354/0.07549. The fine panels match
ordered four-child subdivisions of the coarse panels to 2.22e-16 m, with the
same Green-function settings. `capytaine_reference_stability` checks these
existing files offline; its six integrity tests pass. This establishes
reference stability across two resolutions, not Topwater agreement.

The prescribed nonlinear controller also passes an independent physical-work
check on the moving canonical canoe with nonflat surface, all six velocities
and imposed accelerations. Without damping, its controller supplies 1.15926 W;
with the weak sponge, controller plus sponge supplies 0.674319 W. Independent
physical-energy differences at two displacement steps approach those values
with residuals below 5e-9 W and 1e-9 W respectively. This checks instantaneous
work, separately from the still-required complete radiation trajectories.

## Spatial-pocket accuracy and construction baseline

Two checkout-only commands test an actual spatial decomposition:

```text
uv run --locked python -m tests.experiments.free_surface.pocket
uv run --locked --group correlation python -m tests.experiments.free_surface.pocket_canoe
```

`pocket.py` joins a local BEM domain to a sparse tetrahedral FEM exterior by
sharing potential nodes and summing weak fluxes over the full interface.
An independent oblique harmonic field on a curved graph supplies the reference.
Moving the artificial partition by pi/3 on the unchanged physical mesh gives
whole-surface flux differences of 11.17% and 4.35% at two refinements, below
the full-FEM errors of 23.05% and 7.85%. The executable check requires decreasing
partition and mixed-minus-FEM errors below that baseline. This tests the coupling
implementation; it does not qualify the remaining physical flux error.

`pocket_canoe.py` uses the same canonical canoe, with a local domain ending 1 m
below the surface and a sparse exterior extending to the 4 m bed. It retains all
six reaction columns, a fresh 704-panel Capytaine instantaneous-inertia reference,
individual diagonal errors and setup/assembly/factor/solve timings. The refined
case differs by at most 3.4% on the six diagonal terms. Finite-box/depth effects
remain mixed with spatial error; this is not an open-water correlation gate.
Physical hull geometry stays fixed during coplanar panel refinement.

The material-layout falsifier separates subdivision before clipping from the
reference's subdivision after clipping:

```text
uv run --locked --group correlation python -m tests.experiments.free_surface.material_layout
```

All layouts use the same physical hull, fluid formulation and exterior, with
fresh Capytaine comparisons at rest and two heave/heel poses. The transverse
material refinement reduces the resting hull from 704 to 640 panels and passes
the existing 3.4% diagonal screen there. It fails a held-out pose at 3.53% surge
error. The reference passes all three poses. The report retains all six modes,
the failures, preparation cost and complete case timings; these first timings
are not warm benchmarks. Neither smaller layout replaces the reference.

The repeated six-RHS timings reuse unchanged operators and factors. They are
not moving-stage timings. Full candidate time includes geometry, construction,
factorization and the repeated solves; independent Capytaine work is outside
that timer. Memory includes sparse storage/fill and lifetime process peak RSS,
which can include earlier cases and the reference solver.

Both static reports remain `qualified=false`; they do not establish moving-canoe
time evolution, outgoing-wave behavior or real-time performance. The intended
deliverable remains full nonlinear six-DOF evolution and Capytaine validation.
The unusable measured runtime makes optimization the immediate priority;
adaptive and partial updates must preserve accuracy and count estimator,
refresh and rejection costs.

## Nonlinear moving canoe

```text
uv run --locked --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.moving_canoe --duration 1
uv run --locked --group correlation --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.validate_moving_canoe --capytaine
```

The coupled experiment advances all six body velocities, world position and
quaternion, free-surface heights and potentials. Hull/free-surface contact nodes
slide on original physical edges and facets with exact shared coordinates.
Triangles are subdivided after clipping; interior wet nodes redistribute within
their physical facets as the waterline moves. Fine waterline heights remain
independent. Sparse geometry Jacobians carry this coupled motion through the
surface and pressure calculations. ALE graph kinematics
and Bernoulli transport use the current surface geometry. Every stage rebuilds
the current BEM/FEM operator, integrates total nonlinear pressure, and solves
body acceleration with the current six acceleration-potential columns. Two
offset field solves differentiate the known moving-node potential; no frozen
added-mass or hydrostatic force is added to the pressure calculation.

On the earlier material-edge chart (`7be7048`), a one-second GPU trajectory with nonflat initial water, nonzero potential and
all six initial body velocities completes on the n=12 mesh. Independent tests
check moving-edge contact, exact uniform flow on tilted water, clipped
hydrostatics, pressure transport and world/body mechanics. The complete stage
retains hydrostatic equilibrium to about 1e-11 in its rates; halving its
directional derivative step changes acceleration by less than 2e-11 relative.
Accepted states, water volume and mechanical energy are saved by the runner.

This remains an unfinished validation baseline. Contact topology is fixed;
reaching an original physical vertex or folding a mesh triangle raises rather than silently changing or clamping the
waterline. The tank is finite and reflecting. Remeshing, physical failure-event
coverage, outgoing waves, space/time/energy convergence and finite-frequency
Capytaine comparisons remain required. The same-hull instantaneous Capytaine
comparison has up to 5.94% diagonal error on n=12; finite-domain and spatial
effects are not separated. Component and short-run timings do not qualify
complete nonlinear performance. Optimization is deferred until validation.

On that earlier chart, at n=24 with material subdivision 3, the same-hull instantaneous comparison
reaches a maximum diagonal difference of 1.35%. Both numerical representations
are refined; this does not isolate either solver's error. The refined nonlinear
moving case also completes. Those numerical results have not been transferred
to the new chart. Current commands use `--wet-subdivisions` for refinement after
clipping, and output names use `wet` to distinguish the layouts. Independent
geometry tests at subdivisions 2 and 3 cover ±1 cm heave and ±1° roll/pitch,
matching clipped area and volume while retaining positive triangle areas.
This establishes those finite poses, not arbitrary motion or remeshing.

`moving_canoe --initial-scale 10` increases every initial pose, velocity,
surface-height and potential perturbation tenfold. Use the same initial scale
when comparing timestep or spatial refinements. These are unforced release
cases; they do not replace wave-excitation validation.

The physical-facet n=12/wet2 release at initial scale 10 completes one second,
with maximum volume change `1.04e-8 m³` and energy change `0.00397 J`, or
`1.16%` of its `0.341 J` initial disturbance energy. The separate
`audit_canoe_energy` command differentiates total energy in the actual nonlinear
RHS direction, independently of time integration. At that initial state it finds
approximately `+0.00993 W` over three derivative steps: a resolved spatial/work
consistency error requiring refinement. A successful trajectory does not yet
qualify nonlinear conservation. The large negative tank gravity baseline is
excluded from the drift normalization.
At n=24/wet3 the same audit increases to approximately `+0.01248 W`, or
`3.60%` per second. This refinement fails to establish energy/work convergence;
the nonlinear baseline remains unqualified while pressure and surface work are
investigated.
Tightening the one-second release from `rtol=1e-4, max_step=0.01` to
`rtol=1e-5, max_step=0.005` leaves energy drift at `1.19%`, while reducing
volume drift to `1.34e-9 m³`. Final position differences are below `3e-8 m`;
this temporal agreement does not resolve the spatial/work failure.
The energy audit now separates fluid/body kinetic and gravity work, subtracting
gravity references per element before summation. Body power agrees with
`velocity · hydrodynamic_load`; the fluid side accounts for the discrepancy.
`audit_canoe_shape_work` independently reproduces the entire coarse discrepancy
as operator shape work minus reconstructed physical-boundary work. Known affine
flow controls improve under refinement, so this does not establish a universal
sign error. Non-affine trace accuracy, velocity reconstruction and interface
work remain under investigation; no conservation correction is applied.
Raising outer quadrature from order 3 to 6 changes the coarse shape defect by
less than `0.08%`. Recovering the auxiliary P0 flux also fails: its hull-normal
RMS difference from the prescribed rigid velocity is `85.4%` on n12/wet2 and
`89.0%` on n24/wet3, while affine controls confirm its sign convention.
Using it in boundary work worsens agreement and would introduce hull slip.
The symmetric operator's imposed weak flux and its auxiliary P0 density are
distinct discrete quantities; substituting one for the other is unsupported.

The separate `canonical_canoe` experiment derives a reference from the same
discrete fluid action, retaining geometry dependence of the body/surface
coupling and momentum maps. It differentiates the complete Hamiltonian at fixed
canonical momenta by finite differences, so it is deliberately a slow reference,
not an adopted timestepper. At n12/wet2 and initial scale 10, two derivative
steps change decoded acceleration by `1.08e-8` and potential rate by `3.36e-8`.
Its height rates match the current kinematics, while heave acceleration changes
from `-0.05327` to `-0.05448 m/s²` and pitch acceleration from `0.4551` to
`0.4569 rad/s²`. Tests independently check physical energy, the momentum
map, waterline kinematics and the canoe's free rotational dynamics. Numerical
Hamiltonian closure alone does not qualify pressure, waves or body response.
An independent physical-energy derivative along the decoded reference motion
is `-1.96e-8 W` at the finer directional step, compared with `+0.00993 W`
for the current pressure-based stage. A stationary-action contraction matches
the reference force gradient in all six body directions and two nonconstant
surface-height directions. This permits direct geometry derivatives without
differentiating the implicit field solution. Neither calculation is a validated
time trajectory or an adopted simulation baseline.

```text
uv run --locked --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.canonical_canoe
```

After producing six fitted radiation reports per frequency,
`validate_radiation <reports...> --output artifacts/radiation-family.json`
compares all 36 added-mass and damping entries using the formal nondimensional
family tolerances. `--reference <report>` can reuse a matching Capytaine result
from a short radiation run. Incomplete mode sets or inconsistent mesh/domain
settings are rejected. Passing coefficient gates alone does not establish
amplitude, timestep, spatial, domain, transient or reference convergence.

The radiation runner uses this same nonlinear geometry and pressure path while
prescribing one body motion for validation:

```text
uv run --locked --group correlation --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.canoe_radiation --mode heave --periods 3 --ramp-periods 1 --fit-periods 2 --capytaine
```

All six prescribed modes are supported. The runner fits all six force/moment
components, including signed couplings, after startup; it removes independent
hydrostatic restoring when extracting added mass and damping. A smooth exterior
sponge damps both canonical surface variables and is zero around the pocket.
This is an artificial absorbing tank, not an exact outgoing boundary. The
Capytaine reference uses the same immersed hull and 4 m water depth and saves
the complete reference matrices. The small default run is only a smoke check.

Finite-frequency qualification requires transient decay, amplitude halving,
space/time refinement, and stable results when the sponge and outer domain
change. `--exterior-scale` adds outer cells while preserving pocket coordinates
and water-cell spacing, separating outer-boundary sensitivity from resolution.
No finite-frequency qualification is claimed yet.

## Nonlinear fluid prerequisite

```text
uv run --locked python -m tests.experiments.free_surface.shape_check --pocket
uv run --locked --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.fluid_baseline --n 24 --periods 1 --rtol 1e-5
```

`pocket_evolution.py` assembles the current BEM/FEM operator at every nonlinear
fluid stage. `fluid_baseline.py` advances it with adaptive RK23, saves accepted
states and energy, counts rejected steps, and checks a deformed GPU stage
against CPU plus gauge and volume invariance. `--cpu` needs no CuPy. The CUDA
kernel uses float64 throughout and the same analytic source integrals and
quadrature as the CPU reference. Cholesky failures and nonfinite rates raise.

Independent quadratic wave terms and directional energy derivatives converge
under spatial refinement. For the baseline wave, directional surface-work
discrepancy decreases from 0.396% to 0.102% at n=12/24; this is neither exact
discrete energy conservation nor a per-step energy loss. The n=24 GPU/CPU rate
comparison differs by less than 2e-15 in the projected mass norm. At time
tolerance 1e-5, a one-period run changes energy by -0.00753% and volume by less
than 2e-15 cubic metres. These are short closed-tank checks, not long-duration
qualification.

Reports distinguish integration time from process time including imports,
preparation, CPU/GPU verification, and state output. The verification warms
the GPU. **These are body-free component timings, not full nonlinear canoe
performance.** Moving wetted geometry, six-DOF pressure coupling and outgoing
waves remain required. No operator-reuse or mesh optimization is adopted here.

The production architecture is under review against the real-time requirement in
`docs/numerics.md`. These dense experiments remain numerical references.
The separate mesh-cost study isolates coplanar solver-panel subdivision from
changes to the canonical curved hull geometry:

```text
uv run --locked --group correlation python -m tests.experiments.capytaine_mesh_cost --case 2 0 --case 2 1
```

Omitting `--case` runs three successive panel resolutions of the fixed level-2
canoe; larger predeclared cases are selected explicitly after inspecting cost and
error. Each case uses a fresh worker and Numba cache, records first and five warm
Topwater instantaneous-inertia solves, and retains complete Capytaine finite-
frequency added-mass/damping matrices and total complex excitation at three
frequencies and headings. Excitation includes Froude-Krylov pressure. See
`artifacts/capytaine-mesh-cost.json` and individual case files.

Inspect roll and other individual modes as well as the family-wide tolerance:
its absolute floor can conceal large relative errors in weak modes. Differences
between physical hull refinements also contain remaining solver-panel error.
These are static resolution/cost results, not moving-step benchmarks or finite-
frequency Topwater correlation. They do not qualify near/far coupling or real-time
performance. Partly submerged physical faces are clipped, not marked wholly wet
or dry; flat-water clipping already implements that operation.

On the development CPU (Python 3.12.8, one BLAS thread), the first fixed-geometry
study produced the following results. All rows use the same 44 physical wetted
triangles; only coplanar numerical subdivisions change.

| Solver panels | Median of five warm Topwater instantaneous-inertia solves | Capytaine roll added-mass change to next row | Roll damping change to next row | Roll excitation change to next row |
| --- | --- | --- | --- | --- |
| 44 | 1.27 ms | 27.3–28.2% | 27.7–29.5% | 15.8% |
| 176 | 7.72 ms | 3.9–4.6% | 6.1–6.8% | 3.8% |
| 704 | 126.6 ms | — | — | — |

Added-mass/damping changes are relative diagonal differences across the three
frequencies; excitation is maximum complex difference divided by the finer
roll family's peak magnitude. The finest row is a comparison, not established
truth. Aggregate family tolerances pass even for 44 versus 176 panels, concealing
the roll errors above. Same-mesh Topwater/Capytaine instantaneous roll added-mass
differences decrease from 2.38% to 0.34% to 0.092%; that agreement does not establish
mesh convergence. These results justify mode-specific resolution studies, not a
blanket 704-panel requirement. Runtime includes assembly/factorization and six
auxiliary right-hand sides on cached clipped panels, and excludes wave evolution,
changing geometry and dynamic pressure. No timestep rate follows from this table.

The additional level-3 geometry case has 140 wetted panels and a 4.66 ms warm
median. Comparing it with the subdivided level-2 hull changes both physical
geometry and remaining panel error; the added-mass and damping family criteria
fail. It therefore cannot establish geometry convergence or supply a faster
equivalent mesh. CI retains all four selected cases; reproduce them with
`--case 2 0 --case 2 1 --case 2 2 --case 3 0`.

This checkout-only experiment tests a concrete path from nodal free-surface
state through moving geometry, a BEM solve and nonlinear time derivatives.
It is not an installed Topwater solver, a completed Stage 1/2, or an open-water
wave-body simulation. The rectangular tank is a fluid verification domain;
`hull_check.py` uses the existing canonical canoe, not another test hull.
`canoe.py` joins its actual wetted shell to a surrounding flat free surface and
tests static acceleration coupling on the complete finite fluid domain.

The earlier centroid-constant prototype cannot simply be wrapped in nodal
evolution: averaging vertex potential onto triangles exactly erases nonconstant
three-color patterns. Consistent flux projection improved smooth modes but could
not recover information lost upstream. A direct dual-patch collocation probe
removed those nulls but produced complex eigenvalues. Neither is a runtime option
here. `verify.py` retains the inexpensive transfer counterexample.

## Data and solve ownership

`Tank` owns the reference closed fluid mesh, node/face index maps and fixed
horizontal P1 mass/gradient/quadrature data, including one reusable mass
factorization. State is two vectors, `eta` and `psi`,
at the same surface vertices. No surface-potential averaging precedes the solve.

Each call to `rhs` moves the top vertices to `eta` and stretches side columns
vertically from the fixed bottom. Connectivity is validated once; stage calls
reject nonfinite state or a graph touching the bottom and build only the actual
stage panels and closed-boundary Laplace energy matrix, eliminates zero-flux
wall/bottom nodes,
and returns both state rates, integrated elevation rate and fluid energy.
No operator, factorization or potential mean is cached across moving stages.
The flat linear verification can reuse a frozen operator because its geometry
does not change. Assembly consumes the owned single-layer buffer after forming
the curl contribution and accumulates the remaining contribution into that same
energy buffer. With `V = L @ L.T`, the remaining contribution is formed as
`Z.T @ Z`, where `Z = solve_triangular(L, H)`. This is the same Schur complement
with one triangular solve instead of two. Dense assembly and factorization still
dominate moving stages.

The potential space is globally continuous P1; auxiliary flux uses P0 panels.
Use the existing positive Rankine kernel and outward **fluid** normals. The
symmetric Steklov operator is

```text
H = K + 0.5*mixed_mass
S = W + H.T @ solve(V, H)
W = sum_axis(surface_curl_axis.T @ V @ surface_curl_axis)
```

`V` is the Galerkin single layer, `K` the mixed P0-test/P1-trial double layer,
and `W` the weak hypersingular operator. See the
[boundary-operator definitions](https://bempp.com/handbook/api/boundary_operators.html)
and [symmetric Steklov identity, equation 2.6](https://link.springer.com/article/10.1186/s13661-016-0570-2).
The surface-curl formula requires the complete closed boundary and global
continuity. It cannot be copied onto an isolated open hull/free-surface patch.
Facet-interior Galerkin integration does not need vertex collocation jumps.

Source integrals and their first moments use a fused scalar-scratch kernel;
outer target integration uses Gauss-Duffy quadrature. `V` averages its two
source/target quadratures of the same symmetric integral, reporting their
discrepancy before averaging. This is not a certified quadrature error estimate.
Uniform target quadrature retains the closed-boundary constant identity;
pair-dependent target orders did not. No row-sum or flux correction is applied.

## Nonlinear graph equations

Let `M` be the fixed projected-area P1 mass matrix. After eliminating wall
potentials, `S_FS @ psi` is the weak physical flux load `q_normal*dS`.

```text
eta_t = Gpsi = solve(M, S_FS @ psi)
psi_t = -g*eta - |grad_xy(psi)|^2/2
        + (Gpsi + grad_xy(eta).grad_xy(psi))^2
          / (2*(1 + |grad_xy(eta)|^2))
```

The dynamic equation is integrated with P1 test functions before solving `M`.
Horizontal gradients are facet gradients, not averaged nodal gradients. Graph
flux is interpolated at quadrature points. Projected area matters: physical
normal velocity alone is not the elevation rate on a sloped graph. The selected
atmospheric gauge is retained, including constant potential-rate terms.

This weak dynamic equation is a consistent continuous formulation, but it is
not exactly the shape gradient of the discrete BEM energy. `shape_check.py`
measures that difference independently. Spectral positivity of the symmetric
construction does not by itself prove physical accuracy.

## Body-coupling constraint

The current evolution has no moving body. A coupled extension must use the same
body flux map for the boundary condition and pressure loads. With six body
motions `v`, define `B[i,j] = integral_body(N_i * Q_j * dS)`, where `Q` uses
fluid-outward normals and rotations about the declared load reference. This is
a nodal weak load, not a separate panel-centroid pressure rule.

Partition nodes into prescribed free-surface potential `F` and unknown potential
`U`. Natural wall/bottom flux is zero. The mixed solve and remaining flux are

```text
S_UU @ phi_U = B_U @ v - S_UF @ phi_F
free_surface_flux = S_FU @ phi_U + S_FF @ phi_F - B_F @ v
generalized_potential_load = B_U.T @ phi_U + B_F.T @ phi_F
```

Retain the full `B`, including waterline nodes in `F`. Dropping its fixed rows
when computing the instantaneous inertia matrix alone is valid for zero
free-surface acceleration potential. Those rows are still needed to recover
free-surface acceleration flux, and for a general moving body. For that auxiliary
solve, `A = density * B_U.T @ solve(S_UU, B_U)` and the acceleration reaction is
`-A @ acceleration`. One factorization handles all six right-hand sides.

This algebra constrains shared ownership; it does not qualify the quadrature,
waterline state transfer, nonlinear pressure or open-water domain treatment.
Finite rigid-wall tanks and Capytaine's infinite-depth open-water problem have
different boundary conditions. Their agreement cannot replace separate outer
mesh, domain-size and quadrature studies on identical wetted hull geometry.

`canoe.py` exercises this auxiliary solve. Its caller owns heave equilibrium and
the physical wet geometry; domain construction owns only the joined mesh and
body/free-surface face indices. One exact P1 load map drives both forcing and
pressure integration. Four additional right-hand sides, `phi=x,y,z,1`, exercise
nonzero free-surface potential, prescribed body translation and matching outer
wall flux. All ten columns share the same factorization. The `phi=z` column has
independently known potential, Dirichlet energy and hydrostatic load. The Dirichlet energy is
twice unit-density kinetic energy and equals fluid volume for this field.
The other three columns add no assembly or factorization.

`coupling.py` checks their surface elevation rates against the independently
known values `0,0,1,0`. Removing body and wall load rows at shared surface nodes
would corrupt those rates. The check also compares `B.T @ phi` against independent
triangle integration, both for the exact trace and for the solved potential.
For prescribed pressure `p=-phi`, the force/moment is the negative of that
potential-load vector. This is not Bernoulli pressure of uniform flow.

The mesher supports the canonical single convex horizontal waterline at `z=0`,
strictly inside the box. It uses unconstrained Delaunay triangulation, checks
convexity and rejects lost shared edges through full closed-mesh validation.
It neither welds coordinates nor supplies general moving-waterline remeshing.
Two hull refinements, off-grid translation, yaw and heel verify exact closure,
fluid volume, free-surface area and pressure/moment integration. The rectangular
`Tank.rhs` still cannot evolve this perforated surface: moving connectivity and
fluid-state transfer remain required.

## Reproduce

Run from a locked checkout; outputs go to ignored `artifacts/`:

```powershell
uv run --locked python -m tests.experiments.free_surface.verify
uv run --locked python -m tests.experiments.free_surface.linear
uv run --locked python -m tests.experiments.free_surface.evolution --fine
uv run --locked python -m tests.experiments.free_surface.integrate
uv run --locked python -m tests.experiments.free_surface.integrate --rtol 1e-5
uv run --locked python -m tests.experiments.free_surface.shape_check
uv run --locked python -m tests.experiments.free_surface.coupling
uv run --locked python -m tests.experiments.free_surface.hull_check
uv run --locked python -m tests.experiments.free_surface.benchmark
uv run --locked python -m tests.experiments.free_surface.benchmark --canoe
uv run --locked --group correlation python -m tests.experiments.free_surface.canoe --study
```

`verify` reports `qualified=false` and remaining work; command completion is not
a qualification gate. `integrate` rebuilds the BEM at every RK23 stage and takes
about a minute for the finer short case on the development host. The benchmark
uses fresh workers, empty compilation caches and five warm samples. Its scope
is geometry through the flat solve and independent diagnostics, **not** a moving
timestep, whose actual cost is recorded by `integrate`.

The canoe comparison generates pinned Capytaine references directly from the
supplied physical triangles and CG, outside candidate timing. It never reads
another test's report. CI retains its six-case outer-mesh/domain/quadrature study
beside the existing Capytaine reports; `qualified=false` remains explicit.
The canoe benchmark uses the same fresh-worker harness, including equilibrium,
clipping, joined-domain geometry, assembly, ten-RHS solve, pressure and affine
diagnostics. Validated input construction, imports and Capytaine are excluded.

On the Windows development host (Intel family 6 model 170, 22 logical CPUs,
Python 3.12.8/NumPy 2.5.4/SciPy 1.18.1/Numba 0.68.0, one BLAS thread), the
complete flat verification path measured:

| Tank panels | Empty-cache first call | Five-run warm median | Fresh-worker peak RSS |
| ---: | ---: | ---: | ---: |
| 320 | 1.57 s | 0.115 s | 155 MiB |
| 1280 | 3.30 s | 1.89 s | 213 MiB |

The report retains individual samples, hardware/thread information and explicit
matrix-buffer byte counts. These timings include independent diagnostic work;
they are not a simulation throughput or accuracy-qualified performance claim.

The subsequent ownership change averages `V` in place, consumes its buffer for
the factorization after forming `W`, and accumulates the final energy into `W`.
The mixed trace `H` survives for its left multiplication. In four fresh workers
at 1280 panels, old/new/new/old peak RSS was 212.6/194.7/193.9/212.6 MiB.
Corresponding five-run warm medians were 1.73/1.71/1.67/1.80 s. Alternating
same-process measurements overlap, so this establishes a memory reduction without
a claimed runtime speedup. The compared energy matrices are bit-for-bit identical.
The buffer byte report describes `V` after factorization and the final energy
buffer; hypersingular and energy storage are no longer distinct allocations.

## Evidence and limits

With 460 wetted canoe triangles, increasing outer horizontal resolution from
`(8,4,4)` to `(16,8,4)` changes diagonal inertia by several percentage points.
Doubling domain width, length and depth at the finer horizontal spacing changes
those diagonals by at most 0.15 percentage point. Doubling vertical resolution
alone changes them by less than 0.0002 percentage point. Third- to sixth-order
quadrature on the coarse outer mesh changes them by at most 0.086 percentage
point. Within these probes, outer horizontal discretization dominates the
remaining difference from the open-water reference; this is not qualification.

For that canoe, refining the outer mesh halves the independent affine potential
error, while relative solved Dirichlet-energy error decreases from 0.0198% to 0.00849%.
The solved affine pressure load also converges toward independent hydrostatics;
tests check the complete force/moment vector, not a symmetry-repaired matrix.

For the four nonzero-surface/body-motion fields, whole-surface continuous P1
elevation-rate RMS errors decrease from `[0.00318, 0.00267, 0.00426]` to
`[0.00154, 0.00181, 0.00156]` for the three coordinate fields. Each has unit
velocity scale. The constant field remains at roundoff, and integrated flux
matches the independently known free-surface area times vertical velocity.
These manufactured problems prescribe the required outer-wall flux: they do
not establish a wave solution with impermeable stationary walls.

Fresh-worker canoe-path warm medians at 864/1820 complete fluid faces are
0.772/3.53 s, with 220/287 MiB peak RSS on the host described above. Empty-cache
first calls are 7.38/10.30 s. These include all four coupling-field diagnostics;
they are static auxiliary solves, not coupled timestep throughput. The dense
method has not demonstrated the cost or accuracy required for production.

The source first-moment vector is computed once per target point/source face
and reused for all three P1 corners. Two warmed, interleaved old/new comparisons
on those canoe domains gave complete-path medians of `0.779 -> 0.688 s` and
`3.140 -> 2.898 s`, then `0.931 -> 0.807 s` and `3.813 -> 3.526 s`. Raw single-
and double-layer matrices, inertia and affine outputs were bit-identical.
These paired warm results support the arithmetic change; separate fresh-run
times vary, and no cold-start or memory improvement is established.

On the 320/1280-panel tank, third-order outer quadrature gives interior **nodal** velocity
errors of 6.16%/0.518% for the axial mode and 14.0%/1.32% for the oblique mode.
The comparison region is fixed physically. Sixth-order quadrature at 1280 panels
changes these to 0.505%/1.30%; spatial error dominates that difference. Independent
whole-surface quadrature gives substantially larger continuous-field errors:
21.2%/6.07% for the axial mode and 25.8%/7.95% for the oblique mode. Nodal accuracy
does not remove the P1 interpolation error between nodes.
Constant potential, net flux and the integrated graph volume rate are preserved
to roundoff without projection. Independent linear-field energy errors decrease
with quadrature and mesh refinement. The generalized spectrum has one constant
zero mode and positive remaining modes.

Actual frozen-flat BEM evolution over five periods reduces continuous elevation
space/time RMS error from 31.5% to 6.31% under that mesh refinement. The separate
comparison against a P1-interpolated nodal reference gives 26.2%/3.17% and is
explicitly labelled as such in the report. Tightening RK tolerance
from `1e-7` to `1e-9` leaves the finer spatial error unchanged while reducing
relative energy drift from about `1.7e-5` to `1.7e-7`.

Moving-graph nonlinear terms are checked independently by expanding
`G1(eta)psi = -G0(eta*G0psi) - div(eta*grad(psi))` and the quadratic dynamic
terms. For the non-cancelling `kd=1` standing mode, spatial refinement reduces
quadratic kinematic error from 6.65% to 0.667% and dynamic error from 1.73% to
0.353% **at nodes in the fixed interior comparison region**. Reports also include
independent continuous whole-surface errors. The `kd=4` case is retained; its quadratic kinematic term nearly
cancels and is reported both relatively and against a natural dimensional scale.
For `kd=1`, continuous whole-surface quadratic errors over three refinements are:

| Tank faces | Surface nodes | Kinematic error | Dynamic error |
| ---: | ---: | ---: | ---: |
| 320 | 45 | 6.68% | 5.05% |
| 1280 | 153 | 2.41% | 1.57% |
| 5120 | 561 | 0.868% | 0.513% |

The finest case reaches the one-percent whole-field target for these nonzero
quadratic coefficients. Its integrated elevation rates are below `1e-15 m^3/s`.
It does not qualify the general field, long nonlinear evolution or open water.
At `a=0.02 m`, `k=0.25 m^-1`, the complete positive/negative-state evaluation
measured 70.4 s and 793 MiB peak RSS on the development host above. That single
invocation includes fixed data construction, both moving geometries, both BEM
assemblies/solves, nonlinear rates and independent continuous diagnostics. It is
not a five-sample throughput benchmark. Wall time remains a major limitation.
`evolution --fine` retains this third refinement with the other cases in CI;
the existing physical test checks whole-field convergence and the one-percent
target, separately from interior-node accuracy. Reports remain `qualified=false`.
Halving amplitude separates perturbation truncation from spatial error.

The independent kinetic-energy shape derivative differs from the dynamic RHS
by 36.9%/3.69% on those two meshes for a smooth second-harmonic perturbation;
halving the finite-difference displacement does not change the result. In the
short quarter-period nonlinear evolution, tightening RK tolerance from `1e-3`
to `1e-5` reduces the fine energy change from -0.0951% to -0.00820%; the coarse
change at `1e-5` is -0.0499%. Volume changes remain below `1e-15 m^3`. Nonlinear
energy consistency therefore still needs spatial refinement and long runs.

On the canonical canoe, bounded Neumann nodal potential RMS errors at
128/512/2048 panels decrease for linear, quadratic and external-source fields.
At 2048 they are approximately 0.0026%, 0.170% and 0.0628%. The close exterior
source remains harder: fine potential error is 1.95% and energy error 3.13%,
with a nonmonotone coarse energy sequence. This is not a general field gate pass.
The canoe closure is only a mathematical domain for these checks.

Remaining requirements include tolerance-controlled self/near/far Galerkin
integration with cancellation handling, allocation budgets and conditioning
checks; longer nonlinear space/time/energy studies; open deep-water generation,
absorption and domain convergence; instantaneous wetted-body coupling and
finite-frequency Capytaine validation. There is no breaking/flooding solver,
checkpoint format or public simulation API in this experiment.

## Canonical SoCal incident-field cost

`incident_baseline_benchmark.py` evaluates the selected observed SoCal spectrum
with seed 42 and the 24 ft by 9 ft boat's scale. It reports measured-band
frequency/derivative quantiles, reconstructed directional spread, and nine
analytic linear incident outputs at 256/1,024/4,096 query points. It compares
48/96/192/384-component realizations with per-band variance retained. Different
component counts have different phases/directional samples; they are not
pointwise approximations of one fixed history.

```powershell
$env:OPENBLAS_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'
uv run --locked --python 3.12 --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.incident_baseline_benchmark
```

The bounded driver refuses to overwrite
`artifacts/canonical-sea-incident-benchmark.json`. GPU timing is optional if its
runtime is unavailable. Warmed device-resident and transfer-inclusive times are
reported separately; compilation/import time is separate. On this machine,
384 components at 4,096 points took about 2.6 ms including transfers, so reducing
the incident spectrum is unnecessary for this component's cost. The independent
CPU control checks potential gradients and linear free-surface identities.

Query clouds are workload controls, not fluid/hull meshes. These timings exclude
the fluid operator, total pressure, body motion and initialization. Neither the
energy bandwidth nor this kernel benchmark establishes spatial resolution,
force/motion accuracy, safe conditions or end-to-end simulation throughput.

## Current boat BEM factor sharing

`boat_current_factor.py` uses the actual 2500 kg boat's flat-water boundary,
prescribed nonzero six-axis velocity and the observed sea's potential trace.
This finite-box snapshot is an operator/shape-kernel check, not a nonlinear
sea initialization or deep-water simulation. The retained `prepare_energy`
result owns the current energy, single-layer factor and mixed trace; matching
force evaluation reuses that setup. Exact coordinate, face and quadrature checks
reject stale geometry. The retained host energy is read-only. Legacy
`assemble_energy` and default density 1000 retain their existing behavior;
the boat check explicitly uses density 1025.

```powershell
$env:OPENBLAS_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'
uv run --locked --python 3.12 --with meshpy==2026.1.1 --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.boat_current_factor --output artifacts/boat-current-factor-reuse.json
```

At 492 panels/248 vertices, fresh and shared-setup gradients agree within
1.3e-13 relative norm. Two directional energy differences give relative errors
1.3e-8 and 5.3e-9, with fixed outer geometry and nonzero near-hull deformation.
Density scaling and stale-setup rejection pass. The saved mixed field has
relative residual 2.3e-14. One warmed force pair costs 0.500 s fresh versus
0.445 s retained; these two observations are not a throughput distribution.
Warm energy assembly costs approximately 0.016 s and the mixed solve 0.004 s.
Retained GPU factor/trace storage is 2.91 MB, excluding runtime and other arrays.

A subsequent CUDA-event profile places 0.402 s of a 0.449 s force call in the
18-lane pair derivative kernel; host geometry Jacobians take about 0.044 s.
This motivates testing the existing acceleration-potential/Bernoulli pressure
route, which needs directional field updates instead of every vertex derivative.
The old hybrid canoe pressure route failed spatial work convergence, so a pure
BEM boat pressure path must earn force/work convergence. Neither route is
qualified by this current kinetic-kernel check. Saved evidence:
`artifacts/boat-current-factor-reuse.{json,npz}` and
`artifacts/boat-current-force-profile.json`.

### Full boat in the observed nonflat sea: geometry reference

```powershell
python -m tests.experiments.free_surface.boat_world_boundary --output artifacts/boat-world-boundary-seed42.json
```

This bounded exact-arithmetic reference constructs the complete closed fluid
boundary around the 2500 kg boat at its mean-water equilibrium pose, using the
unscaled seed-42 SoCal elevations on fixed world P1 parents. Child free faces
inherit parent height and potential coefficients exactly; the saved potential
trace samples the incident Airy potential at the actual graph height. It is a
prescribed snapshot, not an initialized nonlinear wave/body field.

The checked run (`artifacts/boat-world-boundary-seed42-checked.{json,npz}`)
completed in 1.89 s with 185 vertices and 366 triangles. Rim clearance was
0.272 m; independent volume, first-moment and affine rigid-work errors were
below 2e-13 in their respective units. The driver checks numerical acceptance,
rim immersion and loss of contact, and preserves partial failures. Exact
transverse symmetry in the authored mass moments prevents roundoff from creating
sub-ulp sliver intersections at the centerline. This Fraction-based mesher is a
correctness reference; its time is not a production solver benchmark. The 2 m
bed and finite walls remain diagnostic boundaries, not open-ocean physics.

### Current nonflat fluid field and fixed-world surface rates

```powershell
uv run --locked --python 3.12 --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.boat_world_field --snapshot artifacts/boat-world-boundary-seed42-checked.json --output artifacts/boat-world-field-seed42.json
```

One current mixed factor supplies the prescribed field and six acceleration
responses. Offset field evaluations can omit those six columns. Free-surface
reaction subtracts the prescribed body/wall/bed loads even at shared contact
nodes. Recover fine cut-mesh flux before the nonlinear Bernoulli calculation;
only then project the rates onto exposed fixed-world parent coefficients.
Both mass matrices must pass the declared algebraic conditioning budget.

The checked run (`artifacts/boat-world-field-seed42-checked.{json,npz}`) has
relative mixed residual 1.7e-14 and added-mass reciprocity error 3.8e-17. Warm
energy assembly takes 12.7 ms, field plus six responses 21.0 ms, single field
15.9 ms and surface rates 1.2 ms; cold setup is saved separately. These are
individual component observations, not a throughput distribution or full solve.

The coarse 15-parent-coefficient state omits 33% of the fine kinematic rate in
projected-area mass norm on this snapshot. This is an explicit approximation
warning requiring enrichment and force/motion convergence, not an accepted
production resolution. Independent affine harmonic controls exercise sloping
and partially exposed parents; an unresolved local-flux control ensures that
projected-out motion still contributes to Bernoulli kinetic head.

### Material-point pressure derivative audit

```powershell
uv run --locked --python 3.12 --with "cupy-cuda12x[ctk]==14.2.0" python -m tests.experiments.free_surface.boat_pressure_audit --output artifacts/boat-pressure-seed42.json
```

The audit reconstructs four independent offset boundaries at `+/-h` and
`+/-h/2`, evolves the same fixed-world parent coefficients and incident boundary
time, and samples identical original-hull interior quadrature points. Missing
wet coverage fails explicitly. Only the base field computes six acceleration
responses. Each field must pass its mixed residual check; all material traces
and partial geometry checkpoints are preserved.

The first run (`artifacts/boat-pressure-seed42-h1e4.*`) completed in 11.2 s but
FAILED derivative stability. Halving h from 1e-4 to 5e-5 s changed the surge
load by 6970 N, despite field residuals below 3e-14. The norm of the odd
material trace difference stayed 0.0424 versus 0.0413 instead of halving.
The one-sided meshes had different persistent triangulations (400 versus 396
faces, with 366 at the base). This is evidence of a discrete mesh-limit jump;
it does not qualify the pressure path. The initial driver source is frozen
beside the report. Subsequent audits explicitly report a directional-load
budget of 0.1% of weight and weight times length. Passing that numerical budget
would still leave spatial, quadrature, surface-basis and ocean checks open.

The centroid regression then identified finite fan-center jumps of up to
0.215 m caused by counting a nearly coincident cutpoint as a full extra vertex.
Exact polygon area centroids remove that error. The next audit
(`artifacts/boat-pressure-area-centroid-h1e4.*`) reduced the surge-load step
difference to 908 N but still failed (3.70% of weight). Offline extrapolation
shows the one-sided finite-area triangles now share the same geometric limit;
independent cut-node unknowns still coalesce while retaining finite support.

A controlled probe on the four saved offset meshes changed only BEM target
quadrature order 3/6/12. The extrapolated trace-jump norm stayed approximately
0.00684/0.00686/0.00687; ordinary target quadrature refinement is not the cure.
Reports and probe source are saved under
`artifacts/boat-pressure-area-centroid-quadrature-probe.*`. The next candidate
removes submerged background-parent cuts from hull facets while retaining
the actual waterline and original physical facets. No pressure baseline has
been accepted from these component timings.

Exact union of submerged pieces by original physical facet then removes
24 artificial collinear hull points: 29 body polygons become nine, and body
triangles fall from 110 to 46. The full boundary has 153 vertices/302 faces;
each facet's signed area and first moment are preserved exactly. No coalescing
unknown body vertices remain in the saved one-sided limit; remaining coalescing
vertices are prescribed-potential free-surface nodes at the bow.

The union audit (`artifacts/boat-pressure-body-union-h1e4.*`) reduces the
surge-load step change to 13.5 N, but sway changes 50.5 N and yaw 246 Nm.
The h=1e-4 default was not a physical requirement. A separate numerical-noise
control at h=1e-3/5e-4 (`artifacts/boat-pressure-body-union-h1e3.*`) uses
sub-millimetre surface offsets and passes the global directional budget
(0.0246% of weight; force step differences 1.65/6.04/-0.011 N).
This is a finite-difference window, not a smooth-operator proof: yaw torque
still changes from 507 to 536 Nm, approximately 5.8% of that response.
All four step sizes remain part of the evidence; neither the global numerical
budget nor the small algebraic residual qualifies six-mode physical accuracy.

For the latter audit, warmed offset energy assembly costs 10-13 ms, single
field solves 13-17 ms, material trace sampling 2-3 ms and pressure integration
0.6-0.8 ms. The exact-arithmetic reference geometry still costs 1.7-1.8 s per
offset and dominates the complete 12.6 s audit; there is no end-to-end speed
claim. Peak process RSS is 552 MB and peak reserved CuPy pool storage 5.76 MB
(CUDA context/library memory is outside that pool metric). Coarse parent-state
unresolved flux is now 38%, so physical spatial/state convergence and fast
production geometry remain open work.

A frozen union-mesh gauge control (`artifacts/boat-pressure-body-union-gauge-probe.*`)
adds ten to every parent potential without changing physical flux. Body trace
changes are ten to within 1.63e-11; the resulting extrapolated gauge error is
2.48 million times smaller than the remaining derivative trace error. The
union mixed blocks have condition numbers about 42.8, compared with about 1e7
before eliminating artificial hull cuts. A constant-mode/row correction is
therefore unsupported by this evidence. Remaining bow-contact integration
and pressure/state refinement still require checks.

A second private frozen-mesh control enforces the exact zero principal-value
double-layer contribution between triangles with the same physical plane owner.
It changes traces by only 1.6e-11 and leaves the extrapolated remainder
unchanged (`artifacts/boat-coplanar-plane-probe*`). No shader change is adopted.
The remaining auxiliary P0 flux space still contains vanishing integration
triangles; it is distinct from the now-stable mixed potential unknowns.

### Matched-state boundary and state-resolution checks

`boat_pressure_audit --step 0.001 --refinement 1` (or `2`) uniformly subdivides
the SAME boundary and inherits the SAME fixed-world parent fields. Closed fluid
moments, physical facet moments and affine world/material traces are checked.
Each actual refined checkpoint is saved separately from its coarse construction
record and hashed in the report; pressure fields must use that refined mesh.

The three reports `artifacts/boat-pressure-body-union{,-refined1,-refined2}-h1e3.json`
are compared by `boat_pressure_convergence`, which verifies physical input
identity before tabulating all six loads/accelerations and derivative sensitivity.
Summary: `artifacts/boat-pressure-spatial-convergence.json`.

| Faces | Surge N | Sway N | Heave N | Roll Nm | Pitch Nm | Yaw Nm |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 302 | 607 | 323 | 30123 | -308 | 250 | 507 |
| 1208 | 686 | 365 | 29989 | -303 | 86 | 404 |
| 4832 | 736 | 363 | 29946 | -301 | -2 | 427 |

The finer differences contract, but pitch crosses zero and surge/yaw are not
settled. No six-mode physical pass is claimed. The 4832-face audit costs 25.5 s
including both derivative steps, peaks at 723 MB process RSS/1.43 GB reserved
CuPy pool, and takes about 1.9-2.7 s per offset operator assembly. Increasing
uniform resolution is a reference tool, not a practical runtime solution.
The 15-coefficient parent state still omits approximately 37% of fine flux.

Profiling that actual refined state exposed avoidable diagnostic cost: repeated
SVD condition numbers consumed most of the surface-rate time. LAPACK `dpocon`
now estimates the one-norm condition from the existing Cholesky factors, with
the same rank/finite and estimated-roundoff screening. Three warm observations
decrease from 0.281-0.313 s to 0.038-0.044 s, with bit-identical physical rates.
The condition measure is explicitly labeled; it is not the previous two-norm
SVD value or a certified error bound. Original/new source and actual field
replay are preserved in `artifacts/boat-surface-condition-benchmark.*`.

An independent 45-coefficient state preflight samples the SAME explicit sea at
a finer 8-by-4 background grid; this intentionally changes its P1 approximation
and is not part of the matched-state table. It passes geometry/current algebra
checks at 608 panels, with 0.261 m rim clearance and parent mass condition 266.
Omitted fine flux improves to 22.6%, still unqualified. Changing state capacity
must therefore accompany boundary-operator accuracy work. Evidence/source:
`artifacts/boat-world-state-8x4-preflight.*`.

Exact diagonal equilibration of the auxiliary single-layer factor was also
tested on the frozen union meshes: condition estimates become about 86, but
body traces change only 3.6e-15 and derivative sensitivity is unchanged.
That exact algebraic change is not adopted as a cure; probe source/results are
`artifacts/boat-pressure-body-union-equilibration-probe.*`.

# Fixed-world surface state for the canonical boat

Active investigation status and next experiments: [GitHub issue #30](https://github.com/amosborne/topwater/issues/30).
This document retains technical reasoning and historical evidence.

Scope: the 24 ft × 9 ft, 2,500 kg boat, water density 1,025 kg/m³,
gravity 9.81 m/s², and selected observed SoCal sea with seed 42. This is an
architecture/force-audit proposal; no new fluid solve or implementation has
qualified it. The failed canoe volume/chart is not the reference geometry.

The optimization target is this boat/sea problem, following the user's explicit
scope restriction. Substantially larger/smaller boats and substantially more
energetic, faster or slower seas are outside the target. Do not spend runtime
or add abstraction to cover them. Keep full nonlinear six-DOF motion, the
observed directional sea and physical failure checks in timing comparisons.
The target remains at most 20 wall seconds per simulated second on this machine,
with real time preferred; reduced work counts alone do not meet it.

The user has now supplied an output accuracy target: 0.1 degree numerical
convergence in EACH of roll/pitch/yaw, with a practical equivalent 1 cm target
for EACH translation (surge fore/aft, sway sideways, heave vertical), throughout
a full 60-second simulation. Shorter traces only screen candidate changes;
endpoint agreement alone does not satisfy the requirement. A 0.1 degree
rotation over a 4 m lever arm displaces a point about 7 mm. The user's final
scope clarification is ordinary motion tolerable for a person aboard; there
is NO fixed roll-angle cutoff. The earlier 20-30 degree discussion can guide
test coverage, but it is neither an operating limit nor a safety determination.
Motion severity depends on rates, accelerations and combined six-DOF motion,
not angle alone. This intended-use boundary does not define an automatic
human-safety classifier; do not invent numerical safety thresholds from it.
Preserve full six-DOF coupling while prioritizing
the user's roll assessment. Do not equate local quaternion/ODE
tolerances with this end-to-end target. Measure mesh, quadrature, time and
operator approximations against roll over representative wave/roll periods;
short weakly excited roll traces cannot establish it. `boat_motion_metrics.py`
derives roll/pitch/yaw diagnostics without changing quaternion evolution. Roll
uses projected world vertical in the body port/up plane, independent of heading;
periodic angle errors must be wrapped. Do not add robustness requirements for
dangerous/extreme motion, or invent an angle-only human-safety rule. Existing
flooding/contact and incident-wave physical/domain checks remain in force.
For comparisons, project horizontal position differences along/across the SAME
reference heading, retaining world-up heave; do not compare candidate/reference
translations in independently rotated frames. ZYX pitch is positive bow-down
in the existing bow/port/up convention; correlation must retain that sign.

**Recommendation:** investigate a surface-only variational BEM with an evolved
fixed-world surface state and freshly reconstructed wetted boundary. First use
fixed piecewise-affine world parents, retaining local surface response. Do not
start by imposing the incident wavelength cutoff on the scattered solution.
This removes the requirement that every fluid tetrahedron remain valid while
the real chine crosses the waterline; it does not prove BEM speed or accuracy.

## State and the restricted action

Let `Q=(a,ξ)` contain six body configuration coordinates and world surface-height
coefficients. Let `X(Q)` describe the reconstructed physical boundary. With
fluid-outward normals, define `D(Q) Qdot` as the integrated nodal boundary-normal
velocity load, and `E(X)` as the symmetric discrete Dirichlet-to-Neumann energy.
Write the full trace as `φ=Z(Q) z`, where `z=(u,χ)` contains algebraic hull/wall
trace unknowns and evolved free-surface potential coefficients. On the free
surface, `φ_F=S(Q)χ`; the values of fixed world functions at moving vertices are
not a constant nodal matrix.

The unreduced spatially discrete Luke action is

```text
L = T_rigid(a,adot) + ρ φᵀ D(Q) Qdot − (ρ/2) φᵀ E(X(Q)) φ − U(Q).
P = ∂L/∂Qdot = ∂T_rigid/∂Qdot + ρ Dᵀ Z z.
```

Eliminate only the stationary algebraic unknowns `u`; obtain the remaining
state/rates from this same action and momentum map. Body position and orientation
remain six physical freedoms. Converting body canonical momentum rates to world
acceleration still requires differentiating the current momentum map, including
the body rotation chart. A body `Pdot` is not by itself a physical force vector.

For `r=Eφ−D Qdot`, the restricted potential equation gives `Zᵀr=0`. It does not
give `r=0` on every original free-surface node. Holding `z` and `Qdot` fixed, the
fluid contribution to a configuration variation contains

```text
δL_fluid = ρ φᵀ(δD)Qdot − (ρ/2)φᵀ(δE)φ − ρ rᵀ(δZ)z.
```

Consequently the additional free-trace term is **`−ρ r_Fᵀ(δS)χ`**. The existing
fixed-nodal-`φ` shape gradient does not include it. The projected residual cannot
be used to discard it. Derivatives of the geometry velocity map (`Jdot`) and
rigid kinetic energy also remain; their force must come from this same action.

At a moving surface node, a fixed-world height/potential field satisfies

```text
δz = R(x,y)δξ + ∇η(x,y)·δXY,
δψ = S(x,y)δχ + ∇ψ(x,y)·δXY.
```

Reevaluating world basis values while using the old constant-`R` pullback omits
both these contributions. Explicitly time-dependent incident liftings introduce
additional time terms; adding a separate incident pressure to a total-potential
pressure calculation would double count excitation.

## When free-surface body transport cancels

For the exact world graph `X_F=(x,y,η(x,y))`,

```text
n·δX dS = (δz−∇η·δXY) dxdy = R δξ dxdy.
Θ_F = ρ χᵀ C(a,ξ) dξ,
C_ji = ∫_{D(a,ξ)} S_j(x,y) R_i(x,y) dxdy.
```

Thus mesh sliding at fixed `ξ` is tangential, and the free-surface `Cv` term is
zero. The projected water domain `D` excludes the current hull footprint. Its
moving-hole/Reynolds derivatives remain in `C`, gravity, body loads and forces.
When `C` is square and nonsingular, `p_ξ=ρ Cᵀχ`; differentiating it requires
`χdot=C⁻ᵀ(pdot_ξ/ρ−Cdotᵀχ)`. The body part of the one-form still contains the
hull's normal-motion term. It is not removed with `Cv`.

**The cancellation is not automatically exact for sampled cosines on flat BEM
triangles.** A planar triangle's interpolated gradient generally differs from
the analytic graph gradient at its vertices. Moving its vertices changes the
interpolated world graph even with fixed analytic coefficients. Using curved
graph normals in a planar-panel kernel would instead mix two discretizations.

The simplest consistent P1 construction is a fixed world XY parent mesh:

1. Define the height and potential fields as fixed-parent P1 functions. Cosines,
   if wanted, are interpolated once onto those parents.
2. Clip each physical hull facet against the piecewise-affine graph; preserve
   hull-parent and free-surface-parent identities.
3. Clip/retriangulate each exposed free-surface parent without crossing its
   fixed edges. Each child's graph and potential remain its parent's affine
   functions, so child-vertex sliding is exactly tangential.

This permits `Cv=0` in the discrete model. Otherwise retain the actual discrete
transport and basis derivatives rather than forcing cancellation. Very small
cut triangles still challenge BEM quadrature/conditioning; removing their
evolved independent states does not make their geometric integrals disappear.

The parent geometry mesh need not define one evolved state per parent node.
Functions supported entirely beneath the projected hull have zero exposed-water
mass; nearly covered supports create ill conditioning and cut-cell frequencies.
Choose a fixed trial space whose exposed restriction stays full rank and check
`C` conditioning. Global modes interpolated once onto the parents avoid strictly
hidden individual supports, but their local-response accuracy remains unproved.
Dropping/reinstating covered state variables dynamically would introduce a state
transfer problem; it is not an automatic consequence of clipped geometry.

Only world `η/ψ` coefficients and body canonical variables evolve. Hull/wall
potential coefficients are elliptic algebraic unknowns solved at the current
geometry. A changed wet triangulation therefore needs no invented transfer of
evolved nodal potential history. Nevertheless its discrete energy/force can
change: topology independence of the state is not exact discretization transfer.
Audit the limits on each side of chine/contact events and under refinement.

## Existing kernels and legal caching

`operator.py` constructs `E=W+HᵀV⁻¹H` on a complete closed boundary with P1 trace
and P0 auxiliary flux. Its uniform order-three quadrature and single-layer
symmetrization are experimental, not tolerance-controlled singular integration.
`cuda_shape_gradient.py` differentiates that same quadrature at fixed nodal
potential. Its kernel hardcodes density 1,000: the new case needs consistent
1,025 scaling in kinetic energy, one-form, gravity and pressure.

A small number of evolved surface coefficients does **not** reduce the number
of P0 auxiliary flux panels. Dense `V` factorization remains cubic in panel
count. The current force helper also rebuilds `V/H` and its factor independently;
sharing the already-current primal factor and stationary auxiliary field is a
valid simplification before approximating operators.

For scale only, 500 faces/250 vertices require about 3.5 MB for `V/H/E`; at
2,000/1,000 the same arrays require 56 MB, excluding factors, copies and shape
data. Dense Cholesky work scales from approximately 42 million to 2.67 billion
operations. Pair-gradient work also grows quadratically with 18 coordinate lanes.
These arithmetic envelopes are not measured timings on the new boat. Price
assembly, factor, solve and force together before selecting a panel resolution.

Caching rules:

- Fixed-geometry far/far **primitive Green integrals** can be reused exactly
  while that geometry, basis and quadrature are unchanged.
- An effective far block after eliminating `V` or other near unknowns generally
  changes with near/far coupling. It is not fixed merely because its own panels
  are fixed. Preserve the current coupling/elimination or measure its defect.
- A genuinely moving nonlinear far graph changes its Green integrals. A cached
  flat/linear far DtN is a physical approximation needing overlap/interface and
  force/motion error controls, not an exact geometry cache.
- Interface potential and weak normal flux must be shared consistently. A local
  closed-boundary energy cannot simply lose its artificial faces and still be
  the same operator. Reflection/open-domain errors require separate controls.

The incident field is already inexpensive: all 384 components and nine outputs
at 4,096 points cost about 2.6 ms including GPU transfers. Keep it intact. The
incident `λ≥2B` contract places no cutoff on body-generated radiation, contact
response or nonlinear harmonics; the earlier global-only canoe restriction
omitted 81–91% of the initial surface-rate mass norm. That warns against adopting
the same restriction here; it is not a measured error for the new boat.

## Smallest physical force falsifier before adoption

Use **one geometrically valid prescribed snapshot of this boat in the seed-42
SoCal field, with nonzero surface potential**, followed by a small pitch
displacement that exercises the problematic chine. A completed trajectory or
ramp is not a prerequisite for this one-stage comparison: the current elliptic
solve supplies the interior field from the prescribed surface traces and body
velocity. Retain the explicit pose, surface fields and velocity. A declared
wave-following initial pose is permissible; reducing the supplied sea to obtain
a passing geometry is not. Reject rim immersion/loss of contact. If the candidate
basis cannot represent the supplied traces, report the projection error, compare
both force paths at the SAME represented trace and use nested enrichment to
assess the approximation. Do not call unlike states a matched comparison.

If retaining the variational formulation, derive its full stage, including `δS`,
hole derivatives and the momentum-map derivative. Exact discrete Hamiltonian
structure is not a user requirement or a prerequisite for testing the cheaper
primitive-pressure formulation. Compare **six world accelerations and six
integrated hydrodynamic pressure loads** on the *same* current physical boundary;
physical convergence, rather than bitwise agreement between discretizations,
decides whether a formulation is usable:

1. Reuse one current factor for the six unit-world-velocity potential responses.
2. At two directional step sizes `h,h/2`, reconstruct and solve at the four
   offsets along the actual physical state rate. For the known-pressure part,
   hold the world twist fixed; evolve the world height and potential coefficients
   with their derived rates. Evaluate potentials at matching material hull
   samples, identified through original hull-facet coordinates. Do not difference
   arrays of newly numbered waterline nodes.
3. Recover `φ_t=φdot_ALE−w_mesh·∇φ`, add the six acceleration responses, and
   integrate `p=−ρ(φ_t+|∇φ|²/2+gz)` over the same wetted hull. Solve the joint
   six-body balance, including rigid angular gyroscopic terms and weight once.
   `body_pressure.py` supplies this balance/reconstruction contract, but its
   caller still owns the matched samples and correct time derivative.
4. Compare physical loads/accelerations, directional-step behavior and boundary
   refinement. Variational and direct pressure tractions need not agree to
   machine precision on a coarse mixed P1/P0 discretization; their discrepancy
   must decrease under meaningful quadrature/space refinement. Repeat the same
   snapshot on each side of the small pitch/chine event and verify force limits
   converge rather than jump with mesh numbering.

This costs a base geometry/factor and four offset configurations, not a long
trajectory. Budget those operations explicitly before running. A scalar action
finite difference is a useful implementation check, but cannot replace this
physical load comparison. Neither symmetry, scalar work cancellation nor the
old zero-`ψ` body-force smoke test is adequate. Later adoption still requires
world-state enrichment and an actual nonlinear trajectory/force comparison.

## Optimization checkpoint: application-specific cost

The development target is the 24 ft, 2500 kg boat and observed SoCal spectrum,
with phase variation retained. Substantially different boat scales, wave-energy
levels and frequency regimes are outside scope. Do not build general backend
selection or extensibility to accommodate them.

On the actual boat's 492-panel flat-water diagnostic boundary, a warmed current
BEM energy assembly costs about 16 ms. Reusing its exact current factor in the
shape force gives agreement with fresh assembly to about 1.3e-13 relative error,
but the force still costs about 445 ms. One profiled force attributes 402 ms to
the 18-lane CUDA coordinate derivative and approximately 44 ms to host geometry
Jacobians (`artifacts/boat-current-force-profile.json`). This diagnostic has a
2 m-deep finite box; it is not an open-ocean simulation timing.

The next force candidate therefore uses Bernoulli pressure and acceleration
potentials, avoiding the full coordinate gradient. Compute the six acceleration
responses only at the current geometry; directional offset fields need the one
prescribed world twist, not six new response columns each. Verify temporal-step
and spatial convergence before adopting this simplification. Operator caching
and adaptive updates remain candidates after profiling that physical stage.

The current nonflat pressure audit now runs on the actual seed-42 geometry.
It found and corrected a vertex-count-dependent polygon fan-center jump, then
removed artificial submerged parent cuts while preserving exact physical facet
areas and first moments. Body triangles fall from 110 to 46. This avoids the
coalescing independent unknowns introduced by those bookkeeping cuts; the
corrected mixed blocks have condition numbers around 43 rather than 1e7.

Pressure integration itself costs less than 1 ms in the bounded audit, with
offset field assembly/solve approximately 10-13/13-17 ms. Fraction reference
geometry remains 1.7-1.8 s per offset, so these are component observations, not
a fast coupled simulation. A four-step directional check exposes a noise floor:
1e-3/5e-4 s passes the weight-scaled numerical budget, but yaw torque still
differs by 5.8% of its response. Frozen controls rejected target-quadrature
refinement for the pre-union jump and gauge correction for the post-union
sensitivity. Coarse parent rates omit
38% of the fine flux norm. The pressure formulation therefore remains a
candidate: physical/state/quadrature convergence, production geometry,
near/far coupling and actual six-DOF trajectory validation are still required.

Subsequent controlled allocation measurements retain the same physical state
and prescribed parent direction while refining the hull or free boundary.
Hull level two uses 1046 panels rather than uniform level two's 4832; current
energy assembly drops from 2.14 s to 0.117 s, with unresolved differences in
sway/roll/yaw and no established pitch/surge plateau. Raw six-axis evidence is
in `boat-refinement-allocation.json`; this is an allocation candidate, not an
accuracy certificate.

That candidate now drives a bounded **coupled** nonlinear six-DOF/free-surface
RK23 pilot (`boat_moving.py`). Current guarded geometry reconstructs at each
stage and refreshes exact clipping when necessary; current base and two offset
operators are still assembled. Two 0.2 s pilots compare 0.05/0.1 s step caps;
the latter accepts two steps, needs only three startup charts and has small
physical endpoint differences over this window. Checked current material-face
containment further reduces sampling cost while reproducing the coupled
endpoint to roundoff. `boat-moving-pilot.json` records total startup-inclusive
times separately from post-startup RHS component costs. These observations do
not establish deep-water, temporal, spatial or long-trajectory qualification.
The finite tank and 15-coefficient surface state must still be replaced or
qualified for the actual open-water problem; finer boundary agreement alone
does not settle them.

An offline reconstruction of the accepted 0.05/0.20 s states identifies a
representation cost before introducing operator caches. Only 128 of 752 body
panels retain their material coordinates, although these cover 66.45% of wet
area. Moving facet-fan centers spread waterline motion into 338 of 415 body
nodes (7.08 cm median material displacement among moving nodes). Source and
hashed reconstructions: `artifacts/audit_boat_material_motion.py` and
`artifacts/boat-material-motion-compatible-chart.{json,npz}`. The next geometry
candidate, evaluated below, is a fixed material hull tessellation, clipping only
its waterline fringe and retaining intended material edges. Physical moments,
closure and pressure-offset continuity are required; no generic remeshing
framework is needed.

Rigid bulk/bulk free-space kernel pairs preserve distances, areas, normals and
paired surface-curl products under common six-DOF motion. This permits exact
reuse of corresponding **raw** V/H/W blocks. It does not make their condensed
energy block invariant: `E=W+H.T@inv(V)@H` remains globally coupled to changing
free/fringe panels. With the current fan mesh, invariant bulk/bed pairs occupy
only 1.63% of raw pairs, so direct pair caching alone has little leverage.
Any new tessellation therefore earns its place through an actual geometry and
pressure comparison before implementing partial assembly or factor updates.

The nonlinear graph HOS reference (`hos_graph.py`) separately passed independent
harmonic-field checks, including direct continuation, grid and product-padding
controls. It is a prospective far-field component, not an implemented near/far
coupling. Incoming analytical waves, outgoing body radiation, interface pressure
continuity and long-time periodic-boundary return still need separate checks.

## Operator integration and material mesh checks

The next independent checkpoint separates finite-domain bias from body
discretization using exactly the same flat-startup wetted triangles. At 608
body panels, the finite 2 m-deep tank raises heave instantaneous added mass by
9.05% relative to the unbounded mirrored P1 formulation. Against Capytaine on
those same triangles, deep P1 diagonal errors range from -11.45% to -4.24%;
production halfspace P0 errors range from -1.69% to +8.45%. Neither sequence
has established a six-axis spatial plateau. Finite-domain and spatial errors
can cancel, making the tank appear closer to Capytaine. These are
infinite-frequency/instantaneous inertia checks, not finite-frequency or
trajectory correlation. Raw matrices and provenance are retained in
`boat-operator-integration.json`.

That actual geometry also exposes inadequate default outer quadrature on
extremely slender cut panels. After diagonal equilibration, the single-layer
matrix's smallest eigenvalue is -0.0858 at order three, +0.0890 at order six,
and +0.1292 at order twelve. Most negative-mode weight lies on thin body
triangles. This is integration-induced negative energy, not merely raw matrix
scaling; no diagonal shift or jitter is applied.

A bounded allocation candidate uses order six on complete target rows with
maximum-edge/height greater than ten and order three elsewhere. Keeping one
rule throughout each target row preserves the closed-boundary mixed-trace
constant cancellation. On the checked flat boat, assembly with precomputed
row classification takes about 0.084 s versus 0.209 s with uniform order six
in one three-sample comparison;
all six added-mass diagonals differ by less than 0.07% from uniform order
twelve. On the nonflat observed-sea snapshot, full nonlinear pressure differs
by at most 1.48 N and 2.56 Nm, with similar h/h2 differentiation sensitivity.
The experiment is `boat_adaptive_quadrature.py`. It keeps all interactions and
physical geometry. Threshold switching, moving trajectories, spatial error
and the exterior domain remain unqualified; moving defaults are unchanged.

Matched coupled six-DOF/nonlinear-surface pilots now compare the candidate
with uniform order twelve at the same 0.2 s endpoint, identical initial state
and timestep controls. Each accepts two steps. Candidate-minus-reference
endpoint differences are 0.0203 mm maximum position, 1.19e-5 rad physical
rotation, 0.000351 m/s velocity, 0.000123 rad/s angular velocity, 0.263 mm
parent height and 0.000286 m2/s parent potential. Startup-inclusive times are
14.64 versus 34.09 s. Summing post-startup RHS costs gives 25.49 versus 98.79
wall/sim seconds. Thus the candidate remains above the user's 20x ceiling
in this run, even before ocean qualification. This is an integration-cost
comparison against uniform order twelve, not a speedup over the earlier
lower-accuracy order-three pilot. Sequential GPU clock/thermal variation is
uncontrolled. The short finite-tank window and coarse surface state still
prevent a sustained-runtime claim.

An exact fixed-material hull overlay now tests the proposed mesh change before
building a fast update path. It clips convex wet parent pieces against dyadic
material triangles, then removes artificial world-parent cuts within each cell.
Full cells keep their original three corners. Shared material/waterline cuts
propagate into neighboring body/free faces exactly, without coordinate welding.
Actual moving states exposed disconnected wet components and concave fringe
loops; components remain separate, while nonstar fringe loops use exact ears
instead of an invalid centroid fan. Holes and ambiguous branches still reject.
Nine geometry tests include the actual canonical boat, advancing incident sea,
six-axis motion and the saved 0.2 s nonlinear state. They check closure,
per-group area/moments, parent inheritance and invariant bulk material triangles.

Level two has 272 body/570 total panels at the reference snapshot. Its pitch
load is about 100 Nm above the old body-two reference; changing uniform outer
order six to twelve changes pitch by only 0.30 Nm. Level three has 848 body/
1194 total panels; its six fine-step loads are [736.05, 322.00, 29962.60,
-285.70, 6.48, 449.35] N/Nm, compared with the uniform-order-twelve old mesh's
[735.29, 337.28, 29961.15, -289.76, -0.64, 432.63]. All use the same prescribed
physical state direction, not independently evolved surface rates. Spatial
adequacy and sway/yaw differentiation sensitivity remain unresolved. Cold exact
clipping costs dominate the 106 s level-three audit; this is no runtime baseline.

Across accepted 0.05/0.20 s states, common whole material triangles retain
86-93% of wet area. They account for only 4.3-4.5% of all raw panel pairs,
limiting a cache of those pairs alone. No cache is implemented, and the
condensed energy remains globally coupled. This reference earns a route toward
stable material coordinates and guarded local operator linearization, not
baseline adoption or a claimed large speedup. Full evidence and source hashes:
`boat-material-mesh.json`; experimental implementation:
`tests/experiments/free_surface/boat_material_boundary.py`.

The first operator-reuse falsifier uses accepted 0.05/0.10 s states to predict
the 0.15 s geometry on the same guarded chart (1070 panels). Its rank-one map
projects nodal displacement in a fixed lumped-area metric; the same map is
evaluated on actual pressure-direction offset geometries. Frozen and secant
operators change fine-step heave/pitch loads by [-373 N, -358 Nm] and
[-85 N, -61 Nm], respectively. Fresh h/h2 sensitivity is below 0.001 N/Nm.
Neither approximation is adopted.

Replaying saved matrices separates operator value from its directional slope.
Correcting the center alone still leaves -68 N heave and -65 Nm pitch error;
correcting the slope alone reduces those to -16 N and +3.8 Nm. The exact-center
case has essentially zero current field residual despite wrong pressure.
A refresh controller must therefore check operator variation along the current
motion/surface direction as well as the current field. These exact-slope controls
use held-out operators and are diagnostics, not a runtime policy.

A one-direction GPU action now differentiates the shared numerical pair formula
and reuses the current V/H factor. It applies all product-rule terms of
`E=W+H.T@inv(V)@H`, including both surface curls, without forming dE or dW.
On the actual state, halving h reduces action/finite-difference discrepancy
from 2.08e-6 to 5.21e-7. Independent canonical-hull rigid/scaling identities
also pass. Uniform-rule ownership is explicit; stale geometry or alternate
integration kernels reject. Matched order-six action medians are 0.387/0.674 s
versus 0.394/0.596 s for two assemblies, so there is no consistent speed win.
This is not yet a differentiated pressure solver and remains experimental.
Matrices, all six loads, timing samples and producer hashes:
`boat-operator-prediction.json`.

Two-direction local models now test whether curvature buys useful reuse on this
specific boat/sea case. Four physical perturbations about .10 s fit a separable
quadratic; two additional corner operators fit the mixed term. Training direction
one uses .15 s geometry, so .15 is not a fully held-out check. At the later .20 s
state the models are reused unchanged, with no new training assemblies.

At .20 s, frozen/linear/separable/full-quadratic heave errors are approximately
-721/-105/-27/-19 N relative to fresh operators on the same prescribed pressure
direction. Adding the mixed term worsens pitch error from -3.3 to -9.0 Nm;
it does not uniformly improve loads. A separable refit using the same six samples
has -45 N heave error. This supports measuring interaction, not a conclusion
that curvature alone explains the remaining error.

Four-sample setup costs 5.619 s; adding the two corners and refitting costs another
4.391 s. Normalized coordinates reach 118 and quadratic features 7004 at .20,
far outside the tiny training steps. No amortized speed or coupled trajectory
accuracy is established. Keep these models experimental. Full six-axis errors,
h-pair sensitivities, cost breakdown and verified producer provenance are in
`boat-local-operator.json`.

The first past-only coupled reuse falsifier trains one geometry direction from
the fresh base/plus/minus operators already required for pressure. Separate
body/free geometry and direction guards, age, extrapolation and ownership apply
to the entire RHS. Candidate failure discards all predicted values and rebuilds
the base and both offsets from the fresh surface rate. No mixed-model derivative
is used. Numerical model-fit failure leaves the successful fresh RHS intact.

This policy fails to save work in matched 50 ms pilots with 10 ms maximum RK23
steps and uniform order-six quadrature: all 16 RHS remain fresh (48 assemblies).
There are five body-geometry, five free-geometry and four body-direction
rejections; initial stencil ownership also prevents the first model from
training. Candidate/fresh totals are 28.42/25.93 s; one pair cannot isolate
overhead from runtime variation. Maximum accepted position difference is
1.03e-12 m, but no reused derivative was accepted, so this only verifies
fallback agreement. The small steps were chosen to give conservative reuse a
chance, not to claim simulation throughput. Geometry limits remain unchanged;
reuse stays opt-in and experimental. Full evidence: `boat-operator-reuse.json`.

A current-anchored, past-only two-direction linear model now preserves the
current center/odd pressure chord and fits transverse response from older fresh
operators. At .15/.20 s, heave errors are -16.1/-87.7 N and pitch errors
-9.0/-38.7 Nm. The training-anchor pressure mismatch is below 1e-4 N/Nm, but
matching an odd operator chord alone does not reproduce individual offsets.
The ordinary same-data LS control loses 5.3% of that current odd chord.

One matched .2 s nonlinear six-DOF/evolving-surface pair with uniform order-six
integration accepts seven reused RHS and cuts assemblies 39 to 18. Total time
is 21.89 s fresh versus 18.16 s cached; a single pair does not establish a timing
distribution or sustained throughput. All guards, fallback, fitting and refresh
audits count. Maximum accepted position-component/orientation differences are
4.95 micrometers/5.81 microradians, while angular-velocity difference reaches
1.031e-4 rad/s, approximately the exploratory component atol. Keep opt-in.

A fresh full RHS at the cached run's exact endpoint isolates remaining policy
error: -26.63 N heave and -11.87 Nm pitch, with no trajectory-state difference.
The new linear model can also differ from the returned fresh anchor RHS, so
this history-dependent policy has no established RK/FSAL consistency or temporal
order. Full results, costs and measured producer provenance:
`boat-past-operator.json`.

The longitudinal-curvature experiment now uses `[a1, a2, .5*a1^2]`, enforcing
both current offsets with a full SVD and fitting the remaining nullspace
coefficient to past records. Observability is checked after elimination; a
parabolic similarity path with independent geometry directions but unobservable
curvature correctly rejects. Anchor total/increment errors are checked with an
additive roundoff allowance, including independently assembled near-invariant
translation controls. This is neither a full Hessian nor a force certificate.

At the .15 s holdout, heave/pitch errors improve to -6.02 N/+0.77 Nm, with no
new training operators. The default-h anchor pressure agrees exactly with
fresh; independent h/2 anchor differences are below 5.1e-5 N/Nm. Thus this
model removes the previous direct anchor/FSAL mismatch. Cache decisions remain
history-dependent, so temporal order is still unproven.

The coupled .2 s curved pilot again uses 18 versus 39 fresh assemblies. The
fresh/cached totals are 24.12/24.21 s in this pair; substantial startup and
per-assembly variation prevents an elapsed speed claim. The same-state endpoint
audit improves heave/pitch errors to -3.53 N/+2.54 Nm, and maximum accepted
angular-velocity difference drops to 8.90e-6 rad/s. Keep opt-in. Evidence:
`boat-curved-operator.json`.

The requested longer all-fresh control fails at its .25 s trial (accepted .2 s):
some current interior material pressure samples become dry at a 1 ms offset.
This is a pressure-derivative sampling limitation, not a rim-immersion or
loss-of-contact claim. No longer cached run is timed against a failed control.
`tests/data/boat-state-250ms.json` records that derived RK23 trial, not a
validated trajectory. Next validate nodal ALE rates on a strictly persistent
wet chart, subtracting actual chart mesh velocity in Bernoulli. Such chart
nodes remain wet while their material labels can slide along the hull; no dry
potential extrapolation is warranted. Require harmonic moving-mesh identities
and actual-boat h/h2 checks before returning to longer timing.

The opt-in charted derivative now follows identical wet nodes through one
explicitly identified base chart, using actual chart velocity in the existing
Bernoulli chain rule. An oblique translating harmonic-flow patch with a sliding
interior node recovers exact pressure. At the actual .25 s failed trial, material
coverage still rejects while nodal h/h2 loads differ by at most .00435 N/Nm;
the maximum facet-normal mesh/rigid velocity mismatch drops from 4.73e-7 to
1.18e-7 m/s. At the covered .20 s state, nodal/material loads agree within
.00714 N/Nm. These are derivative checks, not physical validation.

At initial topology degeneracy, neither offset shares the base nodal chart.
`boat_moving --charted-pressure` explicitly retries a wholly fresh material
stencil on chart failure, recording representation and cost. Dry samples still
reject. Refresh audits record both representations; a mixed comparison is not
labelled operator-only error. Adversarial review prompted these distinctions.

The next fresh run advances past .25 s, then fails at .275 s because the current
wet-facet polygon is not star-shaped about its centroid. The exact triangulator
now retains existing valid fans and partitions concave exceptions into exact
ears, preserving all shared boundary cutpoints. Each ear uses the existing
centroid-fan/chart machinery. An independent concave control verifies exact
area, first moments and boundary edges; the actual derived .275 s state verifies
cold construction, chart replay and closed opposite edge winding. Union-component
counts precede ear splitting; additional pieces are reported separately.
The subsequent run accepts .30 s, then finds multiple disconnected wet loops
on a hull facet at the .3375 s trial. Enable the existing exact positive-loop
support: each component retains its boundary/ownership, without bridging the
dry gap. Holes, opposite winding and ambiguous branches still reject. The
saved trial has the same cold/replay/closed-edge regression coverage. Exact
ear clipping is shared with the material-overlay experiment, with no duplicate
triangulation algorithm. No complete one-second timing has yet been earned.
Evidence and archived measured producers: `boat-ale-pressure.json`.

## Current cost and next exact simplifications

Repeated exact intersection-plane predicates were unnecessary work. A bounded
cache stores three sign bitmasks per intersection triple, retaining no large
Fraction values. Two alternating-order matched pairs at each of the actual
.275/.3375 s states give 1.81x/2.10x geometry speedups, with EVERY saved geometry
and ownership array bitwise identical. Enable it in canonical `construct`;
unrelated reference registries keep their previous default. The cache-off path
includes new counters/sign-comparison overhead, so this is its own matched
comparison, not an exact timing replay of earlier uninstrumented code. Evidence:
`boat-geometry-cache.json`.

Classify a failed fan only after all owner/birth/sign/rim/domain guards pass.
At the .275 s trial, an old .20 s chart passes those guards but loses fan 71.
That earns a local-retriangulation experiment; it does not authorize returning
an inverted triangle. Several historical cold rebuilds may therefore be
avoiding a numerical triangulation issue by repeating physical clipping work.

The .35 s trial leaves the chart at 1 ms, .5 ms and .125 ms. At 15.625 us both
h/h2 stencils complete, but load differences reach 2.423 N/Nm. The bounded fresh
driver now halves h for chart validity, with a minimum 1 ms/64, before explicit
material fallback. This is NOT an accuracy-adaptive derivative. Each immutable
RHS shares its current center field/rates across attempts, but not between
candidate and fresh operators. Cached candidates use their actual training h
and representation; refresh diagnostics record both, protecting the prior
anchor-consistency checks. Evidence: `boat-chart-step.json`.

The next bounded all-fresh run accepts .91990 s in 180.12 s before its elapsed
budget expires: 195.8 wall seconds per accepted simulated second. All 98 RHS
complete; 294 operators consume 110.72 s and 29 cold geometries consume 43.42 s.
This is a timed nonlinear prefix, not a completed one-second control or an
acceptable baseline. Exact geometry caching has not solved the dominant cost.

The larger-budget repeat still fails before one second: it accepts .924931 s
in 169.66 s, then the .924942 s trial has no common chart down to the minimum
h and loses material coverage on fallback. The RK steps shrink sharply near
this event; a cause has not been isolated. This is a distinct failed control,
not a timing improvement or stability result. The exact attempted state is
`tests/data/boat-state-925ms.json`; do not rerun the whole trajectory merely to
rediscover the failure. First test the pressure/triangulation event directly.

The next operator experiment keeps the full current center solve, including
six acceleration columns, and uses its mixed factor only to precondition
single-RHS CURRENT offset operator actions. It can avoid forming full offset
energy matrices while retaining current V/H quadrature/factors and actual-time
loads. It must match offset potentials and all six h/h2 pressures before any
adoption; arithmetic savings are not a timing prediction. Concrete falsifier:
`boat-offset-action.md`.

The motion accuracy targets also earn revisiting cheaper body/quadrature meshes
and temporal controls. Previous comparisons concentrated on instantaneous loads
without a user-owned motion-error budget. Treat 0.1 degree / 1 cm as observed output
convergence gate, never a local integrator atol or a promise of physical-model
accuracy. The current numerical event and the finite-tank/open-water limits
remain unresolved.

The larger-cost avenue remains reduction before geometry approximation. At a
fixed geometry, six rigid modes plus fifteen parent potential modes give only
the homogeneous response space: the observed sea also drives time-dependent
wall/bed Neumann loads. An exact decomposition is `phi=R(Q)*[V,psi]+p(Q,t)`;
the reduced stationary functional is affine-quadratic, with time-dependent
linear and constant terms. Dropping the particular field or the discrete defect
of an analytic incident lifting changes the current equations. The 384
components generally have distinct stratified frequencies, so 48 measured bands
do not imply 96 exact harmonic forcing columns.

Before building a reduced geometry expansion, reconstruct full fields and fine
reactions from homogeneous columns plus the exact current particular, checking
nonlinear fine-flux Bernoulli rates and six pressures at correct offset times.
A small reduced matrix alone cannot replace the fine velocity-squared terms.
Existing P4 restricted-spectrum/rank-six experiments do not qualify this moving
boat affine action. A nonlinear spectral exterior coupled to a local hull BEM
is another larger avenue, but earlier FFT/image coupling failed; domain/interface
flux, work and pressure must pass before another trajectory. Neither avenue has
earned a speed claim.

## Actual-state integration defect: next correction target

Evidence: `boat-contact-integration.json`. Cheaper body-refinement-one runs at
q3/q6 fail the requested motion agreement against the unqualified body-two/q6
prefix. Over the shared 0.92493 s, maximum roll differences are 0.431/0.457
degree and surge differences 2.69/2.30 cm. Neither cheaper run completes one
second within its bounded work budget. The finer reference reaches only about
1.53 degrees roll, so even agreement would not qualify representative motion.

A replacement chart covers the saved 0.925 s pressure trial and allows a fresh
whole-center/stencil audit. The opt-in runtime replacement is still experimental:
the subsequent bounded moving run reaches 0.96446 s in 181.07 s and makes ZERO
replacement attempts. It proves neither runtime repair nor a speed improvement.

That accepted endpoint exposes a tiny free triangle (projected area 1.79e-6 m2)
with a quadrature point 11 micrometers from the contact edge. Its recovered
vertical speed is 599 m/s at q6. A compatible affine harmonic control on the
SAME physical boundary proves numerical error independently of any possible
physical corner singularity. Exact affine vertical speed is 0.4 m/s; the current
continuous graph-flux projection independently recovers a 0.4481 m/s peak.
The solved control peaks are 12.23, 3.88 and 0.842 m/s at q6/q12/q24. Distinguish
the recovery method's smaller consistency error from the much larger BEM defect.

Uniform red subdivision of target integration regions, retaining original
panels, unknowns, self jumps and equations, is rejected as the next default.
Composite q6 with four/sixteen children gives affine peaks 4.93/1.19 m/s;
ordinary q12/q24 is better at the SAME 144/576 sample counts. These diagnostic
costs include capture/compilation effects and are not assembly benchmarks.

The raw pair matrices satisfy `H*u - Vraw*q` to about 1e-15 relative for this
affine field. After `Vs=(Vraw+Vraw.T)/2`, the residual is 1.38e-5 at q6.
The complete energy defect is independently decomposed as
`E*u-b = (Ws*u + H.T*q-b) + H.T*Vs^-1*(H*u-Vs*q)`.
Its residual closes to less than 9e-15 absolute across all five rules. Both
the weak-curl integration and the reciprocity correction contribute; removing
symmetrization does not repair the formulation. Improve the underlying
reciprocal/weak integrals before another long trajectory. No clipping, smoothing
or algebraic affine patch is adopted. The actual-state q24 peak remains about
637 m/s; this alone cannot establish physical validity or its numerical cause.

## Strong trace recovery is not yet qualified

Further frozen-state controls are recorded in `boat-contact-recovery.json`.
Uniform q48/q96 recovers the affine-control peak as 0.44810004/0.44809990 m/s,
matching the independent current-recovery baseline 0.44809989. The actual-state
peak remains 628/630 m/s and its vertical-square integral is 11.7984/11.8005
m4/s2. Thus the proven integration error does not explain the whole actual-state
spike. Three-edge averaged, edge-graded q6 rules also fail to beat ordinary q12
(108 versus 144 target samples); stop this target-rule tuning avenue.

The auxiliary P0 normal trace from `Vs^-1*H*phi` is a second DISCRETE
reconstruction, not a reference solution. At q48 its actual-state vertical-square
integral is 9.05095 m4/s2 versus reaction recovery 11.7984. One conforming free
refinement leaves geometry and the 15+15 evolved state exactly unchanged but
changes these integrals to 8.87818 and 14.1912. Quadrature q24/q48 differences on
that refined mesh are much smaller than this spatial change. Local affine
errors remain (q48 peak 0.6079 versus independent 0.4505), so this is evidence
of unresolved recovery, not a completed spatial-convergence study.

An exact saved-field decomposition separates the sources of the free reaction.
Let `L` map constant face fluxes to integrated nodal loads, and `b` be the
prescribed nonfree load. Then, on free nodes,
`E*phi-b = Lfree*qaux + (Lnonfree*qaux-b) + (E*phi-Lall*qaux)`.
Separate body and outer contributions in the middle term. Their vectors close
to 2.8e-17 absolute; the mass Gram matrix retains cancellation cross terms.
At the original q48 peak the signed vertical contributions are approximately
41.63 (auxiliary free), +1308.24 (body residual), +0.003 (outer residual),
-722.20 (consistency), and -0.004 (tangent), summing to 627.67 m/s. This does not
justify subtracting any component. Shared contact-node tests have body and
free support; their weak reaction is not intrinsically a strong free velocity.

A two-level BODY-TARGETED refinement around the original contact node adds only
28 faces (1088 to 1116); conformity also splits adjacent free faces. At the SAME
original physical probe, the recovered velocity changes from +627.67 m/s to
-5.26 at q48 (+6.30 at q24). The sampled maximum moves about 1.8 mm and becomes
860.66 m/s. The original node's physical body/free support ratio falls from
34.8 to 8.70. Global reaction/auxiliary integrals become 11.7231/9.04669; this is
not a converged repair and no refinement policy earns adoption. Point maxima,
fixed probes, quadrature sensitivity, integrated quantities and all cancelling
terms must be considered together.

Next test the surface rates against derivatives of the SAME restricted
stationary energy, including moving contact/basis terms, before squaring a
strong lift of a weak reaction in another long nonlinear trajectory. Do not
declare the auxiliary flux exact, clip velocities, or project away residuals.
Contact-distance bins in the diagnostic use sampled indicators, not exact
geometric strip clipping; global polynomial integrals have the stronger claim.

The scalar falsifier is now implemented in `boat_surface_action_audit`;
evidence is `boat-surface-action.json`. With `b=B*V+b_inc`, define
`F=rho*(phi.T*E*phi/2-phi.T*b)` at the stationary mixed field. Its derivative
must include `rho*rF.T*dS*psi` in addition to the fixed-trace energy shape
derivative and `-rho*phi.T*db`. Complete h/h2 scalar finite differences re-solve
the mixed field. The next analytic implementation must match these controls
before extending to all directions.
The surface one-form is `rho*psi.T*C*eta_dot`, where
`C=S.T*Mxy*S` depends on the moving exposed-water hole. Consequently
`rho*C.T*psi_dot = -F_eta-U_eta + rho*((psi.T*C_eta[k]*eta_dot)_k-C_dot.T*psi)`.
`C_dot` includes both body and surface motion. Compare these WEAK parent loads
with the existing Bernoulli rate, including incident-load and gravity
derivatives; do not compare only the kinetic energy derivative.

Primary directions are interior parent coefficients 7 and 10, explicitly zero
on every free/wall shared node. Both have contact support on this coarse grid;
there is no claimed pure interior direction away from contact. The initial
parent-0 control changes the external-wall top and forcing domain and is kept
as a boundary-work diagnostic, not primary evidence of the interior mismatch.
For parent 7 the action/Bernoulli weak-load discrepancy is 550.23 N at q12 and
554.49 N at q24. Parent 10 gives 396.70/396.76 N. At q24, halving the height
stencil from 10 to 5 micrometers changes these differences by less than
0.0005 N. Independently integrated gravity derivatives `rho*g*C*eta` agree
with volume-first-moment differences within 0.0003 N. Parent-mass construction
also agrees with the existing fine-mass restriction. These checks establish a
repeatable discrete-action/current-Bernoulli mismatch at this actual state;
they do not establish quadrature convergence, physical accuracy, or replacement
rates. In particular the parent-7 mismatch still changes 4.26 N from q12 to q24.

The order-matched derivative chain now passes the frozen-state controls;
evidence is `boat-surface-gradient.json`. Auxiliary reuse requires an EXPLICIT
matching quadrature order, uniform rule, exact current geometry and winding.
That order reaches every geometric-Jacobian seed. The complete functional
derivative includes both boundary-load and moving-free-trace terms. At q12,
the smaller-stencil analytic derivatives match re-solved scalar differences
within 5e-5 N. A single vertex gradient also contracts to the independently
computed directional energy actions within 8e-10 N on the same geometry.
Actual q3/q6 controls and hull translation/rotation/dilation identities pass too.
This validates the discrete chain; it does not enable new surface rates.

The original 18 scalar-lane vertex pass costs about 40 seconds at q12, so
computing it freshly at every stage is not a credible timing architecture.
Three shared-primal derivative lanes now form the experimental helper default;
the same numeric pair formula and stationary scalar retain single ownership.
Four matched calls on this machine, in order width1/3/3/1, cost
39.69/34.73/36.61/46.88 s. Observed paired ratios are 1.14 and 1.28, with the
second pair the cleaner warm comparison; startup and thermal variation are not
eliminated. All 1680 gradient components agree with the validated reference
within 3.81e-8 N, and its four directional contractions within 5.5e-10 N.
The compiler reports 255 registers for both layouts and a larger per-thread
local frame (504 versus 72 bytes) for width3. These attributes do not establish
actual occupancy, peak VRAM or the cause of the limited improvement. Width1
remains an explicit reference; width6 is tested but not timing-qualified.

All 15 candidate parent weak loads are now constructed from the validated
gradient and same-chart geometry/load/trace derivatives, without additional
energy assemblies. Prescribed incident-flux ports require endpoint virtual
work: the discrete moving-port expression agrees with an independently
integrated perimeter load within 0.015 N at this frozen state. A manufactured
affine-potential test also compares the port expression with a closed-form
corner load. The bottom incident flux remains a forcing port, not a closed bed.

At q3/q6 the geometry/load work takes 3.95/3.64 s for 64 geometry evaluations,
excluding the primal operator and current gradient. A fresh constant-potential
gauge shift changes the corrected potential rate by at most 6.7e-8 m2/s2;
the uncorrected weak load changes by about 898 N. The h/h2 rate difference is
below 3.2e-7 m2/s2, while q3/q6 rates differ by up to 0.0152 m2/s2. These
checks establish discrete consistency, not physical quadrature convergence.
The boundary-data refactor reproduces every archived field output bitwise,
including all six acceleration responses. Reports and source provenance are
retained in `boat-parent-action.json`.

Six-DOF work/pressure consistency must still be resolved before replacement
rates earn adoption in moving runs. The next derivative speed
avenue should remove repeated primal integration more substantially (a bounded
sample-local reverse/analytic scalar derivative); expanding forward batch width
alone has not earned a larger investigation. No new trajectory or full-RHS
speed claim follows from this component optimization.

## Body work and pressure comparison

The frozen-state action now includes all six body equations. At fixed world
twist, its known fluid load is d(F_V)/dt - F_Q + the moving-surface one-form
derivative - U_Q, plus external-port work and the rotational correction
-omega cross F_omega. F_V=-rho B.T phi uses the full nodal body map, including
contact rows. Both B and phi are differentiated. The action and pressure routes
share exactly the same four offset fields and candidate surface rates.

The q12 scalar functional checks of translation and roll derivatives agree
within 0.0004 N or N m. Hydrostatic pose derivatives agree with direct pressure
within 0.0007 N or N m. Independently assembled added masses agree within
2e-11 in their corresponding mixed units; action reciprocity is 8.1e-17
relative. These establish discrete identities and signs, not continuum accuracy.

Dynamic action and pressure loads nevertheless disagree. At q12, h/2,
action minus pressure is approximately [176, -862, -422] N and
[534, 1494, 347] N m, producing a roll-acceleration difference of 0.261 rad/s2.
The mismatch persists from q6 to q12 and under h/h2 changes. Subtracting the
shared -rho B.T phi_dot term isolates it in the geometric terms: that temporal
identity agrees to 1e-10, and the remaining centered-product error is below
0.0005 N or N m. A constant potential gauge changes the complete action load
by less than 0.00008 N or N m. The canonical boat's independent uniform-coflow
test recovers equilibrium hydrostatics with the actual sliding chart velocity.

Evidence and exact source/input ownership are in `boat-body-action.json`.
The discrepancy's spatial convergence remains unresolved; neither formulation
has earned physical qualification. The action equations provide a candidate
consistent with the current discrete energy, not a calibrated pressure fix.
## Experimental action-coupled motion

The complete action equations now run behind `boat_moving --action-coupling`,
retaining matched pressure diagnostics and fresh current/offset operators.
The default symmetric initial position fails the owned mass-transport chart
before advancing. The explicit generic position x=0.123 m, y=0.087 m permits
the diagnostic pilot; it changes the physical initial condition and does not
resolve the default failure.

At q3/body1, 20 ms completes in 32.8 s with seven RHS calls. Increasing the
step ceiling to 50 ms completes 0.2 s in 74.3 s with 13 RHS calls. A one-second
attempt with a 0.2 s ceiling exhausts its 40-RHS bound after accepting 0.6463 s
in 205.1 s, with eight accepted steps and 18 cold chart builds. Neither that
partial ratio nor the short complete runs approach the 20x runtime ceiling.
These small motions and unqualified coarse discretization establish execution
and cost only, not motion accuracy.

Geometry/load derivatives now dominate the RHS. Incident-load derivatives
cancel exactly with the corresponding open-port virtual work, leaving
-rho times the current port covector contracted with the geometry derivative.
Skip incident-wave evaluation only in those geometry stencils; every actual
fluid solve retains full incident forcing. The reference derivative mode
continues to evaluate both terms separately. Four alternating frozen-state
calls take 6.20/4.98/5.18/6.30 s (reference/combined/combined/reference), with
all parent rates agreeing within 9e-16 and six body pose loads within 1.3e-7
N or N m. Separate diagnostic terms are not misrepresented as standalone port
work after cancellation.

A matched complete 20 ms coupled run drops from 32.8 to 28.2 s, an observed
1.16x gain. Accepted states differ by at most 5.7e-10 in their respective units;
angular differences stay below 1.2e-10 degree. This earns the cancellation in
the experimental action path, with single-run timing caveats. Provenance,
controls, partial failures, cost breakdown and memory metrics are retained in
`boat-action-moving.json`.

## Analytic current geometry/action derivatives

The experimental action path now differentiates its owned chart analytically.
One batched owner-intersection solve supplies all 21 world-pose/parent-height
directions. Fan centroids, inherited free traces and body-refinement midpoints
follow the same recipes as the primal geometry. Body loads, projected parent
mass and fluid gravity use exact polynomial derivatives; Cdot contracts them
with the current world twist and parent height rate. The current operator and
vertex gradient are still fresh, as are all four time-offset fluid fields.

On the saved nonlinear state at q3/q6/q12, the complete chain takes 0.13-0.15 s
versus approximately 5 s for combined finite differences. All parent potential
rates agree with h/2 controls within 1.1e-7 m2/s2; body pose loads agree within
0.0007 N or N m, dominated by subtractive cancellation in the old gravity FD.
CPU tests check all 21 directions, a mixed direction, B/C/U derivatives and
two refinement levels. Aliases and relevant incident derivatives remain guarded;
the default symmetric initial chart still rejects.

A strict ownership failure after a cold refresh was traced to rational cold
coordinates differing from floating replay by up to 8.9e-16 m. The action
tangent branch now selects the same guarded replay BEFORE building its primal
operator or derivative. Bitwise ownership checks remain. A later guard at
t=0.921077 s isolated a triangular fan's area-quotient derivative cancellation.
Using the exact triangle-centroid differential (the mean vertex Jacobian)
reduces the incidence residual from 1.3e-12 to 4.4e-15 without changing primal
coordinates or guard thresholds. All 21 h/h2 geometry directions replay there;
their FD errors remain conditioning-dependent diagnostics.

The matched 20 ms full coupled run falls from 28.2 to 8.56 s (observed 3.30x).
The corrected one-second run COMPLETES in 107.7 s, with 79 RHS calls, 16 accepted
steps and 31 cold builds. Roll reaches 2.17 degrees; this remains a short, coarse
diagnostic. Against the old FD trajectory on their shared 0.6463 s interval,
body-position differences stay below 1.3e-6 m and angular differences below 0.000182 degree.
This earns the analytic derivative as the experimental action default;
`--no-action-tangent` retains the reference geometry stencils. The pressure
driver remains unchanged by default.

Cold geometry consumes 39.0 s; current gradients consume about 31.3 s, offset
fields 12.4 s, and time-offset geometry 6.5 s. Process peak RSS is about 595 MB;
the reserved CuPy pool is about 54 MB (not total GPU memory). This is still
107.7x simulated time, far above the 20x ceiling. Evidence and source ownership
are in `boat-action-tangent.json`. Next investigate filtered cold construction
and fewer RHS calls; action-side affine consistency, spatial/motion convergence
and Capytaine correlation remain required before physical adoption.

## Time-integration sensitivity

Higher order alone does not earn a faster baseline. On the identical one-second
q3/body1 action case, loose RK45 completes in 83.8 s with 61 RHS evaluations,
eight accepted steps and two rejected attempts. The RK23 control takes 107.7 s,
79 RHS evaluations, 16 accepted steps and ten rejected attempts. Both use the
same componentwise absolute tolerances, rtol=0.001 and 0.2 s step ceiling.
Despite the lower cost, loose RK45 differs from RK23 by up to 0.346 degree roll
and 0.445 degree yaw, so it fails the requested 0.1 degree comparison.

Tenfold tighter RK45 tolerances complete in 185.1 s with 163 RHS evaluations,
18 accepted steps and nine rejected attempts. This run differs from loose
RK45 by 0.312 degree roll and 0.413 degree yaw, but from RK23 by at most
[0.0347, 0.0364, 0.0324] degree in roll/pitch/yaw and 4.9 mm per translation.
Keep RK23 as the experimental baseline. This is supporting temporal evidence
on one short coarse case, not a certified error bound or physical convergence.

New runs save the actual SciPy accepted-step dense polynomial coefficients;
RK45 is not reconstructed using the legacy RK23 Hermite cubic. Comparisons
reject missing non-RK23 coefficients and require explicit notes for changed
driver code or integration controls. `boat-time-integration.json` retains all
three run summaries, ownership and comparisons. No integrator speedup is adopted.

## Primary context

The exact momentum-tangent path differentiates the current mixed solve:
`E_UU phi_dot_U = (b_dot - E_dot phi - E phi_dot_prescribed)_U`.
It uses current analytic chart/map/body-load derivatives and the material
incident-load derivative, retaining zero body acceleration for the known load
and the existing current added-mass responses. This removes all four offset
fluid solves and geometries. Owned early and difficult-state h/h2 audits agree
with the exact derivative; the difficult finite-difference reference converges
toward it when h halves. Evidence is in `boat-momentum-tangent.json`.

The complete one-second candidate takes 117.9 s, 64 RHS/assemblies and 29 cold
builds. It agrees with the archived RK23 trace within 0.0101 degree and 0.066 mm
on the sampled common grid, but costs more than that historical 107.7 s run.
That historical comparison alone did not earn adoption. A contemporary serial
matched pair now completes the FD control in 212.4 s and exact derivative in
130.8 s, an observed 1.624x improvement. Maximum sampled differences remain
0.0101 degree and 0.066 mm. Runtime variation precludes a guaranteed speed
factor; neither trace establishes full-minute accuracy. Independent review
earns the exact derivative as the experimental analytic action-path default;
`--no-momentum-tangent` retains the h/h2 reference. Issue #30 owns active status.

Variational free-surface discretization with moving boundaries is established;
its Hamiltonian can be nonautonomous when boundaries move:
[Brink, Izsák and van der Vegt, Hamiltonian finite elements](https://link.springer.com/article/10.1007/s10915-017-0416-9).
Near-body BEM coupled to an exterior FEM domain is also established; repeated
remeshing and its energy effects are explicit concerns:
[Wu and Eatock Taylor, nonlinear FEM/BEM](https://doi.org/10.1016/S0029-8018(02)00037-9).
Independent pressure from acceleration potential with current wetted-surface
intersection has been used in nonlinear wave tanks:
[Exact hydrodynamic loads with NURBS and acceleration potential](https://doi.org/10.1016/j.enganabound.2017.09.005).
These support the proposed checks, not a speed or accuracy claim for Topwater.

# Fixed-world surface state for the canonical boat

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

## Primary context

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

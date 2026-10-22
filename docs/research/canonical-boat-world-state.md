# Fixed-world surface state for the canonical boat

Scope: the 24 ft × 9 ft, 2,500 kg boat, water density 1,025 kg/m³,
gravity 9.81 m/s², and selected observed SoCal sea with seed 42. This is an
architecture/force-audit proposal; no new fluid solve or implementation has
qualified it. The failed canoe volume/chart is not the reference geometry.

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

The nonlinear graph HOS reference (`hos_graph.py`) separately passed independent
harmonic-field checks, including direct continuation, grid and product-padding
controls. It is a prospective far-field component, not an implemented near/far
coupling. Incoming analytical waves, outgoing body radiation, interface pressure
continuity and long-time periodic-boundary return still need separate checks.

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

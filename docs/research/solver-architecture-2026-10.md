# Architecture decision memo: nonlinear canoe within 20 wall seconds per simulated second

2026-10-22; independent research, first-principles review. Read-only source/data inspection and primary-literature research; no PDE, GPU, or trajectory jobs. This supersedes the priority ordering in the earlier water-meshing memo. That memo concerned the old body-free P1 pocket; it cannot justify continuing the present whole-domain P4 reference.

## Decision

Do not continue long trajectories with freshly assembled/factored whole-domain P4 tetrahedra. Keep this implementation as a bounded physical/action reference. The most credible next architecture is **current-geometry, matrix-free elliptic solves, a reusable low-order/p-coarse preconditioner, and a regular tensor-product exterior with only a small geometrically complicated canoe region**. Determine whether the near region is better represented by curvilinear prisms or high-order BEM with a bounded actual-canoe comparison. Simultaneously diagnose the sub-millisecond accepted steps; operator acceleration alone cannot recover the target if those steps remain necessary.

This is a numerical architecture change, not a change from fully nonlinear potential flow to linear response. Fresh current geometry and nonlinear free-surface/body coupling remain in the physical residual. Reusing a solver preconditioner is the easiest safe reuse; approximation of the physical operator is a separate, more demanding optimization.

The user requires correct nonlinear pressure/rigid-body work, not specifically the present global canonical mesh chart, P4 everywhere, or complete sparse LU. Those are self-imposed implementation choices. Correctness can also be established for an ALE/acceleration-potential pressure formulation or a boundary-integral formulation. A discrete action derivative passing for the current chart is valuable evidence, not a reason to prohibit a better chart or discretization.

## What the current evidence actually says

Inspected root `p4_tetra.py`, `p4_geometry.py`, `p4_canoe_field.py`, `p4_variational_canoe.py`, `p4_canoe_stage.py`, `p4_coupled_solve.py`, saved RHS profiles and rejected-release records, plus the two previous research memos. Root owns active experiments; no helpers were edited here.

Current actual canoe: 7,076 curved tetrahedra; 35 P4 nodes per tet; 82,432 allocated nodes (82,429 active; three orphan input vertices); 3,168 free-surface coordinates; interior fluid matrix 79,261 rows and 5,501,437 nonzeros. Volume quadrature is an 8x8x8 Duffy tensor rule per tetrahedron. The optimized instrumented profile reports prepare about50.0s and force evaluation11.0s; cold geometry construction is additional20.9s. These are cProfile diagnostics, not accepted-step throughput. The stopped release record reports2752wall seconds for about0.00165s accepted motion, with the accepted time inferred from logs rather than a final history.

Root has already implemented and verified an exact coupled rank-six fluid/body algebra, replacing six auxiliary body-response solves by a joint SPD action. Its latest bounded moved-state probe reports relative current-K change0.00111, four iterations with a cached coupled preconditioner, field difference9.11e-15 from the reference, and nonzero surface potential. Current assembly10.36s remains more expensive than the1.089s coupled solve in that probe. Recommending that as future research would duplicate completed work. It remains useful with every elliptic discretization below.

### Envelope arithmetic

For s expensive RHS evaluations per accepted step, acceptance fraction a, and step dt, an optimistic cost ratio is R=s*C_rhs/(a*dt), before startup, audits and events. Taking s=3 and a=1:

| Accepted dt | RHS budget for R=20 | RHS budget for R=1 |
|---:|---:|---:|
|0.01s|66.7ms|3.33ms|
|0.001s|6.67ms|0.333ms|
|0.0001s|0.667ms|0.0333ms|

These are requirements, not hardware predictions. Against61s/RHS, even the0.01s case asks roughly915x reduction, and sub-millisecond stepping asks much more. A100x assembly improvement alone is insufficient. Faster solve, force integration, mesh representation, and timestep behavior multiply; their speedups must not be added or assumed independent.

The existing quadrature visits7076*512=3,622,912 points. Storing every physical basis gradient would cost7076*512*35*3*8=3.04GB; six symmetric contravariant metric coefficients need only7076*512*6*8=174MB. Reference basis derivatives are common to every cell. This favors compact quadrature data and fused actions, not storing all basis gradients on an8GB device. Even before sum factorization, action B.T*D(q)*B*x avoids materializing35x35 element matrices. The matrix alone is approximately66.3MB with32-bit sparse indices;50 matrix passes move at least3.32GB, ignoring vectors/preconditioner. At an assumed effective bandwidth B GB/s this term is at least3.32/B seconds. B must be measured on this machine; no peak-throughput timing is claimed.

## Ranked candidates and their actual-canoe falsifiers

### 1. Current-action matrix-free solve with reusable low-order/p-coarse setup

**Established, highest probability/ease; a bridge rather than the final mesh.** Apply the same curved P4 quadratic form without global assembly. Reuse the root's coupled rank-six action. Precondition using a lower-order refined mesh, p-coarse operator, or additive-Schwarz/p-multigrid hierarchy; initially reuse existing reference factors only as a controlled comparison, then remove their setup dependency. Current residual and current pressure/action force must use the intended accurate geometry/quadrature. An old factor is not the current operator.

Cheapest falsifier: on the saved actual finite-pose/all-six-velocity canoe, compare matrix-free K*x against the saved assembled K*x for the actual field and several smooth/contact-supported vectors. Compare body velocity, surface rate and scalar-action directional force on that same configuration. Then use one actual nonzero-psi predictor pose from the ongoing root experiment, with the preconditioner held fixed, and count true-residual iterations and full stage cost including metrics, forces, transfers and setup amortization. Stop if repeated current actions cost more than assembly+solve, or preconditioner setup/update dominates. This is not another flat body-free Poisson benchmark.

Accuracy cost: iterative stopping must control output/action errors as well as algebraic residual, especially near poor conditioning. Memory: compact q-data plus vectors; no physical-gradient cache or dense82k inverse. Risks: direct P4 simplex quadrature still expensive, anisotropic cells hurt naive smoothing, GPU gather/scatter may dominate. A matrix-free port without an effective preconditioner is not a solution. [Nonlinear-wave PDC analysis](https://doi.org/10.1002/fld.3873) and [water-wave p-multigrid](https://doi.org/10.1002/fld.5011) establish this direction; neither proves our6DOF contact chart or timings. [libCEED's operator decomposition](https://libceed.org/en/latest/libCEEDapi/) provides the concrete algebra, without requiring adoption of its full software stack.

### 2. Remove timestep pathology; then integrate the stiff linear wave part implicitly/exponentially

**Established methods, potentially as important as operator cost; cause currently unknown.** Do not assume the observed tiny accepted steps are legitimate gravity-wave bandwidth. They could reflect contact-chart conditioning, parasitic high surface modes, force noise, bad tolerances/scaling, invalid trial geometry, or a real short-wave mode. Making an unstable or noisy reference implicit can conceal the defect rather than fix it.

For a pure linear gravity-wave mode, omega^2=g*k*tanh(k*h). RK23's third-order stability polynomial is1+z+z^2/2+z^3/6; on z=i*y its squared magnitude is1-y^4/12+y^6/36, so its imaginary-axis limit is|y|<=sqrt(3). This gives a diagnostic prediction, not a general nonlinear timestep guarantee. Compute the largest relevant eigenvalues of the *actual coupled canoe* linearized surface/contact system and compare predicted stable dt, observed rejection and mode content. A high-frequency mode caused by a nearly singular chart needs a chart fix, not tolerance relaxation.

Cheapest falsifier: freeze the actual calm-release linearized coupled operator and advance its same physical initial perturbation with the current explicit integrator and a linear matrix-exponential/Krylov reference, using the same output norm. This isolates temporal error and contact-mode spectrum without rebuilding nonlinear geometry for every trial. Then run only a handful of fully nonlinear accepted steps using the viable integrator. If a frozen valid linear problem still forces the observed microscopic steps for accurate body outputs, the spatial/contact representation is the first target. If not, diagnose nonlinear chart/force inconsistency before changing physics.

After that: an IMEX/exponential scheme can integrate a rest-state **coupled** linear wave/body generator while retaining the full nonlinear remainder. Unlike a constant-coefficient physical model, L*y+N(y) is an exact decomposition when N(y)=F(y)-L*y. Fast FFT exponentials are straightforward only in the regular far field; the local body/contact part may need a small Krylov or implicit block. Implicit BDF/ALE is established for nonlinear ship waves, but introduces nonlinear iterations and potential phase damping; stability is not accuracy. Use matched trajectory error and accepted-step cost, not dt alone. [Adaptive implicit ship-wave BEM](https://arxiv.org/abs/1203.3384) and [2026 exponential-integrator analysis](https://doi.org/10.1016/j.apor.2026.104951) support investigation; the latter's wave-only perturbative examples do not prove unrestricted6DOF performance.

### 3. Tensor-product/prismatic exterior, graded in depth and localized geometry motion

**Established high-order strategy; largest credible mesh/implementation simplification.** The open-water exterior does not need an arbitrary high-order tetrahedral mesh. Use hexes/quads or extruded triangular prisms with a smooth depth mapping; reserve geometric flexibility near the hull/contact. Matrix-free sum factorization changes tensor-basis application from dense multidimensional contractions to sequences of1D operations, with O(p^4) work in3D rather than dense O(p^6) application. Prism bases retain triangle-by-line structure; arbitrary tetrahedra do not get the same simple tensor kernel.

For one regular brick split into six P4 tetrahedra, shared polynomial-node count is8vertices+19edges*3+18faces*3+6cellinteriors=125. A Q4 hex also has125nodes. The benefit is therefore **not automatically fewer DOFs**. It is one structured cell rather than six, reusable1D basis operations, fewer quadrature geometries, and simpler multigrid/preconditioning. Fully assembling dense125x125 hex matrices could make matters worse; matrix-free is central. These two spaces have different functions/geometry between nodes, so equal node count is not equal accuracy.

Cheapest falsifier: retain the actual canoe near-region and move only an outer regular ring to tensor/prism cells with a common interface. On the actual canoe pose and six loads, compare coupled fluid actions, pressure force/work and dominant wave modes against the reference at matched error. Include a wave packet crossing the interface both directions. Do not begin with a general hex mesher or reimplement a solver framework. If conforming the near patch costs more than it saves, compare an all-prism extrusion or local BEM alternative directly.

Coarsen smooth depth/outer regions by output error, not by distance alone; incident/radiated short waves still need surface resolution. Restrict mesh displacement to a surface/hull neighborhood when a valid mapping permits it, leaving deep-cell metrics exactly fixed. This changes the interior mesh mapping, not the physical boundary. A fixed partition with clearance preserves more cache than a pocket following every infinitesimal motion. Recenter with hysteresis only for sustained travel. Memory is roughly linear in DOFs/q-data, with no sparse LU fill. Caveats: curved contact compatibility, interface resolution and anisotropic preconditioning. [Tensor GPU kernels](https://arxiv.org/abs/1711.00903), [p-multigrid water waves](https://doi.org/10.1002/fld.5011), and [FNPF-SEM2025/2026](https://arxiv.org/html/2506.09435v1) establish the ingredients; the last uses extruded elements and validates waves/fixed structures, not free6DOF canoe contact.

### 4. Fast high-order BEM, potentially only in a local pocket

**Established alternative; should remain a real competitor, not be dismissed because the old P1 experiment lacked accuracy.** Boundary-only representation removes volume DOFs and the global volume-map shape force. Dense all-BEM nevertheless grows O(N_boundary^2) in memory/work per direct interaction and O(N_boundary^3) for generic factorization. At10000unknowns one F64 dense matrix is800MB, so several matrices/factors plus workspaces can consume8GB. Prefer FMM/PFFT/H2 actions and reusable/preconditioned iterative solves when boundary size justifies them; benchmark the crossover against simple batched dense GPU at small local size.

PFFT-QBEM for fully nonlinear moving floating bodies is documented in [Yan's MIT thesis](https://dspace.mit.edu/entities/publication/f539bb9b-d7ea-473c-b33f-e2b3c2adf77e). [Harris et al.2022](https://doi.org/10.1002/fld.5051) uses cubic B-splines and ExaFMM for nonlinear surface-piercing fixed-body waves. [GPU Laplace FMM for evolving bubbles](https://arxiv.org/abs/1905.01341) confirms moving-boundary feasibility, with different geometry/physics validation. This is stronger evidence than generic asymptotic speculation, but none gives a timing prediction for the canoe.

Cheapest falsifier: one current actual-canoe boundary configuration, same finite pose/nonflat trace/body velocity, high-order curved boundary approximation and analytic/accurate near interactions. Compare surface normal velocity, all six pressure reactions and action/work against the reference; measure complete kernel+iterative solve+force costs. Use a second slightly moved/contact-changed configuration to test geometry-update cost. If high-order near-singular accuracy or preconditioning dominates, reject before implementing long-time remeshing. FMM accelerates far interactions, not singular quadrature, body/free-surface junctions, or operator derivatives automatically.

A local BEM plus tensor FEM exterior combines geometric freedom with scalable outer wave propagation. Shared potential and balanced weak flux are appropriate; moving the partition is unnecessary unless clearance/error changes. Existing old P1 partition tests had percent-level differences, so continuity in the equations is not enough to claim negligible numerical reflection. Low-rank/H-matrix acceleration is legitimate but not automatic: [the2024 Capytaine compression paper](https://doi.org/10.3390/en17020372) explicitly presents accuracy/compression results without timing claims because bookkeeping dominated its implementation. It is a frequency-domain primary source, not our time-domain baseline.

### 5. Eliminate most of the exterior volume: nonlinear spectral/DNO far field

**Established family, higher integration risk but potentially large end-to-end benefit.** On a regular free surface, a Fourier/HOS or transformed-field-expansion exterior evolves surface variables and obtains nonlinear normal velocity through FFT-based operations. With N surface points and nonlinear expansion order M, cost can scale like O(M^2*N*logN) for relevant expansions rather than repeatedly solving a full irregular3D volume problem. This is a method family, not an exact cost/timing for all formulations. Expansion order, depth and wave steepness must converge.

The canoe is a moving surface-piercing obstacle with a free-surface hole; a plain periodic FFT DNO cannot simply be applied through it. Use a local actual-hull boundary solver and a rigorously coupled residual/correction exterior. [Spectral-coupled BEM](https://doi.org/10.1016/j.joes.2023.07.003) is a concrete published nonlinear near/far example. Its sequential split has stated coupling approximations and wavelength/transition-width limits. Do not silently adopt a weak-scatterer or one-way incident-wave approximation.

Cheapest falsifier: same actual nonlinear canoe stage, decompose into a regular nonlinear wave background plus body-scattered correction, solve both directions to residual tolerance, and compare the combined actual boundary conditions/forces at two partition widths. Include deliberately nonzero incident pressure and body-generated waves. If the correction is not sufficiently local or iterations erase savings, stop. An error-controlled truncation is an approximation of a nonlinear PDE solve, not an empirical ROM; it still needs refinement evidence. [Wilkening and Vasan](https://arxiv.org/abs/1406.5226) show cancellation/conditioning failures in several naive DNO expansions, supporting stable BIM/TFE alternatives rather than assuming arbitrarily high expansion order is safe.

## Other credible directions, and why they are not first

**Immersed adaptive harmonic-polynomial cells** solve Laplace locally in harmonic bases on a fixed refined Cartesian background, avoiding moving volume meshes. This is a substantive alternative to both FEM and BEM, supported by [Hanssen2023](https://doi.org/10.1002/nme.7115) and [Tong et al.2024](https://doi.org/10.1063/5.0190037). The latter reports3D waves, sloshing and fixed-cylinder diffraction with immersed boundaries/adaptive cells. Attractive if current contact geometry is the bottleneck, but surface velocity recovery/filtering and action-consistent6DOF pressure need fresh verification; its published fixed-body tests do not settle them.

**Unfitted/shifted-boundary high-order waves** similarly avoid fitted-mesh updates. [Visbech et al., August2026preprint](https://arxiv.org/abs/2608.28123) explicitly treats moving nonlinear free surfaces on a Cartesian background with polynomial correction and gradient recovery. Useful current research, not a drop-in fully verified moving-canoe method. Do not replace the PDE with VOF/Navier-Stokes/SPH merely to avoid the chart: those add physical DOFs and generally work against speed here.

**Exact deep-volume condensation/DtN**: below a fixed horizontal plane in constant-depth open water, the Laplace problem can be eliminated by a Fourier/modal map; in a fixed discrete deep subdomain it can be condensed once. This retains nonlinear physics above the plane, unlike linearizing the moving free surface. It could remove deep DOFs and cache exact setup, but the plane operator is nonlocal and must be applied fast; a dense625x625 map may be manageable whereas an arbitrary growing plane is not. This is especially promising when combined with the tensor exterior, rather than a new generic domain-decomposition framework.

**Mixed precision**: use lower precision in preconditioning and compact geometry only after the true F64 current residual/force check passes. Earlier SP BEM tiny-motion derivative failures remain relevant. Ada hardware is not an FP64 tensor-core accelerator; a Hopper paper's timings cannot be transferred to this GPU. The separate GPU research agent owns machine-specific conclusions.

## Adaptive operator updates: established pieces and a possible contribution

**Immediate next falsifier, following the user's latest caching request:** cache/tangent-update the physical operator itself on the actual canoe; do not stop at preconditioner reuse. This is a contained experiment with the current action reference and can precede a mesh rewrite. If successful over the actual motion/tolerance envelope, it changes the cost balance used to rank the larger alternatives above.

### What to approximate

For a conventional physical-coordinate Jacobian F, the reference-domain Laplace metric is D(F)=det(F)*F^-1*F^-T. Cache this small symmetric tensor at quadrature points, or cache the assembled K0 and apply its directional geometry derivative, rather than repeatedly rebuilding every35x35cell matrix. The reference basis never changes while topology/master IDs remain unchanged. At an anchor F0, let A=F0^-1*(F-F0), G0=F0^-1*F0^-T. The first variation is

    deltaD = det(F0)*(trace(A)*G0 - A*G0 - G0*A.T).

This formula uses physical Jacobian columns; transpose consistently for the current implementation's row-Jacobian convention. It is a derivation of the same FEM energy, not an empirical surrogate. A Neumann expansion of (I+A)^-1 supplies higher terms when ||A||<1. That condition is a convergence/admissibility check, **not** a certified pressure/output bound. A constant operator is the zero-order predictor; a tangent operator is first order. Neither preserves full geometry dependence exactly.

Split cells into: unchanged deep/exterior cells cached exactly; smooth mildly changing cells predicted; near-contact/large-deformation cells recomputed accurately. Cache/refresh geometry groups by measured error, not simply by elapsed step count or distance to hull. Do not assume the current harmonic extension makes far/deep cells unchanged; exact caching requires localization of the mapping or an actual equality check. BEM has analogous geometric block updates; rigid hull self-interactions are invariant under a common rigid transformation only while the represented wet patch/normals/geometry relationship is unchanged. Moving waterline clipping invalidates that shortcut where its support changes.

### Force accuracy is the deciding check

If K_T(X)=K0+DK0[X-X0], then kinetic energy0.5*phi.T*K_T*phi is first order in geometry, and its operator error is generally O(||X-X0||^2). But its geometry derivative uses DK0, while the true derivative is DK0+D2K0[X-X0]+..., so the **force error is generally first order**. A tiny operator/energy error can coexist with an unacceptable pressure force. Derive both field and force from this same cached action, keep current boundary velocity action/body loads/gravity exact, and evaluate their chart pullback consistently. A second-order energy/metric predictor normally gives second-order force error, at greater setup/storage cost. A frozen/tangent force cannot be justified by a scalar energy match alone.

Changing the anchor changes the approximate action. Audit the old/new approximate kinetic energy and its derivatives at the SAME current physical state, not at two adjacent trajectory times. The jump0.5*phi.T*(K_new-K_old)*phi is a diagnostic; a constant offset cannot fix state-dependent velocity/force changes. Do not enforce conservation by projecting a trajectory onto an energy shell. Anchor switches need an error budget and restart of affected temporal history if the derivative approximation changes beyond that budget.

### Error and refresh policy

Use separate controls for (a) fluid algebraic error, (b) geometry/operator approximation, (c) shape-force error and (d) time integration. Small body displacement or a converged residual of K_T does not control the true residual of K. Useful inexpensive indicators are the local relative deformation ||A||, a neglected-next-term estimate, predicted change of each cell's field energy, wet-support/topology change, and cached-force directional error. All are indicators until audited.

There is a rigorous route to a restricted *field-energy* bound: if at every quadrature point -epsilon*D_T <= D-D_T <= epsilon*D_T as quadratic forms, with epsilon<1 and positive D_T, then the same relative form bound holds globally for the fixed-basis stiffness. Consequently ||x-x_T||_(K_T) <= epsilon/(1-epsilon)*||x_T||_(K_T) for identical loads. This does **not** bound the geometry derivative, a changed boundary load, contact mapping, or full6DOF pressure. Computing current small3x3metrics for this audit may still be much cheaper than assembling/factoring; measure it. No cheap full bound is claimed without accounting for those omitted terms.

For a linear output l.T*x, exact current adjoint z gives output error z.T*r with r=b-K*x_T. An approximate adjoint produces an estimate requiring its own error control. Raw residual norm, force-coordinate residual norm, or a few sampled rows is not a certified output bound. For nonlinear pressure/action outputs use directional differences and fresh same-state reference solves until an estimator's reliability is established. Refresh affected pieces and retry the same physical stage on failure; reduce dt only for actual temporal error once operator accuracy is restored.

### Bounded actual-canoe experiment

Use the saved calm finite-pose/all-six-velocity seed and the nonflat high-contact seed. For each, test (1) coupled body rotation/translation moving the wet contact, (2) a smooth propagating surface displacement, and (3) a contact-supported height mode. At two or three displacement magnitudes spanning the *observed accepted-step geometry changes*, compare fresh K/shape force, frozen K, tangent K and, only if necessary, quadratic prediction. Include actual nonzero psi rather than only unit added-mass loads. Report operator action error, all six body accelerations, surface potential/elevation rates, same-state scalar-action work, positivity, sampled Jacobians, and complete update+solve+force time. Root's reviewer is already preparing the tangent-action audit; do not duplicate it.

Repeat the comparison with both the current spatial reference and a cheaper candidate mesh, and temporal error targets spanning the intended loose/tight runs. Operator error should remain below a stated fraction of each run's independently established spatial/temporal output error; there is no universal magic displacement threshold or inherited numerical tolerance. Count total refresh/audit cost per accepted physical time, including false-positive and missed refreshes. Kill a predictor if field error is small but force error fails, the safe anchor radius is too small to amortize setup, or tighter time tolerances expose refresh noise. Only after this bounded result earns its cost should it be inserted into a short nonlinear trajectory.

### What literature does and does not establish

Not novel by themselves: lagged preconditioners, residual-based refresh, recycled Krylov spaces, geometry Taylor/shape derivatives, fast local boundary updates and adaptive grids. [Local-geometry fast inverse updates](https://arxiv.org/abs/1706.01414) already use low-rank perturbations and an existing inverse; global free-surface motion does not necessarily have low rank. [Water-wave p-multigrid](https://doi.org/10.1002/fld.5011) already studies reused initial coarse/smoother data. Therefore do not claim that adaptive operator reuse or a near/far decomposition is a new invention.

[Fazioli and Nicholls2010](https://doi.org/10.1016/j.jcp.2009.10.021) explicitly computes the shape variation of the nonlinear water-wave DNO via stable transformed-field expansions. Together with the [DNO-method comparison](https://arxiv.org/abs/1406.5226), it establishes that geometry expansions are credible published physical-operator approximations and that naive recursions can be numerically unstable. It is not a ready six-DOF moving-contact refresh controller. Local Taylor updating of our cell metric is the FEM analogue derived above; no claim of a new Taylor method is warranted.

A potentially useful project-specific contribution is **one measured error/work controller for refreshing variational operator pieces**: retain a physical action approximation H_tilde(q,p), derive velocity and force from that same approximation, and allocate separate tolerances to current-residual solve error, stale geometry/metric error and temporal error. Refresh pieces only when estimated body/surface output or pressure-work error exceeds its budget, using occasional full current actions at the SAME state as audits. Never interpolate an operator and independently freeze its force derivative. Reject/retry a stage on audit failure; remeshing/contact-topology changes invalidate affected history. A differentiable local Taylor/metric predictor can be tested without training a surrogate.

This policy might combine known tools in a useful way for this canoe, but the search here does not establish novelty. It should be sold as an engineering hypothesis, not a new numerical method. Cheapest falsifier: along two actual-canoe directions (finite body rotation with moving waterline, and a contact-supported surface mode), compare prediction errors for K*phi and the scalar-action force over increasing displacement, with eps refinement; include a topology transition. If cheap indicators miss force/work error or audit/refresh cost cancels the savings, discard the policy. First make the exact current action cheap enough to serve as an audit.

## Next decision, not another indefinite survey

The exact 7,076-cell matrix-free stage and matching force audit now pass; its
warm full-stage cost and frozen contact spectrum are recorded below. The next
decision uses **one actual adaptive-step error audit**, with the same full
nonlinear equations and existing tolerances. Identify the physical state blocks
and contact coordinates responsible for the error norm before choosing implicit
integration, revised scaling or contact-mesh changes. The fastest frozen mode
alone cannot choose among them. Tensor/prism exterior and high-order BEM remain
substantive alternatives if the resulting stage/timestep budget is unreachable;
compare bounded actual-canoe probes before committing to another implementation.
Keep the P4 reference to falsify candidates, not to justify another unusable long
run. Do not build a backend registry.

Promotion requires actual pressure/action consistency, all six body channels, moving wet contact, outgoing-wave behavior, and failure-event handling. It also requires complete accepted-time wall ratio<=20 on this machine, with error/refinement evidence. Rim immersion and complete loss of contact remain stop events. Capytaine frequency responses are external references transformed/compared in the time-domain Topwater workflow; they do not replace nonlinear physics. A fast static operator, an energy-conserving but wrong trajectory, or a huge implicit step is not the outcome.

## Circuit-simulator mechanisms worth copying

Addition following the user's SPICE analogy. There is no fundamental1D-versus3D obstacle: a FEM element's local energy/flux contribution can be treated like a device stamp. The useful ideas are local change detection, cheap tangent prediction, selective reload, convergence verification and rejected-trial rollback. Product rankings and commercial speed claims are not evidence for canoe throughput.

### Actual source behavior

* **ngspice diode bypass:** `cdhat=i_old+g_old*delta_v`; reuse old current/conductance only when both voltage change and predicted-current change satisfy absolute-plus-relative gates, outside predictor initialization. Bypass jumps to the ordinary local stamping code. It does not evaluate the fresh nonlinear current to prove the defect is small. See [dioload.c lines142–161](https://github.com/ngspice/ngspice/blob/master/src/spicelib/devices/dio/dioload.c#L142-L161), [limiting lines168–179](https://github.com/ngspice/ngspice/blob/master/src/spicelib/devices/dio/dioload.c#L168-L179), and [nonconvergence/stamping lines369–395](https://github.com/ngspice/ngspice/blob/master/src/spicelib/devices/dio/dioload.c#L369-L395). MOS1 similarly gates several terminal changes and tangent-current predictions: [mos1load.c lines232–326](https://github.com/ngspice/ngspice/blob/master/src/spicelib/devices/mos1/mos1load.c#L232-L326).
* **Local evaluation reuse is distinct from global factor reuse:** this ngspice snapshot clears the global matrix/RHS and loads device contributions in [cktload.c lines49–70](https://github.com/ngspice/ngspice/blob/master/src/spicelib/analysis/cktload.c#L49-L70); its iteration loop reloads and numerically factors the matrix in [niiter.c lines67–146](https://github.com/ngspice/ngspice/blob/master/src/maths/ni/niiter.c#L67-L146). That implementation is evidence for device bypass, not a ceiling on what another simulator can reuse.
* **Xyce/NOX modified Newton:** Xyce exposes the Modified-Newton direction in [ParameterSet.C lines438–471](https://github.com/Xyce/Xyce/blob/master/src/NonlinearSolverPKG/N_NLS_NOX_ParameterSet.C#L438-L471). NOX's inspected implementation evaluates current F, refreshes the Jacobian by age, and flags refresh after a rescued failed linear solve: [ModifiedNewton.C lines55–96](https://github.com/trilinos/Trilinos/blob/master/packages/nox/src/NOX_Direction_ModifiedNewton.C#L55-L96). That file is guarded by WITH_PRERELEASE, so availability depends on the build; this is not a claim about Xyce's default algorithm. Xyce separately tracks valid Jacobian and solver factors in [Group.C lines313–320](https://github.com/Xyce/Xyce/blob/master/src/NonlinearSolverPKG/N_NLS_NOX_Group.C#L313-L320).
* **Convergence is more than small updates:** Xyce checks device convergence before normal residual/update acceptance and corrects update size for line-search damping: [XyceTests.C lines239–258](https://github.com/Xyce/Xyce/blob/master/src/NonlinearSolverPKG/N_NLS_NOX_XyceTests.C#L239-L258), [lines363–402](https://github.com/Xyce/Xyce/blob/master/src/NonlinearSolverPKG/N_NLS_NOX_XyceTests.C#L363-L402). Its diode records when limiting changed the input voltage: [N_DEV_Diode.C lines1056–1125](https://github.com/Xyce/Xyce/blob/master/src/DeviceModelPKG/OpenModels/N_DEV_Diode.C#L1056-L1125). A limited trial is not automatically a converged physical state.
* **State integration and timestep rejection are separate:** ngspice integrates charge histories through method/order coefficients, not just a frozen capacitance times voltage: [niinteg.c lines21–70](https://github.com/ngspice/ngspice/blob/master/src/maths/ni/niinteg.c#L21-L70). Charge-based divided differences control truncation error in [cktterr.c lines34–70](https://github.com/ngspice/ngspice/blob/master/src/spicelib/analysis/cktterr.c#L34-L70). Transient nonlinear failure rolls back/reduces the step, whereas truncation acceptance is a subsequent test: [dctran.c lines710–809](https://github.com/ngspice/ngspice/blob/master/src/spicelib/analysis/dctran.c#L710-L809). Do not copy a particular reduction factor or tolerance as a universal numerical rule.

These are inspected public source snapshots on2026-10-22, not a claim that all releases/devices have identical gates. The linked master branches can move.

### Broader simulator evidence

The [SIMetrix matrix-solver documentation](https://help.simetrix.co.uk/8.0/simetrix/mergedProjects/simulator_reference/topics/simref_convergence_accuracyandperformance_matrixsolver.htm) explains that efficient ordering/factorization and frequent reordering as timestep regime changes matter. It does not publish a device-bypass refresh formula to transplant. Its [convergence documentation](https://help.simetrix.co.uk/9.1/simetrix/user_manual/topics/analysismodes_convergence.htm) also distinguishes numerical-precision options from changes that physically modify a circuit. No undocumented commercial algorithm is inferred.

[PrimeSim XA's primary datasheet](https://www.synopsys.com/content/dam/synopsys/implementation%26signoff/datasheets/primesim-xa-ds.pdf) documents intelligent partitioning, multi-rate integration and dynamic multilevel synchronization. [Cadence's engineering response](https://community.cadence.com/cadence_technology_forums/f/custom-ic-design/65379/strategy-for-accelerating-full-chip-ams-simulation-with-multi-rate-clocks-1khz-100mhz) similarly identifies Spectre FX partitioning/multi-rate for disparate clocks. These support event/activity-based local work and synchronization, not undocumented error guarantees or claimed speed factors for water.

[Engelhardt's LTspice technical article](https://www.analog.com/en/resources/technical-articles/spice-differentiation.html) emphasizes smooth device values/derivatives for Newton convergence, runtime-generated sparse-solve code, and the accuracy consequences of integration choices. Its performance and product-superiority claims are author/vendor claims, not canoe measurements. The concrete transfers are specialized kernels for fixed structure, smooth current residuals, and comparing numerical damping/phase at matched error. [Qorvo's QSPICE account](https://www.qorvo.com/about/news-events/blogs/mike-engelhardt-the-mind-behind-qorvos-qspice-simulator) identifies redesigned timestep control, but does not publish a detailed refresh algorithm there; it is insufficient to reproduce that algorithm. Public methods from several simulators inform the comparison; ngspice is one inspectable implementation rather than a performance ceiling.

### Concrete transfer to this canoe

For each element or small contiguous patch, keep an anchor geometry/metric, local action and its geometry derivative. Preserve exact constant-potential nullspace. Predict its current flux contribution with the old/tangent metric applied to the **current** potential; use the current field in error weighting. Small potential change is not required to reuse geometry-independent cached data, but it affects whether that reuse is accurate for the current load.

Use two gates analogous to terminal-voltage and predicted-current tests:

1. A local dimensionless geometric change/quality gate, such as ||F0^-1*deltaF|| and positive sampled Jacobians/contact validity. Any wet-support/topology change forces the affected local refresh. This gate protects the expansion domain, not output accuracy.
2. A predicted flux/action-force defect gate. For a frozen metric use the first omitted tangent term; for a tangent use the next term or an audited curvature estimate. Compare each *physical-unit* contribution against `absolute_budget + relative_budget*max(old_scale,predicted_scale)`. Include a nonzero absolute scale near calm water; a relative-only criterion there is pathological. Use separate force, torque, surface-height and potential/work budgets; do not compare their raw mixed-unit norm.

Per-cell numerical tolerances cannot be copied from SPICE RELTOL. Allocate local absolute budgets with weights summing to the stage budget, or verify the globally assembled defect; otherwise many individually small missed stamps can sum to a large error. Body-output adjoint weights can help prioritize cells but remain estimates until their residual/inverse stability is controlled. Keep the action-derived force gate described above: an accurate tangent flux does not establish accurate shape force.

Failing cells are recomputed and re-stamped; for a persistent assembled matrix add the local difference `K_new-K_cached`, rather than clearing it and accidentally omitting cached contributions. In matrix-free form select current versus cached q-data per patch. Updating numerical entries generally invalidates an exact global factor, even with unchanged sparsity; a reused factor can instead precondition a solve checked against the selected current/approximated operator. Count screening, scatter, correction and audits, not only avoided element evaluation.

The safe exact analogue of modified Newton keeps evaluating the intended current residual while reusing its derivative/factor as solver machinery. Freezing FEM K itself changes the physical fluid equation, so it is closer to **device bypass**, requiring controlled approximation/error audits, than to harmless Jacobian reuse. Conversely, fresh small3x3metric evaluation for a defect check can be much cheaper than fresh global assembly/LU. This is an opportunity, not a reason to reject caching because the fluid is3D.

Retry the same stage with refreshed cells/setup when the operator/force audit fails. Only shrink the physical timestep for genuine temporal error, unresolved events, or invalid trial geometry after appropriate numerical recovery. Never accept a clipped/limited geometry as the intended state; flooding and loss of contact remain physical stop events. Hold cache anchors stable within an embedded RK trial so error estimates do not primarily measure anchor-switch noise. On rejection restore accepted-state history; a newly computed exact anchor may be retained only with correct identity and explicit handling of changed numerical history.

The immediate falsifier remains one actual-canoe direction at two displacement magnitudes: compare old, tangent and fresh operators **and their action forces at the same physical state**, recording which local gates fired, missed failures, and complete cost. Start with full same-state audits to calibrate gates, then reduce audit frequency only after evidence. A later activity-based near/far multi-rate scheme must preserve the globally coupled elliptic/added-mass solve; circuit partitions with long inactive intervals do not prove that distant fluid can simply stop responding.

## Exact full nonlinear matrix-free stage, 2026-10-22

The complete nonflat canoe seed now passes an exact-current matrix-free RHS
comparison with the independently saved direct-solve action audit. All 3,168
surface modes, finite body pose, six nonzero body velocities, physical boundary
maps and nonlinear shape/chart terms are retained. The old calm fluid/body
factorization is only a preconditioner. No current global volume matrix,
current fluid LU or six auxiliary body-response fields are built.

`artifacts/p4-nonlinear-exact-matrixfree-rhs.json` records five CG iterations,
field relative difference `1.61e-14`, full canonical-rate relative difference
`4.22e-14`, and maximum body momentum-rate difference `2.50e-12` in the
corresponding SI force/torque components. Weak-normal relative residual is
`4.72e-13`; fluid work mismatch is `5.13e-16 J`. Exact kinetic forces use the
same current geometry and quadrature as the matrix action. Persistent GPU
cache arrays occupy 471,097,296 bytes; this excludes allocator reserve and
temporary arrays.

This first cold run spent 2.53 s on current boundary geometry, 44.24 s on GPU
cache construction including cold setup/compilation, 5.27 s on the coupled
field, 3.84 s on kinetic energy/shape force and 6.42 s on remaining action
terms. Geometry startup and old factorization cost another 42.97/30.59 s.
These are one-state diagnostic timings, not trajectory throughput or a warm
speedup claim. A warm full-stage profile is required before selecting the next
kernel change. Accepted timestep and physical error remain decisive.

The subsequent warm diagnostic is preserved in
`artifacts/p4-nonlinear-exact-matrixfree-warm-rhs.{json,npz,prof}`. Its unprofiled
warm full RHS is 7.94 s, including fresh current geometry, metric and all force
terms. The final profiled repetition is 9.00 s; these two observations are not
a throughput distribution. All three repetitions retain full-rate relative
agreement within `5.62e-14`. Warm phase costs are 1.74 s boundary geometry,
0.556 s exact metric cache, 1.10 s current coupled field, 0.651 s kinetic shape
work and 3.86 s remaining action. The profile identifies unnecessary volume
extensions/transposes in boundary-only terms as a concrete exact simplification.

With three RHS calls per accepted RK step, this measured cost would require
`h >= 3*7.94/20 = 1.19 s` just to reach 20x real time, before rejected trials.
At a hypothetical `h=0.01 s` the estimate is 2,382x. Neither is a measured
accepted timestep. The exact metric cache is now only about seven percent of
the unprofiled stage; improving it alone cannot make this discretization usable.
Boundary work, actual temporal error/scaling and the frozen spectral diagnostic
are therefore the next discriminators before another long trajectory.

The boundary-only transport simplification subsequently passed a full actual
force audit against immutable `8a91e9f` action code and the independently saved
stage. Body momentum rates are bitwise identical; the maximum surface momentum
rate change is `6.94e-18`. The new full momentum calculation differs from the
independent saved stage by at most `7.11e-15`. Both use the same saved physical
field and exact GPU kinetic tuple, so this audits the force path without another
fluid factorization; the supplied configuration-rate half is identical by
construction. One paired action measurement is 2.42 s old / 1.33 s new, not a
throughput distribution. Evidence: `artifacts/p4-boundary-transport-full-action`.

The bounded frozen-body spectrum now converges four largest modes in 41 Schur
actions and 4.34 s eigen work, after 79.54 s setup. Reusing the hash-bound calm
`K_II` and assembling only 1,459 surface-touching cells avoids another whole
volume assembly. The fastest mode has `omega=805.58 rad/s`, period 7.80 ms and
the **frozen linear** RK23 imaginary-axis ceiling 2.150 ms. It has 99.9993% of
its consistent-mass norm on contact faces, with three faces containing 90%.
Those faces have horizontal tangent condition numbers up to 136.24. The four
verified eigen residuals are at most `4.63e-13`; an independent archive audit
confirms mass orthonormality within `4.45e-16` and unit positive face-mass sums.

These are actual spatially localized high-frequency modes, but the fastest
mode's initial actual elevation-rate projection is only `2.22e-6` of its total
mass norm. Their existence is neither the full free-body linearization nor
proof that they caused the observed adaptive steps. Do not delete them solely
from this diagnostic. Actual embedded error, physical-output sensitivity,
contact mesh quality and an implicit comparison are the next discriminators.
The first setup-budget-only spectrum attempt remains preserved as incomplete,
with zero eigen actions. Fresh complete evidence and reusable surface blocks:
`artifacts/calm-canoe-frozen-spectrum-surface-blocks.{json,npz,blocks.npz}`.

The independent 80-cell actual nonlinear predictor audit also falsifies an
operator-only reuse gate: halving the predictor interval from `1e-4` to
`5e-5 s` reduces tangent-action error from `1.13e-5` to `2.84e-6`, but full
shape-force error only falls from 0.891% to 0.446%. Both physical Jacobians and
tangent tensors remain positive. The tangent energy has its matching base
shape derivative; it is a separate approximation and is not used in the exact
full-stage result. See `artifacts/metric-shape-tangent-cohort.json`.

## MATLAB's established solvers: transferable architecture rather than a language migration

The [MATLAB ODE Suite paper](https://www.mathworks.com/help/pdf_doc/otherdocs/ode_suite.pdf), Shampine/Reichelt1997, documents simplified Newton/chord iteration, quasi-constant-step multistep organization and Jacobian reuse. Its performance comparisons count function evaluations, derivative evaluations, LU factorizations and solves separately; fewer Jacobians alone do not establish less total work. Reuse of derivative/setup data is standard numerical practice, not something Topwater should reinvent under a new name.

Current official algorithm choices are concrete: [ode15s](https://www.mathworks.com/help/matlab/ref/ode15s.html) uses variable-step/order NDF1–5 with optional BDF; [ode23t](https://www.mathworks.com/help/matlab/ref/ode23t.html) uses trapezoidal integration with a free interpolant and is appropriate for moderately stiff problems when numerical damping is undesirable. The latter's lack of linear-oscillator damping does not guarantee phase accuracy or nonlinear energy conservation at large steps. [ode15i](https://www.mathworks.com/help/matlab/ref/ode15i.html) accepts a fully implicit residual F(t,y,y_dot)=0 and requires consistent initial values/derivatives; supplying Jacobian information is important. These are candidates after identifying stiffness, not reasons to call BDF universally best for weakly damped waves.

MathWorks explicitly notes that retaining a [mass matrix](https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html) avoids an inconvenient and potentially expensive explicit inverse transformation. The analogous bigger choice here is **reduced ODE versus unreduced coupled DAE**:

* Current reduced ODE evaluates F(q,p) by building/eliminating the fluid field each call. An implicit integrator wrapped around that same black box may make several equally expensive F calls per Newton step. Numerically differencing its6,348-state Jacobian would be especially costly and its eliminated dependence is nonlocal.
* An unreduced residual retains fluid phi as an algebraic unknown alongside body/surface differential variables. Its fluid row is `K_II(q)*phi_I + K_IF(q)*psi - B_I(q)*V=0`; the remaining kinematic, momentum and pressure/work equations must come from the same full action/physical chart. Local element fluxes/forces can be evaluated directly, while modified Newton/JFNK uses cached local tangents and a coupled block preconditioner. Fresh residual evaluation no longer requires a fully converged inner field solve at every call. This changes the numerical organization, not the nonlinear equations.

The DAE is not automatically cheaper: it has more unknowns, still requires a good elliptic/body preconditioner, and nonlinear residual/force evaluation must become cheap. The current harmonic geometry extension also induces nonlocal derivatives; apply its forward/adjoint maps rather than asserting the whole Jacobian is sparse. Establish index/constraint regularity and consistent initialization for the chosen unreduced variables before claiming index-1 support. Keep surface mass/one-form operators and small rigid-body mass in their natural blocks rather than forming dense inverse products. This is the closest direct transfer of circuit simulation's differential/algebraic organization.

**Check scaling before blaming physical stiffness.** Root's current release passes one scalar `atol=1e-10` to positions, rotation coordinates, body momenta, surface heights and surface canonical momenta. These have different units, and pi=rho*C_eta(q).T*psi depends on local mesh mass/area. A uniform canonical-coordinate tolerance is not a specified physical error budget; it can over-resolve some combinations and under-resolve others. [odeset](https://www.mathworks.com/help/matlab/ref/odeset.html) explicitly supports componentwise absolute tolerances, error scaling and sparse Jacobian/mass information. Do not equate its exact acceptance formula with SciPy's different norm/formula.

The cheapest decision check needs no new solver implementation: on an actual embedded error estimate, report which state blocks dominate the existing norm and map its variation into CG translation, body angle/velocity, surface elevation and physical potential. Use constant physically named scales or vector tolerances first, then verify decoded physical errors under refinement. A more sophisticated mass/energy norm is justified only if diagonal scaling fails. Do not simply loosen every tolerance or drop contact modes to obtain larger dt.

The next bounded comparison, once cached actions are affordable, should put explicit RK23 and an established trapezoid/DAE or NDF/BDF implementation on the **same actual canoe initial state and physical error targets**. Count accepted/rejected steps, residual evaluations, local operator/force refreshes, Jacobian setups, factorizations/preconditioner refreshes and solves. An implicit method taking five cheap residual iterations can beat three complete fluid eliminations, but five expensive reduced-ODE calls can lose badly. Compare pressure/body phase and amplitude as well as volume/work. No MATLAB port or bespoke general DAE solver is proposed by this memo.

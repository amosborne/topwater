# Numerical contract and verification gates

This is the decision record for the first numerical implementation. It specifies
work to verify; it does not certify an operator, a wave solver or hydrodynamic
correlation that has not been implemented. The canonical canoe and loading in
`tests/reference_hull.py` remain the one shared hull reference.

## Mathematical signs and physical ownership

Use `u = grad(phi)` and the positive free-space Rankine kernel

```text
G(x,y) = 1 / (4*pi*|x-y|),       Laplacian(G) = -delta
grad_y G =  (x-y)/(4*pi*|x-y|^3)
grad_x G = -(x-y)/(4*pi*|x-y|^3)
q = grad(phi) dot n
```

`n` is always outward from the **fluid domain**. At a wetted hull it is the
negative of the stored body-outward normal. At a graph free surface it points
into the air. Outer artificial boundaries point out of the computational water.
Geometry winding is retained; do not destructively reverse the physical hull.

Define `S q = integral(G*q)` and `D phi = integral(dG/dn_y * phi)`. Green's
identity gives `phi = S q - D phi` inside the fluid. The fluid-side trace at the
interior of a flat panel is

```text
H phi = S q,       H = 0.5*I + K,
K_ij = principal-value integral_panel_j(dG(x_i,y)/dn_y) dS_y.
```

Here `K` excludes the jump term; this definition never changes with an API
convenience. `S` has units of length and `K` is dimensionless. The source-normal
double layer is different from the target-normal adjoint double layer.
The [Bempp representation](https://bempp.com/handbook/theory/greens_representation.html)
and [operator definitions](https://bempp.com/handbook/api/boundary_operators.html)
provide independent references for these conventions; Topwater uses collocation,
not Bempp's Galerkin discretization or implementation.

For a bounded closed fluid domain, `K @ ones = -0.5*ones`. For an exterior
obstacle, with fluid normals into the body, it is `+0.5*ones`; exterior solutions
must decay at infinity. Constant potential is not a decaying exterior solution.
Zero total flux is required for a bounded all-Neumann domain, not separately on
each hull patch and not generally on an exterior obstacle with flux at infinity.
An artificially closed canoe can be a manufactured *interior mathematical domain*;
this says nothing about the physical wetted hull or cap pressure.

Hull impermeability uses world velocities:
`q = (V_CG + omega_world cross (x-X_CG)) dot n`, with
`omega_world = R @ omega_body`. Eulerian Bernoulli pressure is

```text
p_gauge = rho * (C(t) - phi_t - |u|^2/2 - g*z)
dF_body = p_gauge * n_fluid * dS
dM_CG = (x-X_CG) cross dF_body.
```

Atmospheric pressure is zero gauge and gravity points down. Total water pressure
is integrated once; no added buoyancy, fixed added-mass load or separate incident
force is placed on top. Changing `phi` by `a(t)` requires changing `C` by `a'(t)`.
Future atmospheric free-surface evolution fixes `C=0` consistently with its
potential history; repeatedly subtracting a potential mean would change pressure.
An area-weighted zero-mean gauge is only a convention for bounded pure-Neumann
verification and must retain its separate constraint and compatibility residuals.

## First operator and boundary solve

Implementation status: `HullMesh` owns validated physical arrays;
`TrianglePanels` caches O(N) numerical geometry; one fused, serial Numba kernel
fills final S/K buffers directly. `LaplaceOperators` owns those arrays and their
measured diagnostics. Independent quadrature remains a separate verification
method, not a second analytic backend. Boundary factorization is the next layer:
it must own its operators, boundary assignments and scaling, reuse factors for
unchanged systems, and avoid keeping a redundant unfactored mixed matrix just
for residuals. Geometry changes require rebuilding the system. No global cache
or extensible backend framework is planned.

Choose a **direct Rankine-kernel BEM**, whose unknowns are physical potential or
normal velocity, rather than an indirect source-strength formulation. The initial
correctness prototype uses constant values on flat triangular numerical panels,
collocated at their centroids. This is an explicit prototype choice, not a promise
that every future geometry triangle must become an independent BEM unknown.

Use analytic constant-triangle single-layer integrals and oriented solid angles
for the double layer as the primary fast path. Use a separately implemented
transformed adaptive quadrature reference for self/near/far verification and a
bounded fallback when numerical cancellation defeats the fast path. A Duffy
transformation removes the self single-layer singularity. For a flat self panel,
`K_ii=0` as a principal value; the jump is added exactly once in `H`.
Never force diagonal row sums, insert a distance epsilon or return a failed
quadrature as a successful evaluation. Actual edge/vertex targets are rejected;
their trace coefficients are not always one half. Centroid self targets are
supported explicitly. Panels must form a nonintersecting boundary; global
intersection detection is not implemented by the panel cache or assembler.

Measure both kernels on the same physical length scale. A pair estimate must
satisfy `error <= atol + rtol*abs(reference)` for `S/panel_diameter` and `K`.
These are requested tolerances, not certified error bounds; independent checks
below establish their adequacy. Defaults are `rtol=1e-8`, `atol=1e-10` and at
most eight fallback order increases. An exhausted work budget is an explicit error.

Scale the mixed system with a declared positive length `L` and potential scale
`Phi`: `phi_hat=phi/Phi`, `q_hat=L*q/Phi`, `S_hat=S/L`. On a Neumann panel the
unknown column is `H[:,j]`; on a Dirichlet panel it is `-S_hat[:,j]`. Thus

```text
A x = S_hat[:,N] q_hat_known[N] - H[:,D] phi_hat_known[D].
```

For the first bounded pure-Neumann verification, check physical flux compatibility
and impose an area-weighted potential gauge through a bordered constrained solve.
Report the gauge multiplier and the residual of the original BIE independently;
discrete compatibility error must converge, not be hidden by a least-squares
projection or by deleting the worst row. Mixed cases with Dirichlet area and
decaying exterior problems do not receive an arbitrary constant-mode constraint.

Use pivoted dense [LU](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html)
and reuse a geometry-stage factorization for multiple right-hand sides. Estimate
reciprocal 1-norm conditioning from the factors, avoiding a fresh SVD on every
right-hand side. Report scaled backward error
`||b-A*x||_1 / (||A||_1*||x||_1 + ||b||_1)`, treating the exact zero system
separately. Reject non-finite data, singular factors and `rcond < 1e-12`.
For bounded pure Neumann these checks apply to the **bordered constrained
system**, not the singular original `H`. The original BIE residual is reported
separately and must converge with discretization; it is not required to reach
`1e-12`. Backward error of the solved system below `1e-12` is necessary, not
evidence of field accuracy.
Geometry, boundary tags or numerical operator changes invalidate a factorization;
matching panel counts alone never justify reuse.

## Interfaces, controls and cost

`topwater.numerics.LaplaceControls` owns the quadrature requests, fallback depth,
target block size and dense-operator byte budget. `DenseOperatorPlan(N, controls)`
computes and checks the storage for two square float64 operators before allocation.
Its version-1 JSON fixture is a **preflight record**, not an assembly result.

The first assembly implementation returns physical-unit `S` and jump-free `K`,
each a finite float64 `(N,N)` array, and a measured diagnostics record. That
record must include panel/target counts, controls, assembly wall time, operator
bytes, process peak resident memory, backend/thread metadata and fallback counts.
An unavailable analytic error estimate is `null`, not a fabricated zero. Actual
operator results and measured diagnostics serialization belong with that
implementation; no empty matrices or invented timings are introduced here.

Two dense operators need exactly `16*N*N` bytes: 0.25, 4 and 64 MiB for the
canonical 128, 512 and 2048-panel cases; 32768 panels require 16 GiB before any
factorization or scratch space. Block pairwise work over target rows. Never
allocate a global `(N,N,3)` displacement tensor. The initial block-row default
128 and operator budget 512 MiB are adjustable development controls, not tuned
optima or total-process memory guarantees. Benchmark block sizes and compiled
kernels at matched accuracy before changing the default.

Wall time is the primary performance objective; peak memory is second. Both are
subject to correctness. Compare implementations at the same measured error and
case, not only the same input panel count. A slower memory-saving change needs a
documented reason. Record cold setup separately from repeated work. For each
case use a warmup and at least five measured repeats; preserve individual samples,
minimum and median wall times. Measure peak RSS in a fresh process per case,
including setup/import baseline, and report operator bytes separately.

Each numerical branch benchmarks implemented panel integration, assembly,
factorization, repeated multi-RHS solves and field evaluation as applicable.
Record geometry/unknown counts, CPU/OS, Python/NumPy/SciPy, BLAS implementation,
actual thread count, block size and controls. Compare one BLAS thread with the
host default where it affects the measured operation. Do not invent a universal
runtime threshold before obtaining a baseline; explain regressions on comparable
hardware, including any accuracy benefit. Benchmark data go under `artifacts/`,
with concise accuracy/time/memory evidence in the PR.

### Initial kernel measurements (2026-10-21)

Windows 11, Intel family 6 model 170 (22 logical CPUs), Python 3.12.8,
NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.0; serial kernel, one BLAS thread,
128-row blocks. Medians of five warm runs, each case in its own process and
with its own empty compilation cache:

| Canonical panels | Assembly | 128 bulk targets, with gradients | Assembly peak RSS |
| ---: | ---: | ---: | ---: |
| 128 | 1.6 ms | 3.8 ms | 153 MiB |
| 512 | 29 ms | 9.0 ms | 157 MiB |
| 2048 | 468 ms | 44 ms | 217 MiB |

Cold first calls, including compilation, took 2.3–2.9 s. RSS includes imports,
compiler and scratch allocations; the 2048-panel operators alone occupy 64 MiB.
The bounded-domain constant-potential identity had maximum error `6.5e-15`;
all measured cases used zero quadrature fallbacks. Independent kernel/gradient
accuracy tests are separate from that identity.

The checked NumPy prototype took 2.82 s with 236 MiB peak RSS at 2048 panels and
128-row blocks. The fused kernel is about six times faster in this comparison;
the duplicate NumPy formula was removed. Smaller cases pay a larger relative
compiler-memory cost. Block sweeps from 16 to 128 rows did not establish a
consistent runtime winner, so the default remains 128. These are panel costs,
not boundary-solve, timestep or Capytaine-correlation results.

## Verification gates and implementation order

These are predeclared prototype targets, not claims that later stages have passed.
Changes require an explicit rationale and rerunning dependent cases. Tightening
quadrature must isolate integration error from panel-discretization error.

| Gate | Independent reference and acceptance |
| --- | --- |
| Panel kernel | Self/near/far triangle integrals against analytic special cases and separate transformed quadrature. `S/d` and `K` error at most `1e-10 + 1e-8*abs(reference)`. Include centroid approaches from both sides at `h/d=1, 1e-1, 1e-3, 1e-6`, a near-edge projection and near-zero derivatives. The dimensionless bulk gradients `grad_x S` and `d*grad_x K` have tolerance `1e-9 + 1e-7*abs(reference)`; boundary traces are tested separately. |
| Assembly signs | Canonical closed interior/exterior constant identities have maximum absolute error `1e-8`, without row-sum repairs. Normal reversal, proper rotations, translations and length scaling agree within `1e-10 + 1e-9*abs(scaled reference)`. |
| Boundary solve | Linear/quadratic harmonic fields on the actual canonical panels, an external point source for the bounded interior, and a strictly internal source for the decaying exterior. Use three successive refinements. Area-weighted RMS potential/flux errors and off-boundary velocity errors must decrease and reach 1% of the declared `Phi`, `Phi/L` scales. Algebraic, gauge and original BIE residuals are separate from those errors. |
| Calm water and geometry | Reproduce existing exact envelope and shell-loading references. Later clipping/pressure checks use the same canoe at three resolutions; displaced-volume, buoyancy and restoring errors reach 1% with decreasing refinement error. Calm-water stillness uses `max(abs(eta))/L` and `max(abs(u))/sqrt(g*L)` below `1e-10` in the dedicated zero-input test; no relative division by zero energy. |
| Waves | Linear dispersion/velocity and a published finite-amplitude reference, then crossing and broad directional cases. After five periods, regular-wave amplitude error and phase error in radians each reach `0.01`; report reflection, volume and energy drift separately, with spatial/time/domain refinement. Reject unsupported steepness/breaking, rather than fitting filters to pass. |
| Pressure and coupling | Still-water loading, prescribed radiation/diffraction and subsequently free response. Check load signs, phase, reference transformations, reciprocity where applicable, work/energy budgets and three spatial/temporal resolutions. Pressure uses the Eulerian potential derivative, including motion/remeshing corrections. No damping fitting precedes these gates. |
| Capytaine correlation | Match actual immersed geometry, CG/rotation centre, loading, water depth, gravity/density, frequency, heading, units and complex time convention. Normalize added mass as `D^-1 A D^-1/(rho*L^3)` and damping as `D^-1 B D^-1/(rho*omega*L^3)`, with `D=diag(1,1,1,L,L,L)`. Normalize excitation by `D^-1 F/(rho*g*amplitude*L^2)`. Each family is the complete nondimensional added-mass matrix, damping matrix or excitation vector over the predeclared frequency/heading/mode set. Use `abs(z-z_ref) <= 0.05*abs(z_ref) + max(1e-10, 0.005*max(abs(z_ref)))`; the floor covers an entirely zero family. Converge both solvers and report phase/couplings, not merely correlation coefficients. |
| Cost | Publish wall time first, peak memory second, and numerical error together. Initial operator cases use canonical 128/512/2048 numerical panels. A kernel-only agreement or small linear residual never qualifies a wave-body solver. |

Later free-surface stepping uses graph elevation and surface potential with
the [Bogacki-Shampine RK3(2) pair](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK23.html):
third-order accepted updates, an embedded second-order estimate and adaptive
rejection. The reason is oscillatory stability, not assumed speed: Heun has
`|R(i*theta)|^2 = 1 + theta^4/4`, whereas the accepted third-order polynomial has
`1 - theta^4/12 + theta^6/36`, with a finite stable imaginary-axis interval.
The extra stage must earn its cost at matched wave error. Record actual fluid
solve counts, rejected steps and wall time per simulated second. Reuse its final
stage only when it exactly matches the next starting state; remeshing, filtering
or projection invalidates that reuse. Fluid/body acceleration
must be solved consistently at each stage using acceleration-potential auxiliary
problems and a coupled six-motion solve. Filtering/absorption is explicit in
energy accounting. Actual integration controls and evolved-state APIs are added
with those stages, not as unused placeholders now.

Later saved results must contain time, actual pose/velocities, loads by origin,
diagnostics, timings and computed free-surface snapshots. A restart requires the
fluid, geometry and integrator state plus input provenance and wave realization.
Rim immersion and complete loss of contact, including invalid initial states,
terminate the run against the local evolved surface before accepting invalid
loads. The artificial cap is never a watertight deck.

Capytaine's [theory manual](https://capytaine.org/stable/theory_manual/theory.html)
uses an indirect source unknown and the opposite fundamental-solution sign.
Compare physical observables after convention conversion, not its source density
to Topwater's `q` or raw matrix entries without conversion. Use
`Re(Phi*exp(-i*omega*t))`, so linear dynamic pressure amplitude is `i*rho*omega*Phi`.
Capytaine is a benchmark dependency, never a substituted runtime load model.

The progression is panel kernels, assembled mixed boundary problems, waves,
moving wetted geometry/failure detection, pressure, then coupled motion. Each
stage must earn its own gates. Numerical verification and comparison with another
potential-flow solver do not constitute experimental validation of the canoe.

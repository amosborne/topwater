# Numerical contract and verification gates

This is the decision record for the first numerical implementation. It specifies
work to verify; it does not certify an operator, a wave solver or hydrodynamic
correlation that has not been implemented. The canonical canoe and loading in
`tests/reference_hull.py` remain the one shared hull reference.

## Production runtime and near/far decomposition

The user's target is 600 simulated seconds in at most 600 wall-clock seconds on
the development machine (Core Ultra 7 155H, approximately 32 GiB RAM, RTX 2000 Ada
Laptop GPU with 8 GiB VRAM). GPU use is allowed. Speed is the first priority,
memory the second, subject to the physical accuracy gates below. Count startup,
geometry updates, all attempted fluid/body stages, pressure recovery, diagnostics
and declared result output. A wave-only or fixed-coefficient load benchmark does
not qualify the coupled solver. Timestep and resolution must earn their accuracy;
they cannot be enlarged merely to meet the clock.

Real time is preferred; the user accepts 10--20 wall seconds per simulated second
for the complete validated case. The current implementation priority is full
nonlinear six-DOF motion and Capytaine validation, followed by speed optimization.
Body-free, static and kernel timings are component diagnostics only.

Dense global free-surface BEM is a verification choice, not a user requirement.
The next architecture candidate is a small detailed near-body potential-flow
domain coupled to efficient nonlinear far-field wave propagation. Do not build a
general backend framework before the coupling earns its accuracy and cost.
Both incoming waves and outgoing radiation/scattering must cross the artificial
interface. Nonlinear incident/disturbed cross terms remain part of the total
near-field flow.

The target marine-gyro hull geometry is basic. Keep the physical geometry mesh
distinct from any numerical subdivision: neither one unknown per geometry face
nor thousands of refined hull panels is a requirement. Start with the smallest
adequate discretization. Vary hull resolution, near-surface resolution, far-field
wave resolution and interface position independently, and justify increases by
load/motion error rather than mesh appearance. The curved canonical canoe remains
the common validation hull; its required geometry resolution need not transfer
to a hull made of a few planar faces.

With opposite fluid-outward interface normals, the coupled solution must satisfy
`phi_near = phi_far + a(t)` and `q_near + q_far = 0`. Potential continuity supplies
tangential velocity continuity; flux matching supplies normal velocity continuity.
Bernoulli pressure also requires consistent time derivatives and gauges:
`C_near = C_far + a'(t)`. A common `C=0` allows only a time-independent offset.
On a moving interface with velocity `w`, recover the Eulerian derivative using
`phi_t = d(phi(X_interface,t))/dt - w dot grad(phi)`. At the free-surface junction,
share elevation and gauge-aligned surface potential at the same location and time.
These are coupled matching conditions, not independent prescriptions of both
Dirichlet and Neumann data to each subdomain. Smooth visual blending alone does
not establish continuity of velocity, pressure or energy transfer.

The first coupling falsifier is two finite-amplitude oblique wave packets crossing
the interface in opposite directions in a body-free domain. Compare with an
unsplit reference at matched accuracy. Report interface potential, velocity,
pressure and elevation errors, artificial reflection, volume and energy budgets,
and complete wall time and peak memory. Refine space/time and move the interface:
its position must not materially change the physical result within the demonstrated
discretization error. Then repeat with the canonical canoe, including radiation,
diffraction, pressure and six-motion coupling and the Capytaine gates below.
Benchmark the full coupled path before claiming real-time capability. No near/far
implementation or real-time result is qualified yet.

The [spectral coupled BEM study](https://doi.org/10.1016/j.joes.2023.07.003)
supports investigating this split, but its sequential spectral/Green-function
formulation neglects the local residual potential outside its BEM region.
That approximation must be measured; a smooth transition is not a substitute
for satisfying the global boundary-value problem. Fixed-body diffraction in
that study does not qualify a moving canoe or this machine's runtime.

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
method, not a second analytic backend. `BoundarySystem` now owns its operators,
boundary assignments, scaling and reusable LU factors. Its configuration is
immutable, and it retains no redundant unfactored mixed matrix just for
residuals. Geometry changes require rebuilding the system. No geometry cache or
extensible backend framework is implemented. The prototype still fails the
boundary field-accuracy gate described below.

The first verification prototype uses a **direct Rankine-kernel BEM**, whose unknowns are physical potential or
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
Well-separated cancellation uses order-refined quadrature. Close scalar pairs
instead recompute the analytic integrals and geometry from the float vertices
using standard-library decimal arithmetic, from 32 through at most 96 digits.
Success requires agreement between successive precisions plus a binary64
rounding allowance within the same requested tolerance. `max_refinement` also
bounds these precision increases. Close gradient cancellation still fails
explicitly. Independent edge-aligned integration verifies the thin image
interactions found while moving the canonical canoe through its waterline;
agreement of global Gaussian orders is insufficient for those narrow peaks.

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

`BoundarySystem.solve` accepts physical full-panel potential and normal-velocity
arrays, with one or several RHS columns; entries where a field is unknown are
ignored. Its unknowns have common units (`phi` and `L*q`); a common potential
scale cancels from the relative diagnostics, so there is no redundant API knob.
Bounded all-Neumann data must satisfy discrete area-weighted flux compatibility;
the constrained solve imposes an area-weighted zero-mean potential. Exterior
data may have nonzero net flux balanced at infinity. Diagnostics expose actual
net flux, compatibility, gauge multiplier, original BIE residual and backward
error separately.

### Boundary prototype limitation

`uv run --locked python -m tests.verify_laplace` writes independent field errors
to `artifacts/laplace-verification.json` and deliberately exits 1 while gates are
unmet. For the bounded linear field `phi=x/3.5`, the dimensionless flux RMS errors
at 128/512/2048 panels are `0.290`, `0.190`, `0.125`, with declared `Phi=1` and
`L=3.5`. At 2048 panels the all-Neumann potential error is `0.00717`, but its bulk
velocity error is `0.0259`. Algebraic backward errors are small in both cases;
they do not excuse those field errors. Source fields also expose insufficient
flux recovery. The all-Neumann bounded point-source test now uses independently
integrated panel-average flux; surface velocity recovery is not implemented.

Constant panels are therefore a verified kernel and algebraic prototype, not a
qualified wave-solver discretization. The next method decision must improve
nonlinear harmonic and point-source fields, near-rim behavior and off-boundary
velocity at matched complete cost. Exact reproduction of a polynomial used to
construct a method is insufficient. Stage 1 remains open.

### Instantaneous flat-water fluid inertia

At startup from rest with a flat surface, atmospheric Bernoulli gives the
acceleration potential `psi=phi_t=0` on that surface. The same spatial problem
is the infinite-frequency linear limit. Eliminate the plane exactly with
`G_D(x,y)=G(x,y)-G(R(x),y)`, where `R` reflects the target across the water level.
This is the reflected Rankine limit in
[Capytaine's infinite-depth Green function](https://capytaine.org/v3.0.0/theory_manual/inf_depth_green_function.html).
It introduces no finite-domain boundary, wave-frequency table or reflected
geometry cache. `assemble_dirichlet_halfspace` reuses the existing kernel and
cached source panels with bounded target-block scratch.

The smooth image self term remains on `K`'s diagonal. It is not the trace jump;
`H=K+.5I` still adds that jump exactly once. `BoundarySystem` now computes the
actual `||H||_1` during its tiled factor-buffer fill, including this diagonal.

For physical wetted panels with fluid normals, define each row of
`Q=[n, (x-reference) cross n]`. The six simultaneous auxiliary solves give
`H*Psi=S*Q`. Pressure then gives `A=rho*Q.T*diag(area)*Psi` and fluid reaction
`-A*acceleration`. Rows and columns are world translations followed by rotations
about the stated reference. `instantaneous_added_mass` returns the raw matrix
without symmetrization, auxiliary potentials and conditioning/reciprocity/energy
diagnostics. Its diagnostic normalization length is the recorded maximum wetted
world-axis extent; it can change under yaw. No separate inertia load is added on
top of the same dynamic pressure in a future coupled calculation.

Independent checks include target versus source reflection, rationalized
image-kernel quadrature near the waterline, the smooth self diagonal and single
jump, source-plus-mirror bulk convergence, zero plane potential/tangential
gradient, yaw/translation covariance, density scaling and moment-reference
transformations. The bulk-gradient reflection includes its vertical chain-rule
sign, not just a difference of gradient arrays.

Capytaine 3.0.0's independent indirect radiation solve at `omega=inf` checks the
complete six-motion matrix on 140/460/1644 identical physical wetted triangles,
with the same CG, density and water depth. This uses the predeclared added-mass
normalization and 5% plus family-floor comparison above. Fine-mesh diagonal
differences are below 2%; reciprocity and positive energy are checked separately.
The infinite-frequency limit cancels frequency factors as described in
[Capytaine's force conventions](https://capytaine.org/v3.0.0/theory_manual/theory.html).
Both solvers' raw matrices are retained in CI artifacts. Exact coordinate
deduplication supplies shared reference-mesh indices without moving any triangle.

**Same-mesh correlation is not continuous-hull convergence.** The heave entries
are `293.91/285.76/281.47 kg` in Topwater and `297.62/288.84/283.11 kg` in Capytaine
over those three refinements. Some other modes still change by more than 1%.
This qualifies neither a mesh-independent inertia value nor finite-frequency
radiation, damping, excitation, evolving-surface coupling or the broader boundary
field gate. This flat-plane calculation is verification/initialization, not a
fixed-coefficient replacement for the nonlinear fluid solve.

`uv run --locked python -m tests.benchmark_inertia` times heave equilibrium,
clipping, numerical-panel construction, image assembly, one LU, six RHS and
matrix integration together. On the development Windows host with one BLAS
thread, warm medians at 140/460/1644 wetted panels were `0.0082/0.0529/0.7240 s`;
peak process RSS was `202/207/268 MiB`. At 1644 panels, assembly took `0.582 s`,
factor setup `0.117 s`, and six RHS `0.0147 s`; operators/factors occupy
`41.2/20.6 MiB`. First full calls with empty JIT caches took `7.45/7.54/8.15 s`.
Reports retain all five warm samples, first-call stages, actual library threads
and memory. These are initialization costs, not timestep timings.

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
host default where it affects the measured operation. The complete solver must
meet the development-machine real-time target above; component baselines explain
where its budget is spent. Explain regressions on comparable hardware, including
any accuracy benefit. Benchmark data go under `artifacts/`,
with concise accuracy/time/memory evidence alongside the experiment. Work is
organized in branches and commits; pull requests are not part of this workflow.

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

### Boundary setup and reuse measurements (2026-10-21)

Same host and versions as above; five warm repeats after isolated cold
compilation. Setup includes validation, mixed-matrix construction, LU, condition
estimation and norm diagnostics. The six RHS are independent harmonic fields,
not six hydrodynamic radiation results.

| Panels | Complete setup | One RHS | Six RHS together | Process peak RSS |
| ---: | ---: | ---: | ---: | ---: |
| 128 | 0.28 ms | 0.15 ms | 0.50 ms | 158 MiB |
| 512 | 5.0 ms | 0.54 ms | 1.8 ms | 172 MiB |
| 2048 | 235 ms | 6.9 ms | 21 ms | 268 MiB |

At 2048 panels, operators occupy 64 MiB and factors about 32 MiB. Cold assembly
took 2.47 s and first system setup 0.51 s including compilation. Profiling found
column norm scans and mixed-matrix construction consumed substantial setup time;
LAPACK norms on transpose views and a tiled fill into the one factor buffer
reduced setup from 430 ms without adding another dense matrix. Repeated library
discovery also dominated small cases, so already-loaded NumPy/SciPy thread
controls are discovered once.

The host's 22-thread BLAS default took 2.34 s for setup. The default now scopes
LU/condition estimation to one BLAS thread while repeated solves retain ambient
threading. `factorization_threads=None` preserves caller configuration. This
context temporarily changes process-wide limits; concurrent applications should
coordinate threading externally. No persistent thread limit is set at import.
Benchmark with `python -m tests.benchmark_boundary`; use `--threads 0` to retain
the host default for solves and `--factorization-threads 0` to test ambient LU.

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

Boundary verification samples full triangle interiors using independent
Gauss-Duffy quadrature, in addition to reporting collocation errors. A prescribed
field can match at every collocation node while its panel representation varies
incorrectly between them. The report retains the original `Phi=1`, `L=3.5 m`
scales and adds relative surface-potential, normal-velocity and bulk-velocity
RMS errors; changing scales is not a substitute for reducing physical error.
The qualification gate includes both sampled and full-surface errors. A bounded
pure-Neumann reference uses the continuous area-mean potential gauge.

The bounded point-source Neumann case now supplies independently integrated
panel-average flux. Two quadrature orders agree and the integral has zero net
flux without subtraction/projection. The Neumann compatibility tolerance remains
`1e-12`. At 128/512/2048 faces its potential error decreases
`0.01792 → 0.006868 → 0.002198`; its original BIE defect decreases
`0.002757 → 0.0008668 → 0.0002638`. The bordered algebraic solve reaches roundoff;
the original equation defect remains a separate discretization error.

Contained higher-order and auxiliary-source probes have not earned a runtime
API. Independent full-surface checks show that interpolation of the prescribed
fields can dominate their error, while auxiliary-source placement introduces
conditioning and between-node accuracy problems. They do not qualify the current
constant-panel method or justify building a general discretization framework.
Surface-velocity recovery and the full boundary accuracy gate remain unfinished.

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

The first independent comparison is kernel-only: run
`uv run --locked --group correlation python -m tests.test_capytaine_rankine`.
The validation group pins Capytaine 3.0.0; it is not a runtime dependency. A
dedicated CI job requires its native import and retains JSON comparison reports.
Exactly the same 128/512/2048 canonical triangles are used, with source-normal
differentiation and no diagonal trace contribution.

The [tagged Rankine implementation](https://github.com/capytaine/capytaine/blob/v3.0.0/libs/libDelhommeau/src/Green_Rankine.f90)
uses a point approximation beyond seven panel radii, an angle cutoff below
`1e-4*radius`, an absolute edge cutoff `1e-5`, and a tiny-edge threshold. The
comparison classifies those branches from geometry, not from observed agreement.
Its [PI initializer](https://github.com/capytaine/capytaine/blob/v3.0.0/libs/libDelhommeau/src/constants.f90)
also rounds through default real precision: the installed constant is
`3.1415927410125732`. The report retains raw differences and explicitly converts
the reference normalization by `-pi_reference/pi` before strict kernel checks.

At 2048 panels, 651954 of 4194304 pairs fall in the checked analytic region;
maximum converted S/K differences are `4.4e-17`/`6.5e-15`. Every source panel has
nonself coverage. The largest disagreements outside that region are checked
separately against independent quadrature, without weakening Topwater's gate.
These kernel checks establish neither boundary-solve accuracy nor hydrodynamic
force correlation. Instantaneous added mass is compared separately above;
finite-frequency added mass, damping and wave excitation remain future work.

Independent calm-water hydrostatics now clips the original physical hull, never
its artificial cap. Exact triangle pressure and moment integration is checked
against tetrahedral displacement/first moments and finite differences of all six
world-axis load derivatives. A volume-only root solves fixed-attitude heave
balance; it does not solve attitude equilibrium or certify stability. Static
rim contact and loss of water contact terminate evaluation. Evolving-surface
failure detection and still-water time integration remain unimplemented.

The analytic ellipsoid checks volume, buoyancy center, waterplane area and
roll/pitch restoring at 512/2048/8192 total hull faces. Capytaine independently
clips 128/512/2048 original hull faces at equilibrium, at two yaw headings.
Its `Gauss-Legendre 2` rule integrates the quadratic planar-face moments instead
of its default centroid approximation. Displacement, wetted/waterplane areas and
buoyancy center agree within `1e-12 + 1e-10*abs(reference)`; restoring entries
within `1e-8 + 1e-10*abs(reference)` in physical SI units at the same CG, mass,
density and gravity. These are hydrostatic comparisons, not wave-force gates.

The progression is panel kernels, assembled mixed boundary problems, waves,
moving wetted geometry/failure detection, pressure, then coupled motion. Each
stage must earn its own gates. Numerical verification and comparison with another
potential-flow solver do not constitute experimental validation of the canoe.

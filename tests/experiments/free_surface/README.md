# Symmetric free-surface BEM experiment

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
energy buffer. Dense assembly and factorization still dominate moving stages.

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

# Symmetric free-surface BEM experiment

This checkout-only experiment tests a concrete path from nodal free-surface
state through moving geometry, a BEM solve and nonlinear time derivatives.
It is not an installed Topwater solver, a completed Stage 1/2, or an open-water
wave-body simulation. The rectangular tank is a fluid verification domain;
`hull_check.py` uses the existing canonical canoe, not another test hull.

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
does not change. The unoptimized dense temporary ownership remains experimental.

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

## Reproduce

Run from a locked checkout; outputs go to ignored `artifacts/`:

```powershell
uv run --locked python -m tests.experiments.free_surface.verify
uv run --locked python -m tests.experiments.free_surface.linear
uv run --locked python -m tests.experiments.free_surface.evolution
uv run --locked python -m tests.experiments.free_surface.integrate
uv run --locked python -m tests.experiments.free_surface.integrate --rtol 1e-5
uv run --locked python -m tests.experiments.free_surface.shape_check
uv run --locked python -m tests.experiments.free_surface.hull_check
uv run --locked python -m tests.experiments.free_surface.benchmark
```

`verify` reports `qualified=false` and remaining work; command completion is not
a qualification gate. `integrate` rebuilds the BEM at every RK23 stage and takes
about a minute for the finer short case on the development host. The benchmark
uses fresh workers, empty compilation caches and five warm samples. Its scope
is geometry through the flat solve and independent diagnostics, **not** a moving
timestep, whose actual cost is recorded by `integrate`.

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

## Evidence and limits

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
For `kd=1`, the fine continuous quadratic errors are 2.41% kinematic and 1.57%
dynamic; the nodal result is not a one-percent whole-field gate pass.
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

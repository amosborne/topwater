# Nonlinear moving-boat baseline and verification

The experimental baseline has one execution path: adaptive RK23, exact integer
clipping with guarded chart replay, analytic chart and momentum derivatives,
and the sample-local reverse energy gradient. Every RHS uses a fresh current
fluid operator. It evolves all six rigid-body degrees of freedom and the
nonlinear free-surface state.

Physical inputs use `topwater.SimulationInputs`: a convex hull mesh, mass, CG,
inertia tensor, initial pose/velocities and explicit directional wave components.
The reference boat is a fixture, not a required mesh topology. Coplanar facets
and pressure-bearing groups are derived from input geometry; artificial opening
closures never receive pressure. Another supported boat needs a new input file,
not solver edits. The intended scale remains 18-30 ft boats in ordinary ocean
conditions; incident wavelengths must be at least twice beam and frequencies
at most 0.6 Hz. Breaking, slamming, flooding and loss of contact are unsupported.
The breaking/slamming detection is not yet qualified.

## Run the current baseline

Export the example physical data, then run the same driver used for other inputs:

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
uv run --locked python -m tests.reference_boat --output artifacts/reference-boat.json
uv run --locked --with 'cupy-cuda12x[ctk]==14.2.0' python -m tests.experiments.free_surface.boat_moving --inputs artifacts/reference-boat.json --duration 1 --max-step .2 --max-rhs 180 --budget-seconds 240 --output artifacts/boat-baseline.json
```

Refinement, quadrature order and ODE tolerances are convergence controls.
There are no runtime switches for old formulations, derivative implementations,
operator predictors or integrators. Failed physical or geometry guards stop
the case. Clipping treats the supplied coordinates exactly and rejects
nonconvex geometry. Dry artificial closures currently require a conservative
global height certificate; failure to certify is a geometry rejection, not
proof of flooding. This can reject locally valid cases.

Reports own serialized physical inputs, actual producer hashes and
accepted RK dense polynomials; they are not restart files.

## Verification and limits

Independent forward derivatives, scalar finite differences and h/h2 momentum
controls remain in their CUDA and action audit harnesses. Their implementations
do not form alternate moving-driver branches. Shared panel geometry and pair
integrals each have one numerical owner. The old pressure fallback, operator
reuse and rejected material-mesh moving paths remain in Git history.

`boat_motion_compare` compares current-format reports using their owned accepted
polynomials and physical-input identities. Historical schemas and reconstruction
fallbacks are not part of the baseline comparison path. Earlier measured results
retain their original source archives and the commits that produced them.

The current case still uses a 12 x 6 x 2 m diagnostic tank, not an open-ocean
boundary. Pressure/action and spatial discrepancies remain unresolved. No run
has established 0.1 degree per rotation and 1 cm per translation throughout
60 seconds, or the required <=20 wall seconds per simulated second.
`hos_graph.py` is a separately verified nonlinear spectral component; it is
not yet a coupled spatial near/far solver.

[Issue #30](https://github.com/amosborne/topwater/issues/30) owns active targets,
experiments and adoption decisions. Detailed evidence is in `docs/research`;
previous versions of this catalog describe the archived investigations.

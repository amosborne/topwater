# Independent hydrostatics quadrature review, 2026-10-24

Bounded coordinator assignment: independently derive exact triangular moments,
review producer's hydrostatics-only correction, and identify decisive regression
expectations. No correction implementation or numerical launch by this reviewer.

## Source observation

Installed Capytaine 3.0.0 uses default face-centroid integration in
`meshes/meshes.py:472-486`. `meshes/surface_integrals.py:79-87` evaluates center of
buoyancy through integrals of `n_i*x_i^2/(2V)`. The rigid-body formulas in
`bodies/hydrostatics.py:190-218` require quadratic waterplane moments. A centroid
rule is exact for constants and linear functions, not those quadratic terms.
The same physical surface represented with an asymmetric triangle partition can
therefore produce asymmetric numerical moments without physical asymmetry.

Native `Mesh.faces` pads triangles `[a,b,c]` to `[a,b,c,c]`.
`meshes/geometry.py:88-121` supplies four tensor-product two-point Gauss-Legendre
samples with the mapped Jacobian. On that collapsed square, write

`r = c + u*(v*(a-c)+(1-v)*(b-c))`, `u,v in [0,1]`.

The Jacobian is `u*norm(cross(a-c,b-c))`. A physical quadratic in coordinates
times this Jacobian has degree at most three in `u` and two in `v`; two-point
Gauss-Legendre is exact through degree three in each variable. Therefore native
`with_quadrature('Gauss-Legendre 2')` exactly integrates the required planar
triangle moments, apart from floating-point roundoff. This conclusion derives
from the source mapping, not its method name alone.

## Independent analytic reference

For triangle vertices `a,b,c`, area `A`, `s=a+b+c`, and
`Q=a*a' + b*b' + c*c'`, barycentric integration gives

- `integral 1 dS = A`;
- `integral r dS = A*s/3`;
- `integral r*r' dS = A*(s*s'+Q)/12`.

These follow from `integral lambda_i^2 = A/6` and
`integral lambda_i*lambda_j = A/12` for distinct indices. They are an independent
closed-form expectation for all ten Cartesian monomials through degree two.

For the unit right triangle `(0,0,0),(1,0,0),(0,1,0)`, the expected integrals are
area `1/2`, `x=y=1/6`, `x^2=y^2=1/12`, and `xy=1/24`. Centroid integration instead
gives `x^2=1/18`, so this single case decisively rejects the old integration.
Also test a translated/tilted scalene triangle and cyclic/reversed vertex
permutations to avoid accidentally testing only a coordinate-plane special case.

For the actual outward-wound wet surface, the independent moments yield

`V_i = sum_faces n_i * integral x_i dS`,

`CB_i = sum_faces n_i * integral x_i^2 dS / (2V)`.

The missing waterplane closure contributes zero to these expressions at `z=0`.
Waterplane moments of `p(x,y)` follow from
`integral_waterplane p dA = -sum_faces n_z * integral_face p dS`.
Use the full analytic moments, not expected symmetry alone, to verify the native
quadrature result.

## Physical invariants and limits

For the canonical lateral mirror symmetry and `CG_y=0`, exact integration gives
`CB_y=0`, heave/roll coupling zero, roll/pitch coupling zero, and pitch/yaw
coupling zero. Roll/yaw coupling zero additionally depends on the equilibrium
condition `CB_x=CG_x`; it does not follow from fore/aft symmetry, which this hull
does not possess. A nonzero heave/pitch coupling is permitted because the
waterplane centroid need not lie at the CG. Dry products of inertia also must
not be zeroed merely to enforce a convenient diagonal matrix.

With rotations about CG, the dimensional nonzero hydrostatic terms are

`C33 = rho*g*Aw`,

`C34 = rho*g*integral_waterplane(y-yG)`,

`C35 = -rho*g*integral_waterplane(x-xG)`,

`C44 = rho*g*(integral_waterplane(y-yG)^2 + V*(zB-zG))`,

`C45 = -rho*g*integral_waterplane((x-xG)*(y-yG))`,

`C55 = rho*g*(integral_waterplane(x-xG)^2 + V*(zB-zG))`.

The selected source also includes equilibrium-residual terms
`C46 = -rho*g*V*(xB-xG)` and `C56 = -rho*g*V*(yB-yG)`; check these against
independent moments rather than manually setting them to zero. At valid
equilibrium the residual terms vanish to the independently established pose
tolerance and the expected restoring matrix is symmetric.

The correction must operate only on hydrostatic integration/export. Preserve the
17-panel frequency-dependent hydrodynamic coefficients and exact producer
identity; do not change BEM quadrature or silently claim a new hydrodynamic mesh.
Source review and executed artifact equality/regression evidence are distinct.

Read-only implementation review of `reference/capytaine/correct_hydrostatics.py`
passes intended isolation: a separate native GL2 hydrostatic body is built on the
same 17 triangles; original frequency-domain variables/coordinates are compared
exactly after export; source manifests and dependency hashes are checked; no BEM
solve is launched by the correction.

Read-only review of `reference/capytaine/check_hydrostatics.py` passes independent
expectation coverage: unit and tilted/translated triangles, all six vertex
permutations, old centroid-rule rejection, actual 17-face analytic V/CB/full K,
physical lateral symmetry and equilibrium. The expected moments do not call the
native quadrature implementation. The tight comparisons are polynomial-roundoff
checks, not a claim of mesh-converged motion accuracy.

## Final independent executed signoff

PASS for the hydrostatics correction within this minimum 17-panel reference.
Independently inspected the authoritative evidence under
`artifacts/reference/minimum-20261024/native-gl2`:

- `hydrostatics-checks.json`: all eight checks pass; all twelve unit/skew triangle
  permutation cases agree at roundoff; actual hull full stiffness differs from
  independent analytic moments by at most `4.3655745685100555e-11`.
- `manifest.json`: complete, no new BEM solves, every retained native
  frequency-domain variable/coordinate compared exactly after export. Corrected
  CB is approximately `[4.05e-16,1.80e-17,-0.1300704594] m`. Physical heave/pitch
  coupling `-10123.376648462312` remains.
- Recomputed hashes of all six declared output files match that manifest.
  Correction source snapshot and current source match frozen hash
  `41b85ca58ba7ad63d2276399c68c751d66d1c373ff4ff4dd9cbb298e3a898ec0`;
  check source matches proposal hash
  `53c88c500a98e9a80b0c123197c447f749d35f0541b72742c02b5fbec734b44c`.
- `../hydrostatics-run.log` and the execution ledger record correction and check
  terminal exit zero. The lease was absent at this review. This reviewer did
  not rerun numerical work.

The correction author and this reviewer are distinct. Physical reliance may use
these corrected planar hydrostatic moments and preserved coefficient identity.
This signoff does not establish convergence or accuracy of finite-frequency
hydrodynamic coefficients, body motion, packet startup or video.

# Closed fluid boundary and pressure coupling: bounded checkpoint

Scope: complete the minimum watertight fluid-boundary construction and the
acceleration-potential/body coupling on the selected P0 Rankine boundary solver.
Use the canonical hull and preserved mild-case surface maps. The coupled
nonlinear time-domain baseline and reference qualification remain the outcome;
this checkpoint supplies necessary missing components, not a replacement goal.

The fluid boundary has physical wet hull, free surface, outer wall and bottom
tags. Share every contact vertex/edge, preserve the exact physical hull facets
and include no artificial hull closure in the pressure-bearing surface. Triangulate
each wet facet using its actual contact contour, retaining concave portions and
all contact vertices. Require a connected closed consistently oriented manifold,
positive volume, no unused/duplicate vertices and no collapsed panels. Whole
surface/hull overlap checks and input identities remain required.

Use fluid-outward normals. On the hull they are the negative of body normals.
Let N=[n,(x-CG) cross n]. For Eulerian chi=phi_t, the body-normal condition is
N*a + known terms. Given fluid velocity u and Hessian H at body panels, the known
term is n dot (Omega cross (Omega cross r)) - n dot H*v_body
+ (v_body-u) dot (Omega cross n). This includes the material normal derivative;
linear/angular acceleration about the world CG is solved, not prescribed.
The free-surface dynamic condition is chi=-g*z-|u|^2/2, and body pressure is
p=-rho*(chi+|u|^2/2+g*z). The pressure owner must not silently substitute zero
velocity/Hessian for missing nonlinear data.

For one immutable mixed-boundary factorization solve seven chi right-hand sides:
one known-data column and six body-acceleration columns. Retain raw
A=rho*N^T*area*chi_aux, with no symmetry fitting. Solve (M_dry+A)*a=known fluid
load+weight-gyroscopic load. Integrate hydrostatic force AND moment exactly on
the actual planar wet triangles; gravity appears once. Acceleration pressure
and kinetic terms use the declared P0 trace representation. All six couplings
must be returned and independently reconstructed from the solved pressure.

Acceptance before any numerical job: geometry checks above; manufactured
constant-pressure/rigid-rotation/known-harmonic derivative tests at 1e-10 SI-scaled
tolerance; direct pressure-load versus augmented-system balance at 1e-10 scaled;
and explicit failures for missing/nonfinite derivative data or incompatible tags.
Residual/energy/symmetry alone is not physical acceptance. Preserve full mass,
axes, units and source/input hashes; independently review signs and omitted terms.

The 24 m cylinder is a candidate computational boundary, not a seabed. No
zero-flux-bottom/infinite-depth equivalence is claimed (at0.1 Hz exp(-kH)=0.3807).
Outer/bottom chi data must be explicitly supplied. A homogeneous auxiliary
condition is a declared finite-domain approximation whose domain error remains
unqualified. Do not inject future Capytaine body pressure, force or trajectory.
Only physical initial fields and prescribed incident-wave data may enter runtime.

Bounded canonical probe: use the exact calm equilibrium, flat free surface and
empty wave arrays for the pressure solve; then u=H=0 is exact, not a substitution
for missing mild-state derivatives. Also construct the preserved full-amplitude
mild map's closed boundary as geometry-only evidence without a pressure call.
Archive actual producers, both physical packets and field-map bundle
before launch. Run serially with one native thread and a 180 s wall limit checked
between stages. Count a tiny manufactured tetrahedron numerical warm-up separately
as preparation, once, before measuring the calm pressure component; do not multiply
cold compilation/cache loading by time steps. Stop after the calm factorization/pressure solve and mild closure
or a concrete rejection; no one-second integration. Report panel count and all
costs separately. Pressure-only cost is a lower bound, never complete-RHS timing.
For the provisional explicit-midpoint dt=.01 s route, 200 pressure evaluations
per simulated second are needed: if this component alone exceeds200 s, reject
that route. Falling below does not establish practicality; gradients/Hessians,
state evolution, setup and validation remain necessary.

Exclusions: invented zero nonlinear terms in the mild case, a qualified RHS or
restart without derivative/compatibility evidence, reference-gate promotion,
long simulation, alternative solver formulation and major optimization.
Deliverables: implementation and focused physical/manifold tests, immutable
canonical probe or rejected candidate, independent review, durable evidence,
issue #30 update and fresh-context handoff. Stop at this checkpoint, retaining
all unresolved source/field/domain/motion gates and the active broader goal.

# Closed-interior mixed trace regularization for a near/far Taylor experiment

Status: derivation and source review only. No implementation or numerical jobs.
Scope is the current CLOSED diagnostic fluid boundary with fluid-outward normals,
normal_sign=+1, a conforming shared-node P1 trace, and interior target quadrature
samples. It is not automatically valid for the mirrored-hull exterior/infinite
domain added-mass operators with normal_sign=-1 or a decay-at-infinity gauge.

## The identity, with the current sign

Let G(x,y)=1/(4*pi*|x-y|), k(x,y)=dG/dn_y, and K the principal-value source
double-layer operator. With the current interior convention, at a target point
strictly inside a planar boundary face,

```
integral_closed_boundary k(x,y) dy = -1/2
H phi(x) = phi(x)/2 + integral k(x,y)*phi(y) dy
         = integral k(x,y)*[phi(y)-phi(x)] dy .
```

The negative sign follows directly from div_y[(x-y)/|x-y|^3]=-4*pi*delta(x-y).
Equivalently the exact oriented source solid angles sum to -2*pi on the
boundary. The diagonal planar source face has zero smooth kernel because its
height is identically zero; the other source faces supply the half jump through
their closed sum. The regularized formulation therefore has NO explicit self
jump. Adding the old area_i/6 term again would double-count it.

This is an analytic Green-identity subtraction, not a post-assembly row or gauge
projection. Constant subtraction alone does not regularize every single-layer,
hypersingular, contact or geometry derivative in the complete formulation.
In particular it is not a claim that near-contact interactions become smooth
enough to approximate cheaply.

## Exact compatibility with the existing finite target rule

At sample x_iq, use the actual target P1 value

```
phi_i(x_iq) = sum_a b_qa*phi_ia
```

with the same reference barycentric sample as the current target quadrature;
do not use the target centroid value. Let h_ijc be the three existing SMOOTH
source P1 coefficients, excluding the self jump, and define

```
t_ija = -area_i*sum_q w_q*b_qa*sign_j*angle_ijq/(4*pi).
Hreg_ij phi = sum_c h_ijc*phi_jc + sum_a t_ija*phi_ia .
```

At each target sample, exact analytic inner source integrals obey the closed
solid-angle identity. Summing the target terms over all sources gives the
current half-jump mass contribution. The current production target orders
3/6/12 integrate linear target barycentrics exactly, so this equals the current
area_i/6 jump for every target corner, up to floating arithmetic. A target rule
not exact for these linears need not match the historical analytic jump; do not
extend this finite-rule equivalence to an arbitrary low-order rule by assumption.
Self smooth source and target terms are both zero.

Required geometric assumptions are stronger than a generic mesh flag: actual
fluid-interior orientation, no missing/overlapping boundary patches, common
physical trace at shared nodes, valid target interiors and no inverted faces.
Corners/edges would have a different local solid angle, but the current Gauss
samples lie inside faces. Open-water truncations need their actual enclosing
boundary or the correct infinity contribution before applying this identity.

## What derivative equivalence means

On an admissible closed shared-node geometry, full Hreg(X)=Hraw(X). Consequently
their GLOBAL coordinate derivatives also agree if the mesh remains a closed
oriented indexed manifold under those coordinate variations. At a minimum the
current admissible chart contractions dX=JX*dq must agree. A nodal P1 value held
fixed has fixed sample value b_q*phi_i when its target point moves with the
reference barycentric sample; the target area and source/target kernel geometry
still vary.

Individual raw-versus-regularized PAIR coefficients/18 virtual derivatives do
NOT agree: they allocate the jump differently, and independent source-face
perturbations generally destroy the closed-surface identity. Do not use the
previous raw-pair equality diagnostic as an expected-zero comparison between
these representations. Use it only to validate each representation's own
forward/reverse derivative, followed by assembled global/action comparisons.
If a mesh is merely geometrically closed through constrained coincident nodes,
restrict equivalence to closure-preserving variations unless indexed closure
has been independently established.

## Critical limit: full Taylor is not improved by this rewrite

On a shared-vertex closed manifold, regularized and raw GLOBAL H0 and DH0 are
the same. Thus caching a first-order Taylor polynomial of EVERY pair yields
the same assembled H0+DH0*deltaX, to floating arithmetic. Regularization alone
cannot cure the earlier fully linearized 10.53/141.46 Nm held-out body-pose
errors. With merely constrained closure, both first-order chart derivatives
still agree; any discrepancy along a curved finite chart displacement is a
second-order extension effect, not a demonstrated accuracy improvement.

The old full raw Taylor audit already had H*1 about 7.2e-16/8.0e-16 on its
10/20 ms held-out states. Gauge leakage was not the demonstrated cause of
those large force errors. The relevant report is
`artifacts/worktrees/kernel-linearization/artifacts/kernel-anchor-audit-generic-q3-body1.json`.

The actual opportunity is a HYBRID: keep difficult near interactions exact at
the CURRENT geometry and Taylor-expand well-separated far interactions.
Ordinary raw near/far mixing can disturb the global solid-angle cancellation;
regularized PAIRS each annihilate constants algebraically, so this mixing
preserves that property without a global row correction. It still approximates
the operator on nonconstant traces, and needs held-out load tests.

## Five independent mixed coefficients, not necessarily six

The six coefficients satisfy sum(h_source)+sum(t_target)=0 analytically for
EACH pair. A direct fixed-pivot representation uses p=phi_i0:

```
Hreg_ij phi = sum_c h_ijc*(phi_jc-p)
             + t_ij1*(phi_i1-p) + t_ij2*(phi_i2-p).
```

This follows by changing basis INSIDE the regularized integral, not by replacing
a computed row after the fact. It needs three source and two target coefficients.
Constant nodal inputs give zero differences directly. Materialize H, H^T and
their derivatives from this same five-coefficient basis. Floating matrix
assembly/BLAS may still leave a small H*1 residual; do not claim bitwise gauge
annihilation from a mathematical basis identity. Direct difference-form H*phi
provides the stronger numerical constant-input check.

With V, the straightforward cache has seven outputs (V+3+3), or six using this
pivot basis (V+3+2). Both include all eighteen local geometry slopes and a base
value. Full directed-pair retained storage is nfaces^2*noutputs*19*8 bytes:

| Faces | Seven outputs | Six outputs |
|---:|---:|---:|
| 472 | 237,042,176 B (226.06 MiB) | 203,179,008 B (193.77 MiB) |
| 1088 | 1,259,503,616 B (1201.16 MiB) | 1,079,574,528 B (1029.56 MiB) |

These exclude pair indices, exact-near data, factors, E/H outputs, geometry
Jacobians and temporary workspace. The previous 256 MiB retained-cache limit
was an AGENT-CHOSEN prototype guard, NOT a user requirement. Speed is the user's
first priority and the device has 8 GB VRAM. A roughly 1.1--1.3 GB single-anchor
payload is a legitimate candidate to price, not a reason to reject the method
before measuring the full footprint. The hypothetical 256 MiB cap would fit at
most roughly 21%/25% far pairs at 1088 faces, limiting even free-far pair work
to about 1.27/1.33 speedup if per-pair costs are comparable; that conditional
ceiling is an argument AGAINST treating that arbitrary cap as a design limit.

No new memory measurements were taken during this derivation. Existing nearby
sizes give only partial accounting: a 1088-square float64 V array is 9.47 MB,
H at 1088 by 546 is 4.75 MB, and the 546-square E is 2.39 MB. Raw/symmetric V,
factor, outputs and workspaces may coexist; the CuPy pool may retain old
allocations, and derivative kernels can add substantial local/work storage.
The falsifier must measure peak allocated AND reserved VRAM plus process RSS,
including anchor construction, before deciding whether the highest-speed
single-anchor candidate fits. Do not report retained tensor payload as a
measured full footprint. Compression or a small-case fallback is optional;
neither is required merely to satisfy an invented 256 MiB constraint.

## Same approximate energy and its force

For the first falsifier hold the semantic chart, anchor and near/far mask fixed.
Use the current exact near pairs and the stored far coefficient Taylor
polynomials to assemble Vraw_hat and Hreg_hat. Symmetrize V exactly as the owner
does, keep current curls C(X), and recompute z_hat=Vhat^-1*Hhat*phi CURRENTLY:

```
Ehat = C.T*Vhat*C + Hhat.T*Vhat^-1*Hhat
dK = rho*[(dC phi).T*Vhat*(C phi)
          + .5*(C phi).T*dVhat*(C phi)
          + z_hat.T*dHhat*phi - .5*z_hat.T*dVhat*z_hat].
```

Exact-near derivatives are current; far derivatives are the SAME cached
coefficient slopes. No stale anchor phi/z/gradient or exact-current half jump
is spliced into the approximation. Current physical geometry, full incident
forcing, moving trace maps and all six acceleration columns remain current.
Cholesky failure or ownership changes reject/refresh; no jitter or row repair.

For bdot/momentum/pressure, the directional Ehat action must differentiate this
SAME hybrid energy. The existing exact-current directional action cannot be
borrowed by pretending this approximate assembly is its original uniform raw
assembly. A cheap body-pose gradient screen alone does not qualify the full
known hydrodynamic load or acceleration.

If blend weights later depend on X, include dw*(exact-Taylor) in the derivative.
A changing Boolean mask makes a piecewise model; simply switching masks while
ignoring that model change is not the derivative of one smooth energy. A fixed
first mask avoids this ambiguity. Geometric separation/bounding balls can
identify candidate far pairs, but their threshold is not by itself a load-error
bound, especially with long/slender triangles and waterline-adjacent panels.

## Refresh and FSAL ownership

An anchor/mask change changes the approximate RHS even at an unchanged physical
state. Give the approximation one immutable version during a whole RK trial.
If an audit/topology guard requires refresh, discard the trial, refresh from
the last accepted state and recompute its initial RHS. In particular a cached
FSAL endpoint derivative from the old version cannot seed the next step under
the new version. Accepted old dense polynomials describe the old integrated
step and cannot be retroactively relabeled. Solver restart or explicit current
f recomputation is required; do not silently change a fun closure mid-step.
Operator-error control remains separate from the integrator's temporal error.

## One bounded implementation falsifier, only after math review

Reuse the actual 472-face generic anchor and held-out 10/20 ms states. Fix one
conservative geometry-based near mask at the anchor: all touching pairs and
contact/slender near neighborhoods exact, with a true positive far-separation
margin over the held-out displacement budget. Record exact-near/far fractions
and predicted arithmetic cost before allocation, without rejecting a fast
large-cache candidate solely for exceeding 256 MiB. No threshold sweep or new
backend framework.

First require delta=0 regularized global H/E, admissible shape contractions,
all fifteen parent rates and six body-pose loads to match the fresh raw owner.
As a negative control, ALL-Taylor regularized must reproduce the previous raw
ALL-Taylor behavior; apparent large improvement there signals an ownership or
derivative mismatch to investigate. Then test the one fixed hybrid mask on
both held-out sizes, with current phi/z and derivative of the same energy.

Record whole fresh-reverse versus hybrid costs, including anchor build,
retained bytes, exact-near work and amortization at the observed roughly 2.5
RHS uses/chart. The earlier 0.7 s slope build can erase savings relative to the
new 0.09 s whole-gradient component. Also record gradient/action errors: a tiny
potential or field residual was previously compatible with 141 Nm load error.
Stop if the mask leaves too much exact work for a meaningful gain, amortized
cost loses, delta=0 controls fail, or both held-out load errors do not improve
substantially. Passing this cheap screen only authorizes full consistent
directional-action/known-load/pressure/acceleration controls before a trajectory.

The later root controls adopted the exact integer/reverse components at
`73bbdad`: 64 RHS evaluations cost roughly 69--89 wall seconds per simulated
second, with about 8.75 s gradient, 27 s cold geometry and 33 s other work in
the reported component breakdown. Deleting the gradient alone would yield
only about 1.11--1.15x overall speedup across that wall-time range. A cache must
therefore earn its extra construction/maintenance cost against the NEW reverse
baseline; the old 10/141 Nm rejected Taylor result is not an attractive speed
baseline anymore. At roughly 2.5 RHS uses/chart, a bounded mask count and an
honest warm-build amortization estimate should precede any cache implementation.
The component breakdown does not price all possible E/directional-action
savings, but it prevents interpreting a gradient-only gain as the remaining
order-of-magnitude whole-solver improvement. Issue 30 owns subsequent status;
this note records a reviewed derivation, not an implementation/adoption decision.

## Primary context

Analytic subtraction to regularize boundary integrals is established; the
authors' [non-singular fluid-mechanics formulation](https://arxiv.org/abs/1902.03342)
includes potential-flow Laplace problems, and their
[potential-problem formulation](https://arxiv.org/abs/1902.04168) develops broader
non-singular Laplace equations. The specific discrete coefficients, sign,
global-versus-pair derivative limitation and Taylor equivalence above are
derived here from the actual `cuda_pair.cuh` owner, not attributed as a method
or performance result from those papers.

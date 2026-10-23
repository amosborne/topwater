# Exact stationary pair derivative: bounded next experiment

Status: sample-source refactor, contracted forward gate, and sample-local
reverse prototype implemented only in the isolated development worktree. The existing
q3/q6/q12 derivative controls qualify the discrete gradient identities, not the
physical resolution or motion accuracy.

The next candidate is a shared sample primitive, followed by a sample-local
analytic pullback. A contracted forward scalar is the cheapest correctness
gate. Wider forward lanes alone have already encountered register pressure.

## Exact scalar and ownership

Write `z=V_sym^-1 H phi`, evaluated at the CURRENT geometry and potential,
`c_i=(C phi)_i`, `A=(c_i dot c_j-z_i*z_j)/2`,
`g=sum_c phi_jc*grad_jc`, `m=sum_c phi_jc/3`, and
`B=m+g dot (x-center_j-h*n_j)`. For each off-diagonal target sample the existing
`cuda_pair.cuh::energy_pair` contracts exactly to

```
rho*area_i*w/(4*pi) * [
  sum_e (A*d_e-z_i*h*sign_j*(outward_e dot g))*log_e
  + angle*(z_i*sign_j*B-A*h)] .
```

Here `d_e=-r_e dot outward_e`. This follows by contracting the existing three
mixed-trace columns before forming them, not by changing quadrature or a
boundary equation. For `i==j`, retain the existing explicit `h=0, angle=0`
program branch and smooth double-layer zero. Add the jump exactly once per
PAIR: `rho*z_i*area_i*sum_c(phi_ic)/6`. Both target/source local-coordinate
derivatives on a self pair accumulate to the same physical vertices.

The auxiliary `z` and nodal `phi` are stationary constants inside the gradient.
Both current curls in `A` and the current P1 gradients in `g` remain active.
They depend on the same vertices but are distinct feature dependencies; neither
can be dropped. The symmetric coefficient multiplying raw `dV_ij` makes the
directed raw-pair contraction equal to contraction with `d(V+V.T)/2`.

## Small shared-source refactor

Extract the EXISTING quadrature-loop body to one templated `energy_sample`
function returning `single`, `double`, and `moment[3]` before target weighting.
Keep the current `energy_pair` signature as a raw V/H adapter and its self jump.
A second scalar adapter uses

```
area_i*w * [A*single + z_i*(m*double + g dot moment)]
```

and the same jump. The edge logs, solid angle, source geometry and branch
formulas have ONE owner in `energy_sample`; the two adapters only perform the
different linear contractions. Existing primal, directional, shape and Taylor
cache APIs remain compatible. First differentiate this scalar with the existing
one-lane Dual and compare against the existing raw V/H stationary scalar.
This is a validation gate, not a promised large speedup.

## Pullback after the gate

Implement a specific analytic pullback of that one sample, not a dynamic AD
tape. The primal sample runs once. Its log/angle/sqrt values and branch choices
are reused while propagating their scalar adjoints. Keep the workspace to one
sample, stream the three edge pullbacks, and reduce sample contributions before
advancing quadrature. Do not retain a tape for the whole quadrature loop or an
`nf*nf*nq` array.

Useful contracted seeds, with `kappa=rho*area_i*w/(4*pi)`, are

```
K_e = A*d_e-z_i*h*sign*(outward_e dot g)
T   = z_i*sign*B-A*h
bar(log_e) = kappa*K_e
bar(angle) = kappa*T
bar(A) = kappa*(sum_e d_e*log_e-h*angle)
bar(g) = kappa*z_i*sign*((x-center-h*n)*angle-h*sum_e outward_e*log_e)
```

The remaining direct feature adjoints follow from these same local products.
Propagate `bar(A)` through both curls, `bar(g)` through source P1 gradients,
sample-position adjoints through target barycentric coordinates, and geometry
features through the existing analytic coordinate Jacobians or a small direct
triangle pullback. The first prototype should preserve those Jacobian owners;
replacing their chain rule is a separate optimization.

Retaining eighteen coordinate accumulators already costs eighteen doubles per
thread. Keeping all feature adjoints PLUS all primal temporaries live can again
spill. Use scoped edge primitives and promptly reduce their feature adjoints;
do not infer occupancy from source array counts. A two-kernel feature-adjoint
reduction is a fallback if the bounded one-pass kernel spills badly, not the
first implementation: it adds global reductions/atomics and another ownership
interface.

For current default width3 the nonlinear sample evaluation repeats six times
per pair. A successful reverse pass needs one such primal evaluation plus its
pullback. This removes repeated transcendental evaluations; it is an operation
count, not a sixfold timing claim. Gradient packing, memory traffic, atomics,
factor solves and low GPU occupancy remain outside that count.

## Exact branch requirements

- Self `h=0` and `angle=0` are program constants, not a generic formula limit.
- Match existing Dual `fabs(0)` and `fmax` tie conventions exactly.
- Keep source normal-sign ownership and the complete area derivative.
- Initially match the current atan2 pullback `(x*dy-y*dx)/(x*x+y*y)`;
  numerical stabilization would be a distinct, explicitly tested change.
- No float32, fast-math, row/gauge projection, velocity cap or geometry freeze.

## Cheapest falsifier and stop conditions

1. Refactor only the sample source; reproduce raw V/H, current E, and the
   existing directional action on a small closed control and an actual saved
   generic boat state. Stop if same-quadrature differences exceed ordinary
   floating-point reassociation relative to existing derivative controls.
2. One-lane contracted scalar must reproduce the raw-pair stationary value and
   all eighteen local derivatives, including self, adjacent/coplanar and
   slender near-contact pairs. Then compare the complete current gradient on
   archived canonical q3/q6/q12 states and arbitrary-direction contractions.
   FD alone is insufficient on the most slender pairs.
3. Only then implement the sample-local pullback. Compare the same full
   gradients and all fifteen parent / six body action contractions. Measure
   WHOLE gradient call, compilation separately, and kernel registers/local
   storage with serialized old/new/new/old calls. Require a repeatable material
   multiplicative gain (proposed screening bar: at least 2x whole-gradient on
   the target q3 boat case) before integration. If resource pressure erases the
   gain, stop and report the profile rather than widen lanes again.

## Toolchain alternatives

Enzyme's documented C++ path is an LLVM/Clang plugin, and its CUDA guide uses
Clang plugin compilation. NVRTC exposes CUDA C++ runtime compilation and
PTX/cubin generation; its documented options provide no Clang/LLVM plugin hook.
Therefore Enzyme is a separate toolchain experiment, not a drop-in replacement
inside this CuPy/NVRTC path. This is an inference from the published interfaces,
not a claim that Enzyme cannot differentiate CUDA. The CUDA guide itself is
dated; do not use its historical version limits as current support claims.

Primary references: [Enzyme installation](https://enzyme.mit.edu/Installation/),
[Enzyme CUDA guide](https://enzyme.mit.edu/getting_started/CUDAGuide/),
[NVIDIA NVRTC documentation](https://docs.nvidia.com/cuda/nvrtc/index.html).

Local evidence: `artifacts/boat-surface-gradient-batch3-q12-964ms.json` records
width1 39.69/46.88 s versus width3 34.73/36.61 s, both 255 registers, local frame
72 versus 504 bytes. Those compiler attributes are not an occupancy measurement;
the four calls do not establish a stable hardware speed ratio. Recent bounded
q3 current-state audits show the full gradient still dominates an exact current
field-rate evaluation, making it the relevant next target.

## Stage-one result, 2026-10-22

Five GPU primitive cases passed (19.83 s), including q3/q6 contracted gradients
and retained raw-kernel-cache compatibility. A separate actual 472-face q3
boat audit loads the SHA-verified pre-refactor producer archive as its control:
new raw E relative difference 7.06e-17, maximum difference 2.22e-16; raw H
maximum difference 5.55e-17. The contracted width1 gradient differs from frozen
raw width3 by 1.565e-13 relative /4.03e-11 N maximum. All fifteen parent psi-rate
differences are <=2.22e-15 and six body pose-load differences <=1.64e-11 N/Nm.

This validates the shared source and scalar contraction on that bounded case;
it does not establish full eighteen independent local-pair derivatives or
q12/contact coverage. Width1 contracted timings 0.980/0.886 s versus frozen
width3 0.760/0.672 s show no speed adoption case. Both use 255 registers;
reported local frame is 80 versus 504 bytes. The immutable report and producer
archive are `artifacts/scalar-gradient-q3-generic472.json` and its `-source`
directory in the isolated worktree. Reverse work must retain this control.

## Stage-two bounded reverse result

The shared primitive now accepts a fixed typed observer; a null observer leaves
raw/forward paths unchanged. The reverse observer records one sample's primal
values. A specific hand-derived pullback streams edge feature adjoints through
the existing geometry Jacobians into nine target and nine source coordinate
accumulators. No runtime tape or second numerical sample formula was added.

Five primitive GPU cases passed (21.44 s), including q3/q6 reverse/full-forward
gradients and independent eighteen-lane diagnostics. On the same actual
472-face q3 boat, the frozen PRE-observer producer comparison gives raw E
relative difference 6.46e-17 and raw H maximum difference 5.55e-17. Reverse
gradient relative difference is <=1.59e-13, maximum <=4.13e-11 N; random-direction
contraction difference <=8.17e-11 N. All fifteen parent psi-rate differences
are <=4.89e-15 and six body pose-load differences <=2.19e-11 N/Nm.

An explicitly opt-in compile macro writes only eight selected actual pair
outputs (18 virtual derivatives plus value); production output remains global
vertex accumulation. Body/free self, contact-edge adjacency, coplanar wall
adjacency and near-slender pairs in both directions give maximum difference
1.14e-13, relative <=1.38e-15. No global pair-gradient tensor is allocated.

Serialized old/new/new/old whole-gradient calls were 0.714/0.137/0.0938/0.686 s.
The corresponding observed paired ratios are 5.2 and 7.3; four calls on one
state do not establish a stable hardware speed guarantee. Compilation/startup
can affect first calls. Reverse still uses 254 registers and a reported 336-byte
local frame, versus old width3 255 registers/504 bytes. CPU geometry packing,
auxiliary reuse/solve and GPU synchronization are included in these calls.

Report and hashed producer archive: `artifacts/reverse-gradient-q3-generic472.json`
and its `-source` directory in the isolated worktree. There is no driver
integration, motion or physical-resolution qualification.

## Harder canonical controls

The reverse producer was checkpointed separately at `c6c6080`, after the
source/forward checkpoint `3eaeda6`. Subsequent serialized controls used the
same SHA-verified PRE-observer producer archive, current physical potential,
all fifteen parent psi-rate / six body pose-load contractions, and eight tiny
selected-pair value plus eighteen-virtual-derivative diagnostics. The late
state is 0.9644647283463897 s on explicit chart 33: its entire input hash
manifest, chart token, and bitwise fine geometry replay were checked.

| Case | Maximum gradient relative difference | Maximum body pose-load difference N/Nm | Maximum parent psi-rate difference | Whole-gradient old/new/new/old seconds |
|---|---:|---:|---:|---|
| Generic 472 faces, q6 | 4.21e-13 | 2.55e-11 | 6.22e-15 | 2.253 / .284 / .259 / 2.140 |
| Late contact 1088 faces, q3 | 4.18e-10 | 8.52e-10 | 1.78e-13 | 2.408 / .342 / .200 / 1.368 |
| Late contact 1088 faces, q6 | 1.66e-10 | 1.85e-9 | 1.09e-13 | 9.294 / 1.069 / 1.231 / 10.684 |

Late-state coordinate-gradient absolute differences reach 3.49e-7 N at q3 and
1.03e-7 N at q6; do not describe those entire arrays as bitwise or machine-epsilon
matches. The physical parent/body contractions above are substantially closer.
Selected-pair value/derivative differences remain <=4.55e-13 absolute and
<=4.16e-15 relative; the selected slender face quality `area/max(edge_length²)`
is 2.615e-4. Frozen raw E relative differences are <=7.83e-17 and raw H maximum
differences <=8.33e-17 across these three controls.

These are component timings and derivative controls on two archived physical
states, not trajectory speed/accuracy qualification. Startup/clock variability
is visible, especially the first late q3 raw call. Canonical q12 remains
untested. No attempt was made to tune register limits or repair rows/gauge.
The exact reverse tests no longer depend on the earlier approximate raw-kernel
Taylor-cache prototype; that unrelated compatibility check remains historical
evidence in its own preserved test/report set.

Report hashes, source archive manifests and input hashes are recorded in
`docs/research/boat-stationary-pair-reverse.json`. Complete reports and matching
source archives are also copied into the root worktree's `artifacts` directory.

# Exact current-operator action for pressure offsets

Target only the canonical 24 ft / 2,500 kg boat and selected observed SoCal sea.
This is the next bounded experiment, not an adopted solver or speed claim.

Keep the center's existing full energy matrix, mixed solve and six acceleration
potentials. Each pressure offset needs one potential solve and no acceleration
columns. Use the center mixed Cholesky only as a preconditioner; evaluate every
offset action and final residual using that offset's current geometry.

The existing CUDA assembly uses `E = W + H.T V^-1 H`. Its exact action is

```text
E X = sum_alpha C_alpha.T V (C_alpha X) + H.T V^-1 (H X)
```

Here `C_alpha` is the existing P1 face-curl map. Retain the same symmetrized
single-layer `V`, its current Cholesky, mixed trace `H`, face curls, winding and
uniform q6 quadrature. Avoid forming `W`, `L^-1 H`, the Gram matrix, full offset
`E`, and the offset mixed Cholesky. The current pair kernel and current `V`
factor remain necessary. The shape-force auxiliary method's order-three
ownership guard must stay intact; it does not own this q6 action.

Extract current boundary data preparation once from `boat_world_field.py`:
fixed/unknown node partitions, parent map and its partition-of-unity check,
body load, actual-time incident outer load, and full free-trace lifting `g`.
Share potential/reaction reconstruction with the direct path. No general
operator hierarchy is needed. Offset potential obeys

```text
E_II phi_I = (right - E g)_I;  g_F = S psi, g_I = 0.
reaction = E phi - right.
```

Use only matching explicit chart/node/unknown ownership for the first test.
Recompute each offset's map, load and incident forcing. No additional gauge
pin, flux projection, eigenvalue shift or incident-field subtraction is allowed.
Preserve fine free reactions before fine-mass inversion and nonlinear Bernoulli
projection wherever that output is consumed. Pressure's ALE correction remains
in its existing owner.

First compare the action with explicit multiplication on the actual potential,
six auxiliary fields, constants and independent vectors. Then compare both
offset potentials and all six pressure loads at h and h/2, with correct incident
times and the same derivative representation. A residual alone cannot qualify
pressure because differencing divides field error by h.

For an illustrative 1,200 faces / 602 vertices, avoiding triangular formation
and Gram multiplication removes about 1.30 billion floating-point operations
per offset. One single-column action is about `8*Nf^2 + 4*Nf*Nv`, or 14.4 million
operations before sparse curl work. The arithmetic crossover is about 90
actions, not a predicted GPU speedup. Seven-column or 21-mode response
construction is unnecessary for this first offset-only falsifier.

Time pair assembly, current V factor, iteration, final residual and complete
pressure separately; include unchanged center work. Retain this organization
only if complete current-field/pressure work is faster at unchanged output
accuracy. If pair integration or iteration erases the saving, reject it and
retain the result as evidence for the next architecture choice.

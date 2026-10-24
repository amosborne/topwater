# Reference study 002: unqualified

Source checkpoint 6fae539; full source hashes and contract/input identities are
in the [complete evidence bundle](../validation/evidence/reference-study-002.zip).
It contains both reference/freeze JSON and NPZ files, both producer archives,
exact input JSON, prelaunch memory/physical-equality evidence and live memory
observations. [Summary and bundle SHA256](../validation/evidence/reference-study-002-summary.json).

All four graded meshes and all132 physics calls completed. No elapsed/resource
failure occurred. The original physical arrays and input JSON remained exact.
Every convergence family remains unqualified under the unchanged contract.

| Family | Increment failures | Trend failures | Entries |
|---|---:|---:|---:|
| A-infinity |17|4|36|
| finite A |31|19|108|
| B/omega |22|12|108|
| diffraction |6|8|54|
| Froude-Krylov |0|3|54|
| total excitation |3|9|54|
| complex RAO |10|19|54|

This substantially improves the rejected midpoint family but does not qualify
the reference. Latest A-infinity increments still reach6.27 times budget.
Froude-Krylov's three trend failures are independent of the BEM/Green function:
Capytaine integrates analytic incident pressure at panel centers. Thus mesh
redistribution/centroid integration contributes non-smooth refinement behavior.

## Next bounded diagnostic: cyclic vertex order

Installed Capytaine3.0.0 meshes/geometry.py computes each face's characteristic
radius using only its first vertex's distance to the centroid. Equivalent
cyclic orderings can therefore change radii supplied to the native kernel.
This source observation does not yet prove a significant response error.
Independent reviewer reference_adversarial_review confirmed the observation
and recommended a bounded representation-invariance diagnostic before another
mesh-family rewrite.

Scope: tests.validation.diagnose_order uses the frozen2464-panel mesh from
study002. Compare original order, cyclic shifts1/2 and farthest-vertex-first.
Verify unchanged centers, areas, normals and winding; retain radii and all36
radiation entries at infinity and0.1Hz (48 physics calls). Predeclare the same
entrywise normalized budget,1e-6+0.001*abs(farthest-first response), to judge
material representation sensitivity. No excitation/RAO/qualification claim.
Archive installed inspected dependency sources alongside its producers.

Stop after these48 calls or a declared validity/resource failure. If significant
sensitivity is demonstrated, review a conservative cyclic-normalization control
before adoption, then rerun the full convergence protocol. A fixed-topology
refinement family remains a separate possible correction for spatial and
centroid-integration error. No case dropping, tolerance change, Richardson
extrapolation, new runtime formulation or optimization is authorized.

## Resource observation limitation

The final monitor attached after PID31160 had already exited. Last saved live
observation (start of the final frequency) reports2.950GB peak working set and
4.384GB peak committed memory. Those are lower bounds on the full-run peaks,
not a final peak measurement. Do not use them as a certified resource bound.

## Completed vertex-order diagnostic

All 48 declared radiation calls completed. Equivalent cyclic orderings differ
from the farthest-first convention by up to 280.27 reference budgets for
A-infinity, 76.52 for finite A and 16.59 for B/omega. The largest differences
occur in heave/roll couplings that should vanish by reflection symmetry.
Original ordering also differs from farthest-first by 1.103 budgets for one
A-infinity entry. Centers, areas and normals passed roundoff invariance checks.

[Complete diagnostic evidence](../validation/evidence/vertex-order-diagnostic-001.zip)
includes all matrices, radii, triangle variants, producer archive and inspected
installed dependency sources. Its bound freeze is retained in the study002
bundle. [Summary and SHA256](../validation/evidence/vertex-order-diagnostic-001-summary.json).
This is a completed diagnostic with a failed representation-invariance gate,
not a qualified reference.

The next contract revision freezes every reference triangle in an
orientation-preserving farthest-vertex-first cyclic order. This makes the
first-vertex characteristic radius the maximum vertex-to-centroid distance
without modifying Capytaine or any physical facet. The geometry test checks
cyclic input invariance, winding/area moments and the actual maximum radius.
All 33 gate/geometry tests pass. This control corrects a demonstrated numerical
representation sensitivity; it does not resolve the remaining spatial or
centroid-integration convergence failures. No full known-failing mesh family
is rerun solely to demonstrate this local correction.

Independent reviewer approved adoption at source checkpoint d160eae. The
[fresh freeze/control bundle](../validation/evidence/frozen-graded-003.zip) and
[proof summary](../validation/evidence/frozen-graded-003-summary.json) show exact
equality of the original physical arrays and input JSON, and verify Capytaine's
actual radii equal maximum vertex distances within 3.1e-16 m on every panel.
This evidence does not imply new hydrodynamic convergence. The approved next
checkpoint is design/geometry preflight of a fixed-topology graded family,
before another complete hydrodynamic study. See checkpoint-handoff.md.

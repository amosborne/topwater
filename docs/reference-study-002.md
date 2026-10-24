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

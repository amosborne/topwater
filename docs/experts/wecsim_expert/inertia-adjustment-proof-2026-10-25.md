# Simscape inertia and added-mass adjustment evidence

The active rigid-body correction now has measured before/after evidence. This
is a Simscape test, unrelated to the nonworking Topwater-native solver. The
held fork already includes the ordering correction; no new pin was promoted.

A 7 kg body uses rotational inertia
`[10 1.1 -1.3; 1.1 12 1.7; -1.3 1.7 14] kg*m^2`.
Under the same 1 N*m roll torque from rest, with gravity disabled:

| Initial angular acceleration (rad/s^2) | Release ordering | Corrected / expected |
|---|---:|---:|
| Roll | 0.104094244641134 | 0.102615384041695 |
| Pitch | -0.015746134614160 | -0.010944563702830 |
| Yaw | 0.010903090436361 | 0.010857554110644 |

Maximum error falls from 0.004801570911330492 to 1.3877787807814457e-17.
Eight layout controls pass. Evidence: `artifacts/workflow/rigid-inertia-trial-attempt3/`.

Six independent wrench columns identify each native 6x6 matrix from measured
accelerations. Twelve simulations cover both orderings and two independent
components: dry and production-adjusted. The actual residual is then added.

| Matrix check | Maximum entry error before | After |
|---|---:|---:|
| Native dry Solid versus M | 0.6000000000000012 | 5.329070518200751e-15 |
| Native adjusted Solid plus residual versus M+A | 0.5849999999999986 | 3.552713678800501e-15 |

Public products and residuals are unchanged between modes. Residual nonsymmetry
is retained (norm 0.017029386365926404). Evidence:
`artifacts/workflow/rigid-inertia-matrix-attempt1/`, including twelve per-basis
MAT measurements, report, acceptance and terminal receipt. Independent source
and measured-evidence review accepted these narrow claims. Frozen identities
and exact process termination were verified; idle read back and leases released.
Failed earlier expectation/link trials remain preserved.

The fixture seeds synthetic added mass through `constAddedMassAndDamping`, then
calls production `adjustMassMatrix`. Residual composition is **algebraic**. This
does not test actual H5/CIC loading, dynamic residual forces, current boat impact,
full motion, hydrodynamic accuracy, Drag or Flex. The corrected total matches M+A
for these fixtures without loss or double counting; wider claims remain unproved.

Added mass is water inertia, not viscous damping. The example already uses added
mass and radiation memory; viscous damping stays off. No damping was introduced.

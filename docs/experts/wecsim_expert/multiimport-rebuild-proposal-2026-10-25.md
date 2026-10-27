# Queued aggregation and dispatch rebuild

Read-only numerical-methods review for own-fork issue
[#2](../history/wecsim-retired-workflows-2026-10-27.md#issue-2). This is a proposal, not an
implementation or numerical allocation. Inertia authoring remains the current
WEC expert assignment.

Prefer one aggregate `hydroForce.hf#.userDefinedFe` history and one active
imported-excitation branch. Keep the wave array intact. Evaluate each imported
signal with its own direction, then sum freshly prepared arrays on the common
simulation output grid. Different raw-file time grids are acceptable when each
record can be evaluated consistently; incompatible output grids must fail.
Reuse the reviewed history/ramp/zero-extension semantics from PR #8 with an
explicit dependency, instead of duplicating the historical combined fix.

The release already rejects mixed wave types in `initializeWecSim.m`, lines
322-327. Moving that check earlier retains the public contract. The same file
rejects imported elevation with nonlinear hydrodynamics, Morison forces and
variable hydrodynamics. Those combinations remain outside scope. Native
imported-source expressions read `hf1`; no new variable-hydrodynamics support
follows from this proposal.

The branch reconciler must touch only recognized native numbered clones and
their exact positive sum ports. Check source ownership and wiring before
removal; preserve unrelated inputs/custom blocks and fail clearly on ambiguous
modified clones. Repeat initialization must be idempotent. Imported/regular
transitions must restore the appropriate native branch topology.

The old dispatch receipt is graph evidence only: its
`body_simulated:false` does not establish the applied force. The new gate must
simulate the actual copied native excitation parent and capture its downstream
sum/applied-force bus. Two distinct analytic six-channel signals must produce
`F1+F2` once, excluding `F1`, `2*F1` and `2*(F1+F2)`. Exercise reversed order,
repeat initialization, one/two/one imported components, imported/regular
transitions, grid compatibility, protected custom wiring and an unchanged regular
control. No Topwater force-table rebinding is permitted. Full body dynamics are
not needed to test this signal boundary.

If Flex dispatch changes, test its actual excitation parent with more than six
DOFs, or explicitly defer/reject unsupported Flex behavior. Isolated signal
boundary tests do not qualify full Flex dynamics.

Inspected base: official v7.1.0 `a23d2993807410fe78cc63546c1ba9d81912024f`.
Source identities: `bodyClass.m` blob `b1dee271ab7160a6af5cb5ee319ca07af73dfc63`,
initializer `f425407bc76860996fba759e62c589b009cd6849`, body library
`34655aed0de2eb10a4033d5f580293180379ce91`. Relevant archive parents are
`system_8_4071.xml` (rigid) and `system_15_5103.xml` (Flex); imported sources
are in `system_8_4071_6338.xml` and `system_15_5196.xml`. Old `2f43f54` and
`f962103` were reviewed as provenance, not automatically replayed. No files,
models, numerical jobs or upstream state were changed by the reviewer.

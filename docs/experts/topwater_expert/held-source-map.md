# Topwater held-state evidence map

Setup snapshot: 2026-10-24. Native custom activation: agent type
`topwater_expert`, actual agent ID `/root/topwater_setup`. The runtime supplied
the Topwater developer role instructions; reading the TOML is corroborating
configuration evidence, not the basis for claiming native activation. Model and
reasoning settings were inherited without overrides.

## Scope and evidence limit

This is a documentation-only map. No implementation source, tests, MATLAB,
simulation, numerical reproduction or performance experiment was executed.
Only this expert's knowledge files were changed. Current code behavior and
installed dependency versions have not been independently revalidated here.
Durable context is these files, not permanent model learning.

The runtime role requires bounded assignments, coordinator-owned acceptance and
serial numerical execution leases, immutable producers/inputs, independent
review of physical corrections, preservation of existing work/failures and no
Topwater implementation during the external-reference phase. All are recognized.

## Governing equations and ownership recorded in project evidence

The documented selected implementation is the P0 Rankine-panel core extracted
from `4e8f511`, on `baseline/fluid-body`; that extraction is provenance, not a
claim about current source identity. `docs/validation.md` maps
`topwater.inertia.instantaneous_added_mass` to `bem`, `boundary` and `rankine`.
The prior expert record also maps `src/topwater/fluid_boundary.py` and
`src/topwater/pressure.py` to closed geometry and pressure ownership. These paths
were not opened during setup.

`docs/mixed-pressure-protocol.md` records fluid-outward normals, hull generalized
normal `N=[n,(x-CG) cross n]`, Eulerian acceleration potential `chi=phi_t`, free
surface condition `chi=-g*z-|u|^2/2`, and Bernoulli pressure
`p=-rho*(chi+|u|^2/2+g*z)`. Its known hull normal-derivative data includes fluid
velocity/Hessian and rotational terms. Seven mixed solves represent one known
column and six acceleration columns; raw `A=rho*N^T*area*chi_aux` enters
`(M_dry+A)*a=known fluid load+weight-gyroscopic load`. These are documented
formulation requirements, not a new derivation or source audit.

## Evidence and its physical limits

| Record/version | Established in the record | Limit |
|---|---|---|
| `docs/validation.md`, extraction `4e8f511` | Hydrostatics and instantaneous added mass; generic mathematical coverage retained | Explicitly unqualified; no complete coupled wave-response implementation claimed |
| `docs/mixed-pressure-001.md`, pressure `6213d96`, attribution `902c2be` | Closed conforming fluid geometry, seven-column pressure/body coupling, calm pressure algebra and geometry checks; saved generic-test evidence | Both diagnostic artifacts have `passed=false`; no wave/motion/runtime qualification |
| Same pressure record | Calm acceleration at most `3.12e-15` SI and direct pressure-load reconstruction | Calm algebra does not establish fluid-response accuracy; six auxiliary traces were not archived |
| Same pressure record | All 144 retained inertia comparisons reproduce four failed original diagnostic outcomes | Multiple discretization/domain contributions unresolved; no same-mesh or finest-reference qualification |
| `docs/reference-convergence.md` | Study 003 failed original coefficient/trend rules | Superseded campaign remains failed; do not reinterpret under new motion policy |
| `docs/branches.md` | Optimization graph preserved at `archive/optimization-2026-10-23`, `4e8f511` | Historical work is not implementation authorization |

The pressure record's provisional `36.17` wall-s/s estimate excludes derivative
recovery and remaining RHS work. It cannot establish complete runtime. Archived
P1/canoe evidence cannot qualify the canonical 24 ft boat.

## Missing or unqualified capabilities

The records identify complete coupled wave-response RHS/state evolution and
boundary-gradient/Hessian recovery as unfinished. Known-incident representation,
initial compatibility, finite-domain controls, nonlinear/time-step controls,
matched six-motion validation and complete runtime remain unqualified. Existing
recorded geometry/algebra checks cannot substitute for them. Before any later
implementation, recheck actual source identity and current behavior rather than
treating this historical map as timeless.

`docs/time-domain-baseline-goal.md` is explicitly HELD. Its old numerical goals
and `docs/validation.md`'s stale active-goal/coefficient instructions do not
override `docs/capytaine-wecsim-workflow.md`, shared contract, current user setup
scope or native role boundary. External reference loads/trajectories must never
replace independently calculated Topwater fluid physics.

## Exact read set

- `.codex/agents/topwater_expert.toml`
- `docs/experts/topwater_expert/README.md`
- `docs/experts/topwater_expert/current-state.md`
- `docs/experts/shared-contract.md`
- `docs/experts/coordination.md`
- `docs/validation-reset.md`
- `docs/mixed-pressure-001.md`
- `docs/time-domain-baseline-goal.md`
- `docs/reference-convergence.md`
- `docs/mixed-pressure-protocol.md`
- `docs/validation.md`
- `docs/branches.md`
- `docs/capytaine-wecsim-workflow.md`

A bounded `rg` documentation search additionally inspected matching lines in
`docs/fixed-asset-preflight-001.md`, `docs/pressure-discrepancy-protocol.md`,
`docs/reset-review.md` and `docs/experts/setup-instructions.md`; these were not
fully read. No code files or evidence-bundle payloads were opened.

## Stopping checkpoint

Native setup confirmed; implementation remains HELD. No execution lease was
requested or possessed. Next action is to await a bounded coordinator assignment;
the immediate project milestone remains external-reference minimum-mesh 10 s
incoming-wave motion/video followed by motion variation with mesh. No technical
setup blocker was identified for this role.

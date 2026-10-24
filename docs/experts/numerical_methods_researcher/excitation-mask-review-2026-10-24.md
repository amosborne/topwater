# Independent native excitation-mask review, 2026-10-24

Bounded assignment: inspect the actual selected native library mask and imported
branch, derive independent applied-force expectations, and review the owner's
minimal correction and regression. No implementation or numerical launch.

Earlier retained-history/component-sum regressions exercised native
`bodyClass.hydroForcePre` preprocessing. They did not exercise the compiled
Simulink excitation mask; their signoff cannot establish applied-force behavior.

## Actual library source observations

Read the selected
`artifacts/dependencies/topwater-wec-sim/source/lib/WEC-Sim/WECSim_Lib_Body_Elements.slx`
as a ZIP/XML archive, without editing it or loading a model.

- `simulink/systems/system_8_4071.xml`, block SID `8:4071:6252`, parent
  `Wave Diffraction and Excitation Force Calculation`: mask initialization sets
  `nWaves=size(waveGroup,2)`, clones a `Linear Wave Excitation Force Variant
  Subsystem` for each additional component, and adds positive input ports to
  `Add`.
- `simulink/systems/system_8_4071_6338.xml`, imported-elevation branch:
  FromWorkspace reads
  `[Nwave.waveAmpTime(:,1),body.hydroForce.hf1.userDefinedFe]`. Only the time column
  depends on `Nwave`; the force array is the same precomputed total for every
  cloned imported branch.
- A second masked path SID `15:5109` has the analogous clone/read pattern. The
  owner must identify which actual path the first body uses and bound any
  correction's claimed coverage accordingly.

Independent expectation: homogeneous imported components are already combined
in preprocessing, so the total must be applied ONCE. With three native cloned
branches each reading that same total, the source predicts applied force
`3*Ftotal`, not `Ftotal`. This is a source-derived prediction pending executed
actual-mask reproduction.

## Minimal correction and regression expectations

Using one native excitation branch for a homogeneous imported group is a valid
fresh-mask policy; ordinary regular/irregular per-component branches must retain
their prior semantics. Reject mixed wave-type groups explicitly as unsupported
by the native case contract.

Critically, changing only `nWaves=1` may prevent new duplicates without removing
clones and Add ports already persisted by earlier initialization. Test repeated
initialization and three-ordinary-waves to three-imported-waves reconfiguration,
or fail explicitly on an unsupported persisted state. Do not silently leave
duplicate applied forces.

The decisive regression must instantiate the actual selected library mask and
branch, using proper Simulink APIs, with a known six-component force history.
Assert exactly one applied total for three imports and unchanged ordinary
multiple-wave behavior; test row/column wave-group handling if the correction
claims that support. Preprocessing-only sums or hand-constructed substitute
force paths cannot prove the integration correction.

## Actual native-model inventory

Owner stock attempt four compiled/updated the genuine native model successfully,
then its shallow branch-count query failed. That failed query proves neither
presence nor absence of duplication. Subsequent read-only all-mask/all-link
inventory of the retained stopped updated model is authoritative:
`artifacts/reference-imported-dispatch-stock-attempt4/structure.json` records
three actual linear-excitation branches, three identical imported total-force
expressions, `waveGroup=waves`, three prepared waves and Add inputs `+++++`.

The five Add inputs do not mean five imported copies: selected XML shows the
original three ports are linear excitation, nonlinear FK and second-order
excitation. The two cloned linear branches append ports four and five. For the
first linear case, nonlinear FK and second-order terms are disabled.

Subsequent authoritative `wiring.json` and inspection source establish the
actual three-copy defect. Each child output resolves independently to Add ports
1, 4 and 5; all signs are positive. `Simulink.match.activeVariants` finds three
active imported-history blocks, with `sv_udfWaves.Condition` (`typeNum>=40`)
true, all reading the same prepared total-force array. Native setup gives the
same simulation times to each wave; corrected preprocessing also requires exact
common output times. The actual updated graph therefore applies three copies
of the total, where one is required. This is an algebraic graph proof, not a
body-motion simulation.

Recomputed `inspect_dispatch_structure.m` hash
`0a31ad2ff4ac83034043014bd2f3f7011dcf88c4776e6011622f3b2eaa82f6f7`
matches the frozen wiring manifest. The log and execution ledger record the
read-only inspection's terminal exit zero. Earlier failed shallow queries and
mask-initialization attempts remain historical failures, not defect proof.

Status: correction source and fixed-regression review pending. No body result
may rely on this path until resolved.

## Prospective correction source review

Read-only review passes the proposed selected Rigid Body correction:
`source/functions/simulink/model/importedExcitationBranches.m` in the maintained
fork returns one branch for a homogeneous imported group, explicitly rejects
mixed imported/non-imported groups, and preserves the native `size(group,2)`
count for other modes. It detects persisted numbered native clones, validates
their wave references and unique positive Add destinations, then removes their
output lines, blocks and corresponding sum ports. An already-clean imported
group has no clones to remove, so repeated initialization is idempotent.

`reference/wecsim/fix_imported_dispatch.m` edits the selected authored library
mask through `Simulink.Mask` API, changing only the count header to call that
helper and retaining the native non-imported branch-creation loop verbatim.
This is not a string replacement of a generated case model. The separate
SID-15 library path is untouched and outside this correction signoff.

Executed fresh, repeated, ordinary-three-wave to imported-three-wave transition,
ordinary-mode preservation and mixed-group rejection tests on the actual mask
remain required. Source review alone is not final physical reliance signoff.

The proposed `check_imported_dispatch.m` regression source passes independent
review: it builds and updates a genuine fresh native case; verifies one active
imported history connected to a positive sum port; changes the group expression
to force repeated initialization; switches to three ordinary `regularCIC`
objects and verifies three positive branches with no active import; restores
three imports and verifies stale clones/sum ports are removed; explicitly tests
mixed rejection and final restoration, with precomputed total force unchanged.
The ordinary-mode check verifies replication topology, not ordinary-wave
hydrodynamic accuracy. Body time integration is not part of this regression.

The first executed fixed suite passed fresh and repeated imported-case updates,
then failed in the ordinary-wave test fixture: numeric `Inf` passed as the BEM
water-depth argument produced an invalid infinite visualization frame. The
selected `waveClass.setWaterDepth` source explicitly recognizes the H5 string
`'infinite'`, sets deep-water physics and a finite 200 m visualization depth.
The owner changed only that fixture to pass its actual H5 water-depth token as
native `initializeWecSim` does. This fixture correction passes source review;
the failed attempt is preserved and does not establish the remaining tests.

## Final independent executed mask signoff

PASS for the selected Rigid Body linear imported-elevation dispatch correction.
Independently inspected
`artifacts/reference-imported-dispatch-fixed-attempt2/result.json`, its execution
log and exact terminal ledger: fresh, repeated and final imported groups each
have one active imported history on positive Add port one; ordinary three-wave
mode retains three positive branches on ports 1, 4 and 5 with no active imported
history; returning to imports removes stale clones and restores `+++`; mixed
groups explicitly fail; the prepared total-force array remains unchanged.
Exact session 40609 terminated exit zero, owner model closed, no time integration.
The lease was absent at this review.

Recomputed all fourteen current frozen-input hashes from the attempt-two
manifest, plus archived helper/test/library snapshot hashes: all match. Relevant
identities are helper
`20ce3154543e8fa03bc7e2fe9818358eb3be5eb3e5d77e9a96642314e0078a03`, test
`77320e1ba7273e2aa7e4e86776d134a058b27109096bb578f7aac56c81621b59`, and library
`9251a61c5a23d20eb8a4cbe98a278054983e84defc20b1446e74e1b3db64bc87`.

Independently compared all 146 SLX system XML trees against fork baseline
`fcd7c365a1201bbb56ce7f571a5baf9b2a94d89b`, canonicalizing parameter ordering
changed by the native save operation. Exactly one system changes semantically,
`system_8_4071.xml`, and its sole semantic change is the expected count-header
replacement. All other system semantics and the imported-history branch remain
unchanged. The new helper is separately source-reviewed above.

The author and reviewer are distinct. Physical reliance on this selected-mask
correction is permitted after the coordinator freezes its maintained-fork pin.
This signoff does not establish body dynamics or broadly certify other body
library paths/wave modes. The first actual body run must still demonstrate exact
rest initial state and compare runtime excitation against the prepared total
history applied ONCE. Prior preprocessing-only tests and this mask graph test
remain separately scoped evidence.

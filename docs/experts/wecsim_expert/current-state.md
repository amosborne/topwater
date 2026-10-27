# Current state: WEC-Sim / BEMIO expert

Updated 2026-10-27 UTC. Authority: [handoff](../../checkpoint-handoff.md),
[physical contract](../shared-contract.md), [coordination](../coordination.md).
Root owns acceptance, integration and serial execution. Consult the
[execution ledger](../execution-state.json) for live/terminal status; this page
confers no execution permission.

## Active reference and adoption

The user squash-merged inertia PR9 and spectral PR14. Our release-based reference
is now **6b774228e1f6b31dee473e626204b0ef65e5ec38**, including merged PR6/7/9/14.
Its clean detached checkout is `artifacts/dependencies/wecsim-reference-validation`.
The isolated Topwater standard-input adapter is **9539192bb733663584483f33ef2fb3de273d6fcf**,
branch `work/standard-wecsim-input`, at `artifacts/worktrees/standard-wecsim-input`.
It has not been adopted into the active implementation branch.

The route uses WEC-Sim's EXISTING public `bodyClass(hydroStruct)` input. Official
BEMIO import/normalization and `rebuildHydroStruct` prepare ordinary excitation
fields; consumed Ainf/CIC/restoring inputs remain unchanged. Unused raw A/B/T
are omitted rather than fabricated. No custom coefficient-view/accessor or BEMIO
writer extension remains in the merged fork. WEC-Sim computes forces and motion;
Topwater supplies physical waves/coefficient inputs, never applied force histories.

**New-reference force/motion validation is accepted; default pin9952634 remains held.**
Fresh field/delivery and adopted-default readiness are still pending.
The user authorizes completing validation/adoption, then the
[fork cleanup plan](fork-cleanup-plan-2026-10-27.md). Pin/source/output publication
and cleanup are root-owned later steps, not consequences of source approval.

## Evidence for this standard-input route

- The tiny public-input discriminator passed all three low/off-grid/seam and
  missing/empty-unused-A/B cases. Consumed CIC/Ainf/mass/restoring were exact.
  [Independent acceptance](../../../artifacts/workflow/standard-input-review-v1/acceptance.json).
- Narrowed PR14 passed ten generic native regressions, including the native force
  block. [Acceptance](../../../artifacts/workflow/standard-spectrum-unit-review-v1/acceptance.json).
- Exact6b/953 passed **38 focused Python checks and21 MATLAB cases**
  (13 validators +8 handoff), with no failed/incomplete native cases.
  [Independent acceptance0cffbbc3](../../../artifacts/workflow/reference-handoff-review-v1/acceptance.json)
  binds actual source/import paths and measured evidence. Sparse-bin/CG-CB fixture
  warnings do not qualify angular interpolation or physical accuracy.
- **Full canonical6b/953 motion is independently accepted:**
  [receipt286f3a82](../../../artifacts/workflow/reference-adoption-motion-review-v1/acceptance.json)
  binds [case3 reportec4d8d8d](../../../artifacts/workflow/reference-adoption-case3/report.json).
  All3001x6 force/startup/phase, actual once-only application, rest/negative-history
  certificate, full CIC replay and clock gates pass. Maximum position difference
  is59.758 micrometres and relative SO3 angle0.00156985deg, within1mm/.01deg.
  This is numerical compatibility, not hydrodynamic accuracy. The native allocation
  is terminal with frozen inputs and process/cleanup proof verified.
- Measured controlled runtime1212.407s includes simulation745.503s, native force
  replay305.104s, initialization49.085s and model preparation51.773s; it is not all
  integration time. Stage timings and raw outputs are retained in the case3 report.
  Its final MATLAB log has one coarse-direction coverage warning plus visualization-
  depth/imported-frequency informational messages. Nodal force agreement does not
  establish accuracy of the original angular coefficient grid. Earlier startup/
  guard failures remain preserved under their original case/allocation directories.
- **Fresh TOTAL field assembly is the current root allocation**, using this saved
  accepted motion: no body rerun or new coefficient campaign. Delivery sourceb078f2fe
  is independently reviewed. Same operatorf1f788f9, centroid quadrature and normalized
  producer source240198 remain; raw953 CRLF source528e2e94 is separately bound.
  Old-cache inputs are screened by exact bytes; no old field is reused. Result/PNG,
  full ordinary CLI render and exact-eight-file audit remain pending.
  [Delivery requirements](../../../artifacts/workflow/reference-adoption-delivery-v1/README.md).
- Corrected default-preflight sourcedf6dd55a is independently source-approved,
  **unexecuted**. It requires the adopted main adapter/default6b checkout, owns a
  fresh hidden shared engine, calls actual preflight and requests MATLAB-side
  exit only after PID verification; terminal checks remain mandatory.
  [Default-path check](../../../artifacts/workflow/reference-adoption-preflight-v1/README.md).

## Preserved evidence and limits

The earlier e24 optional-view route's accepted forces, bus, motion and private
8-file delivery are HISTORY, not acceptance of6b/953:
[earlier motion](spectral-motion-qualification-2026-10-27.md),
[earlier end-to-end qualification](spectral-end-to-end-qualification-2026-10-27.md).
Public canonical files remain unchanged. Physical issues31–34 remain unresolved;
mesh/irregular-frequency, field-phase and fixed-mean applicability warnings are
not error bounds or blanket waivers.

[Inertia component proof](inertia-pr9-qualification-2026-10-26.md) covers the active
Rigid Hydro boundary, not Drag/Flex or general hydrodynamic accuracy.
[Engine lifecycle](workflow-engine-lifecycle-2026-10-25.md),
[retired imported-history study](imported-history-rebuild-2026-10-25.md) and
[kernel/rejected-route evidence](kernel-and-inertia-checkpoint-2026-10-26.md)
retain their original versions and policies. Do not resume retired studies from
historical plans. No upstream issue/PR/contact is authorized without explicit user
permission.

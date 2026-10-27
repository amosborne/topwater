# Current state: WEC-Sim / BEMIO expert

Updated 2026-10-27 UTC. Authority: [handoff](../../checkpoint-handoff.md),
[physical contract](../shared-contract.md), [coordination](../coordination.md).
Root owns acceptance, integration and serial execution. Consult the
[execution ledger](../execution-state.json) for live/terminal status; this page
confers no execution permission.

## Adopted reference and workflow

**The default pin is6b774228e1f6b31dee473e626204b0ef65e5ec38.** The user merged
PR6/7/9/14 into the release-based `topwater/reference`. The default checkout is
`artifacts/dependencies/topwater-wec-sim`, clean on that branch; no path override
was used for actual default-path preflight.

Adapter9539192bb733663584483f33ef2fb3de273d6fcf is integrated into the active
`baseline/fluid-body` checkout with its scientific files exact after only CRLF-to-LF normalization. The backend differs
only by the validated6b pin. The source953 worktree and detached6b validation
checkout remain for provenance, not as required runtime overrides.

WEC-Sim computes forces and motion through its EXISTING public
`bodyClass(hydroStruct)` spectral input. Official BEMIO import/normalization and
`rebuildHydroStruct` prepare ordinary excitation fields; consumed Ainf/CIC and
restoring inputs remain unchanged. Unused raw A/B/T are omitted rather than
fabricated. No custom coefficient-view/accessor or BEMIO writer extension remains.
Topwater supplies waves/coefficient data and reconstructs the total surface.

**The accepted eight-file result is now promoted to `examples/simple-boat`.**
Boat/sea are unchanged. Cleanup is the next authorized phase; it is NOT complete.
The legacy995 branch remains until its archive is verified and retirement follows
[the cleanup plan](fork-cleanup-plan-2026-10-27.md). Root owns those actions.

## Accepted evidence

- Tiny public-input discrimination passed three low/off-grid/seam and missing/
  empty-unused-A/B cases; consumed CIC/Ainf/mass/restoring were exact.
  [Acceptance](../../../artifacts/workflow/standard-input-review-v1/acceptance.json).
  Narrowed PR14 also passed ten native regressions, including the real force block.
- Exact6b/953 passed38 Python and21 MATLAB cases (13 validators +8 handoff).
  [Acceptance0cffbbc3](../../../artifacts/workflow/reference-handoff-review-v1/acceptance.json).
- Full canonical3001x6 forces/startup/phase, once-only application, rest/negative-
  history certificate, full CIC replay, clocks and motion passed.
  [Acceptance286f3a82](../../../artifacts/workflow/reference-adoption-motion-review-v1/acceptance.json)
  binds [case3 reportec4d8d8d](../../../artifacts/workflow/reference-adoption-case3/report.json).
  Maximum position difference59.758 micrometres and SO3 angle0.00156985deg meet
  1mm/.01deg. Controlled runtime1212.407s includes simulation745.503s and native
  force replay305.104s; it is not all integration time. Failed prior attempts remain.
- Fresh TOTAL incident/diffraction/radiation field, standard Result/PNG, full
  ordinary CLI render and exact-eight audit passed from the saved NEW motion.
  No body rerun, coefficient campaign or old-field substitution.
  [Private acceptance4622b301](../../../artifacts/workflow/reference-adoption-private-review-v1/acceptance.json)
  binds final report2104ca1f,3001 motion/751 surface samples and751 decoded frames.
  Same operatorf1f788f9,17 panels and centroid quadrature were retained.
- The active Topwater checkout passed38 focused tests in2.69s and real default preflight: clean6b
  checkout, no override, hidden fresh R2025b, all three required licenses, zero
  loaded models and verified owned-server shutdown.
  [Default-path acceptance3712ab33](../../../artifacts/workflow/reference-adoption-main-review-v1/acceptance.json).
  This readiness check performed no further body integration or public solve.
- [Promotion report](../../../artifacts/workflow/reference-adoption-promotion1/report.json)
  verifies all eight canonical files byte-for-byte against the accepted candidate:
  Result2b615b56, PNG275ac4c5, videoe5fce0e4. Physical inputs are unchanged and the
  verified rollback copy of all eight original files is retained.

## Limits and preserved history

This is implementation/numerical compatibility acceptance, not physical accuracy.
Issues31?34 remain unresolved. The native coarse-direction warning is a coverage
heuristic for the repacked nodes; matching forces does not establish genuine
angular coefficient resolution. Mesh/irregular-frequency, field-phase and fixed-
mean applicability warnings are not error bounds or blanket waivers.

Earlier e24 optional-view qualification is historical evidence, not the current
architecture: [motion](spectral-motion-qualification-2026-10-27.md),
[end-to-end](spectral-end-to-end-qualification-2026-10-27.md).
[Inertia proof](inertia-pr9-qualification-2026-10-26.md) covers active Rigid Hydro,
not Drag/Flex or general accuracy. [Engine lifecycle](workflow-engine-lifecycle-2026-10-25.md)
and [retired workflows](../history/wecsim-retired-workflows-2026-10-27.md) retain
original evidence/policies. Do not resume retired studies from historical plans.
Generic fixes remain upstream candidates; no upstream issue/PR/contact is
authorized without explicit user permission.

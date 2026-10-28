# Current state: WEC-Sim / BEMIO expert

Updated 2026-10-27 UTC. Authority: [handoff](../../checkpoint-handoff.md),
[physical contract](../shared-contract.md), [coordination](../coordination.md).
Root owns acceptance, integration and serial execution. Consult the
[execution ledger](../execution-state.json) for live/terminal status; this page
confers no execution permission. Latest numerical allocations are terminal; the
coordinator reports an idle ledger and absent lease at this checkpoint.

## Active packet integration (source checkpoint; not yet validated)

Coordinator authorized the new directory solve contract: positional dwell,
rest/ramp defaults10/5s, dwell100 means total130s, packet/native dt.01s,
surface cadence.08s. Canonical boat/sea carrier arrays stay unchanged; phases
refer directly to midpoint65 in the existing negative-time physical convention.
`_wecsim_setup` calls the utility with negative phase and conjugates positive
phasors; each sideband keeps original heading and its own dispersion. Aggregate
1cm envelope-exceedance budget is allocated proportional to carrier amplitude.
Original carrier operating eligibility is retained; generated sidebands above
it are warned about and retain genuine H coverage, without mesh-accuracy claims.

The approved DC reference shift keeps all computational BEM/CG/K geometry fixed.
Native initial heave−d supplies static fullK-heave-column coupling; physical
motion and surface receive+d once. Surface geometry runs in computational frame;
startup tests physical dynamic+d and private replay restores the world offset.
Native net-load/startup, applied-force, full inertia, CIC and equation-of-motion
audits remain. No Topwater force history or fork change was introduced.

Source-only implementation is being independently reviewed by Capytaine and
numerical-methods agents. Focused tests added for phase/direction/coherence,
aggregate budget, fullK/DC initial pose and cache replay. This author has run
no tests, MATLAB or numerical jobs. Root owns checks, preflight, execution,
acceptance and publication; no new physical result is claimed here.

Root's first pure-Python readiness job reports terminal0 with frozen sources
matching, empty Job/absent owned PIDs and a separate idle ledger. Preparation
took23.102s for384carriers (368active):128positive frequency nodes through
.984615Hz,272unique headings,34816dense native pairs. Preserved uniformDC is
.01216723m; summed certified origin envelope exceedance is.00942218m. Full-hull
height/velocity and native net-force gates were NOT run in that allocation.
Canonical input bytes remained identical. These are readiness observations,
not a completed dynamics/surface validation. Subsequent generic source fixes
handle a Nyquist singleton without doubling and accept genuinely empty calm
Sea arrays in the native spectral packing; both await root's focused tests.

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
Boat/sea are unchanged. Adoption checkpoint7e1b854 is pushed. **Cleanup actions are
complete**. Independent current-state audit confirms source/output identity,
remote/local retirement and recovery preservation.
Only fork branches `main`0753b2e and `topwater/reference`6b remain; exact995 is
preserved by the verified tag `archive/topwater-legacy-9952634`, not an active
legacy branch. [Completion/recovery record](fork-cleanup-plan-2026-10-27.md).

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

- [Cleanup execution](../../../artifacts/workflow/wecsim-cleanup-execution1/) verifies
  native archival8/10, deletion of obsolete issues1/2/11/13 and nine remote heads,
  fork-only Actions disable, plus removal of11 old worktrees/13 trial refs.
  Merged PR6/7/9/14 and closed issues3/4/5/12 remain. The shared-store bundle,
  shallow boundary, six dirty and24 ignored snapshots passed recovery checks.
  Production6b, stocka23, the clean common parent/linked validation and source953
  remain; original failed studies are preserved.

## Limits and preserved history

This is implementation/numerical compatibility acceptance, not physical accuracy.
Issues31-34 remain unresolved. The native coarse-direction warning is a coverage
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

[Independent completion audit](../../../artifacts/workflow/reference-adoption-completion-review-v1/audit.json)
records the current source, API and recovery checks.

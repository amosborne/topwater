# Checkpoint handoff

## Current goal: expert-agent setup COMPLETE (independently verified)

All nine setup-checklist items have been verified; the development checkpoint
and issue30 record preserve the result. The downstream simulation remains pending.
The user cleared the simulation goal and explicitly requested the desktop expert
setup. Read docs/experts/setup-instructions.md, README.md, setup-checklist.md,
coordination.md and assignments.md. Stable rules/routing are in AGENTS.md.
Do not automatically resume retired coefficient studies or the interrupted
simulation preparation just because setup is complete.

Four project-scoped definitions are installed under .codex/agents:
capytaine_expert, wecsim_expert, topwater_expert, numerical_methods_researcher.
Each has name/description/developer_instructions only, no model/effort overrides.
Durable owner KBs are docs/experts/<role>/README.md and current-state.md, with
versioned detailed source/reproduction notes linked there.

Native activation is VERIFIED, not inferred from files:
- Codex0.164.0 strict startup-config validation exits0; model-visible native
  spawn agent_type selector lists all4 custom roles.
- Native CLI activation session31290 is TERMINALexit0. Exact native spawn calls
  and canonical returns are extracted in docs/experts/native-spawn-receipts.json;
  actual child session metadata binds each role/path/parent/UUID. Their selected
  model/effort match the native parent, independently checked without overrides.
- All4 returned developer-role/KB/scope confirmations and updated their own KB.
  See docs/experts/native-activation.json and client-verification.json.
- The pre-existing host chat's spawn wrapper lacks agent_type. Do NOT claim it
  hot-reloaded. The installed native CLI's named-role mechanism was exercised;
  reuse instructions are in docs/experts/activation.md. These are reusable roles,
  not guaranteed permanent live background services or permanent model learning.

Owned task status:
- CAP-001: Capytaine producer/interface readiness inspected READ-ONLY; minimum
  native export and joint interface execution queued, no numerical lease.
- WEC-001: WEC-Sim expert owns BEMIO. Fork created and pinned; source map and
  independent four-lead reproduction plans delivered. Actual reproductions/
  corrections remain ASSIGNED/QUEUED, not claimed done. No blind workaround port.
- TOP-001: held-state governing-equation/code/evidence map completed from docs;
  NO Topwater source changes/tests/simulations/optimization permitted.
- NUM-001: reusable proposal template and focused incoming-wave/IRF-padding
  question prepared; no broad literature survey or numerical experiment.
- Independent setup audit PASSED. Non-author verified raw native spawn calls,
  all child role/path/parent/model metadata, role file hashes and clean fork pin.
  No operational setup blocker; setup-checklist records acceptance.

Maintained fork: https://github.com/amosborne/WEC-Sim; local checkout
artifacts/dependencies/topwater-wec-sim; origin user fork, upstream official.
Branch topwater/reference-v7.1.0 is PUBLISHED at exact release/base/HEAD
 a23d2993807410fe78cc63546c1ba9d81912024f (v7.1.0), clean. No engine corrections,
new engine commits, numerical tests, upstream PR/messages or main push. Native
activation note records its earlier local-only state; docs/experts/wecsim-fork.json
records coordinator's subsequent unmodified branch publication and remote check.

Marine-gyro engine HEAD independently verified READ-ONLY as
0753b2e47f2457c078751dcfe5d251d1767b80ab, matching data/case_study.json, clean main
and official origin. Yaw order, unequal grids, angular-bin normalization and
row/column directions remain UNREPRODUCED leads on selectedv7.1.0; expected-math
plans and exact source maps are owned by WEC-Sim expert. See its current-state,
native-setup note/read manifest and marine-gyro-leads.md. No copied fixes.

Preserved work and processes:
- Uncommitted reference/capytaine/export_case.py remains EXACTLY unchanged,
  SHAa1feb6e8581695c1a7562bf1e129bbd609f96b54ab7dcac25da9b5ecdaa92942.
  It is incomplete/unrun (incident_packet/export_controls not frozen). Do not
  commit it as accepted export code, discard it, or launch it without owner review.
- No numerical job was run during setup; execution-state.json is IDLE/unleased.
  MATLAB shared engine CAPYTAINE_REFERENCE checkedUPidle; startup helper8651
  TERMINALexit0. Engine presence is not a numerical execution lease.
- Earlier quadrature36762 was user-interrupted TERMINALexit1, not a numerical
  failure or a task to resume. Local partial data retained unchanged.
- All old studies/failures are historical; no Topwater source or existing test
  was changed in setup. Current development branch baseline/fluid-body.

Downstream purpose remains minimum-mesh Capytaine -> BEMIO -> WEC-Sim10s freely
floating six-DOF reference, isometric video, then six-motion variation over mesh.
Boat starts at rest in locally calm water; waves roll in from upstream, no
imposed temporal ramp/global wave/load step. Force input and displayed incident
water must share physical arrival/direction/phase. Intermediate coefficients are
DIAGNOSTICS, not standalone accuracy gates. Setup proves no boat motion/video.

**Next bounded coordination action (only when setup is complete and reference
work is authorized):** use the reference experts, not main-thread reinvestigation.
Assign one execution owner for the smallest relevant WEC-001 reproductions and
case-path/interface checks while Capytaine producer planning remains read-only.
Freeze numerical expectations/resources/producer/input before any numerical job,
run serially with exact terminal handles. Independently review consequential
confirmed fixes in maintained fork. In parallel only read-only producer/consumer
contract work; then first minimum-mesh end-to-end10s result, not a comprehensive
whole-library audit. No Topwater work. Do not replace the user-visible reference
milestone with indefinite training, coefficient checks or optimization.

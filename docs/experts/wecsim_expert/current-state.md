# Current state: WEC-Sim / BEMIO expert

2026-10-24 native `wecsim_expert` activated as `/root/wecsim_setup`; actual injected developer role recognized (BEMIO/fork, startup/memory, six-DOF, output/video ownership). Durable context is files. Setup-only checkpoint complete; no numerical execution lease or numerical work, no engine fixes.

Maintained fork https://github.com/amosborne/WEC-Sim established after authenticated owner/existing-fork inspection. Local `artifacts/dependencies/topwater-wec-sim`, clean branch `topwater/reference-v7.1.0`, deliberately selected official v7.1.0 base/current `a23d2993807410fe78cc63546c1ba9d81912024f`. Origin fork; upstream https://github.com/WEC-Sim/WEC-Sim.git. Local branch not published; no content pushed. Existing upstream and marine-gyro checkouts preserved. This task has not numerically tested the pin.

Marine-gyro data/case_study.json recorded engine commit independently verified equal to its actual clean WEC-Sim checkout HEAD `0753b2e47f2457c078751dcfe5d251d1767b80ab`, main, official origin. Four reported issues remain hypotheses on selected release.

Read [native setup/source map/bug ledger/reproduction plan](native-setup-2026-10-24.md), [exact read manifest](native-setup-read-manifest.json), [historical leads](marine-gyro-leads.md), shared contract and coordination. Existing MATLAB R2025b/Simulink/Simscape/Multibody setup evidence remains docs/capytaine-wecsim-environment.json; idle CAPYTAINE_REFERENCE is not a lease.

Explicit subsequent assignment: independently reproduce yaw-order, unequal-grid postprocessing, angular-bin normalization and row/column directions against independent mathematical expectations on selected pin. Request sole numerical lease before executing. If independently confirmed, make minimal maintained-fork correction and focused regression; obtain separate independent review before physical reliance. No whole-library audit/global refactor, Topwater work or upstream publication.

BEMIO third argument supports separate hydrostatics directory (suffix with leading separator); sent actual selected importer contract to Capytaine child. Joint executed interface signoff pending. Unfinished/unrun producer missing packet/export controls remains unmodified; per-direction arrival/padding/full inertia/Ainf mapping remain first-case integration questions. Downstream milestone remains minimum-mesh 10s incoming-wave six-motion video, then motion versus mesh; this task stops before simulation.

Coordinator checkpoint update: the unchanged selected-release project branch is now
published to origin/topwater/reference-v7.1.0; see ../wecsim-fork.json. The native
setup note truthfully records its earlier local-only state. No engine fixes or
numerical tests have occurred. Assigned follow-up WEC-001 is in ../assignments.md.

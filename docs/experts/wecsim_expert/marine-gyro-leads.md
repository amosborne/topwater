# Marine-gyro leads: hypotheses until independently reproduced

Read-only comparison project C:/Users/alex/marine-gyro. Its data/case_study.json
records engine_commit0753b2e47f2457c078751dcfe5d251d1767b80ab. The expert independently verified the actual engine checkout HEAD is the SAME
0753b2e47f2457c078751dcfe5d251d1767b80ab, clean main with official upstream origin.
See native-setup-2026-10-24.md and native-setup-read-manifest.json for exact reads.
A default sandbox probe initially lacked permission; approved read-only queries
succeeded. This verifies version/configuration, NOT any of the four diagnoses.

1. reports/wec-sim/run_case.m around25 replaces kinYawTransformNoB2B with
   marineKinYaw by generated Stateflow script replacement. Inspect marineKinYaw.m
   and native quaternion sensor conventions. Derive expected combined yaw/roll/
   pitch and velocity/acceleration transforms independently. Do not copy patch.
2. run_case.m around67 calls waves.calculateElevation on body1_out.time before
   stopWecSim. Reproduce unequal output/wave time grids and determine engine
   versus case ownership before adopting a workaround.
3. docs/branch-review.md around127 reports missing angular-bin-width in full
   directional excitation; generated model substitutes
   Nwave.amplitude .* Nwave.dTheta.' for Nwave.amplitude. Verify selected native
   version against independently integrated spectral energy/force expectation.
4. Same report claims row/column direction samples produce erroneous broadcasting.
   Establish an explicit shape contract and independent expected force, then
   reproduce each orientation. A function matching itself is not independent proof.

Use the desktop expert-agent instructions' full assignment as authoritative.
For each lead record selected version, actual source path, minimal reproduction,
expected/actual behavior, physical consequence, cause confidence, correction if
needed, regression evidence and upstream status. A cleared symptom alone does
not establish cause. Keep unrelated speculative audits out of the first-case
critical path. No numerical tests before execution-owner lease.

# One genuinely milder observed buoy input: source readiness

Latest user steering supersedes energy-bug and incremental-official studies:
select a less energetic credible buoy sea, start locally calm with boat at rest,
incoming waves without a ramp/history seed, run full60s and produce the centered
video. Wider physical wave field and camera framing are separate. Topwater solver
implementation, tests, optimization, BEM and dynamics remain outside this task.

## Selected observation, not an amplitude-scaled old state

Coordinator selected `ndbc-46086-20230319T0840Z`, station46086,2023-03-19T08:40Z.
It is source California row21 in the frozen marine-sea-library/1 payload at
`C:/Users/alex/marine-gyro/artifacts/sea-direction-check/sea_states.json`, SHA256
`62162b862ed338e4e1dcac9e1476220b8d0550c0cbb42f981bc5c404129b1e56`.

Read-only metadata inspection finds Hm0=0.5351635264103863m,
Tm02=7.4269289105858185s, peak-density-band period11.4285714286s, all48 measured
usable bands .030-.495Hz, complete paired four int16 direction moments. It is the
only observation satisfying prospective Hm0.45-.55m/Tm026-8s mild-input criteria;
there was no response/RAO/phase fit or numerical wave realization in selection.
Observed energy-weighted propagation is23.3193deg CCW from east, first resultant
0.46278. This is broad directional spread, not a single directional wave. The
actual seed42 realization determines the default bow toward its arriving mean,
not this observation-only direction number. No rotation of only the rendered box.

The existing nine-record package's milder records are still Hm0~0.96m; the selected
record therefore comes from the complete source library using the public
`ObservedSeaState` constructor and unchanged `realize` API. The source library
stores lossless float64 band amplitudes and int16 moments. The packaged catalog's
five2023 NDBC source-file hashes/URLs accompany the selected record. The main
catalog and all old realizations are unchanged.

## Source and outputs

- Controls `reference/sea_states/calm-buoy-export.json`, SHA
  `8b9a69beabc3f9482913477fad93f2012a18d7026c44a31d46c7199f9673c8f5`.
- Producer `reference/sea_states/export_calm_buoy.py`, SHA
  `22bcf2c0e95f6f021eb7ffb9c1d3b5b34ffd0e1b954bbcd4e1d21ddde38e9621`.
- Archive preparer `reference/sea_states/prepare_calm_buoy_export.py`, SHA
  `53c8b17bf6b576a6f3e2dadd5f16e49fb6d19a29724dd37625d38c0951a33d78`.
- Ready proposal
  `artifacts/reference/calm-buoy-20261024/export-attempt1-source/proposal.json`, SHA
  `072f1fb5d03be38dfd63b621d5c3997e915057c5a8d6a2f47c20d4acdf22e3b9`.

Preparer only hashes/copies source and dependencies; no waves or numerical case
have been generated. Exact library bytes, physical input, catalog provenance,
actual generator/input modules, guard and numpy/scipy implementation files are
frozen. Python syntax compilation and git whitespace checks pass. Source requires
independent review, coordinator allocation and the serial claim/separate verify
procedure before running:

```
.venv/Scripts/python.exe reference/sea_states/export_calm_buoy.py --config reference/sea_states/calm-buoy-export.json --output artifacts/reference/calm-buoy-20261024/export-attempt1
```

OneCPU/noMATLAB, expected1-10s, internal30s/external60s. No automatic retry or
overwrite. Numerical generation retains seed42,8components/band and all measured
band energies. It checks exact JSON replay, full measured/realized band-variance
closure, physical frequency support and global/body phase-coordinate closure.
Outputs: observation.json, realization-global.json, components.npz with original
band/component IDs and both global/body arrays, case.json with complete R/t,
summary.json with actual realized direction, variance and output hashes.

## Physical consequence and stopping point

These outputs are stationary source carriers, not a calm-start fluid field.
WEC and researcher own the separately reviewed physical upstream finite-train
construction. Preserve carrier amplitudes, random phases and MEM directions;
declare any spectral broadening/window redistribution and distinguish the finite
train's local statistics from stationary observed Hm0. No artificial global
force ramp, stationary held release, phase fit or claim of unchanged finite-train
band energies merely because carrier amplitudes match. Same actual derived field
must drive full-hull initial eta/velocity/force checks and wider-region video.

Stop at reviewed neutral-export readiness pending coordinator execution grant.
No new sea realization/export/field/body solve or video has run in this task.

## Executed neutral export receipt

Coordinator subsequently allocated exactly this independently source-reviewed
proposal072f1fb5, unchanged source/settings. ONE export completed1.110s, all
energy/replay/envelope/frame checks PASS. Actual source/realized variance is
0.017900000000000017m2; Hm0.5351635264m/Tm027.42692891s unchanged. All384carrier
IDs0..383 retain band0..47 membership,8perband. Actual realized mean propagation
is21.3081431853deg CCW from east, first resultant0.4540120716; default bow is
201.3081431853deg, yaw3.51348990965954rad. Full R/t and corresponding body phase
mapping are stored in case.json/components.npz; frame elevation error8.33e-17m.

Output `artifacts/reference/calm-buoy-20261024/export-attempt1`:

- realization-global.json SHA
  `f08dfee04dc16d54f331cc4a5f991b4fab2987e6a4f8abe4f212ebadf51f7b21`.
- components.npz SHA
  `814644c993a957528432a6b5aed92d16d1ee94d52119e65e2e7dc26916b7eaa3`.
- case.json SHA
  `b63d82d280b5fd3a70ef5dab93e707e300d1ca7cdbf156efcf66931fc15b26ee`.
- observation.json SHA
  `7d2b13dcd51a4978c8a4bdd39d63d25dcba792d56511eae348793443726bd1c6`.

No wave-field/body/video simulation was run. Source carriers remain stationary;
their physical incoming-train construction/derived statistics are WEC/researcher
responsibility and require separate acceptance/allocation.

Authoritative prior-client absenceec7f87 preceded lease claim2fbec0; SEPARATE
verify2b3207 passed. Exactproducer2568/guard29816 terminal0/chunk50fd75; elevated
CIMe42561 confirmed both/children absent. All11inputs/dependencies/output hashes
verified and idle recorded UNDERlease90461a, exact terminal proof a7d34c,
SEPARATE idle readback7e3600, FINAL release84725c. No ledger writes after release.
WEC and coordinator received actual default yaw/carrier artifacts. Stop reached;
no automatic further generator/field/BEM/dynamics job.

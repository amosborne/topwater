# One stronger actual buoy input: source readiness

Latest user, after the complete mild60s video: use a more energetic sea and
report diagnostic checks such as radiation energy loss. Coordinator retained
the approximate5-10deg roll objective, without a response guarantee or amplitude
fit. Existing physical incoming train, locally calm boatrest/zero history/no ramp,
full60s and video remain required. This assignment owns neutral generator data
only, not Topwater solver, field, BEM or dynamics.

## Selection from observations

Selected `ndbc-46086-20231223T0040Z`,2023-12-23T00:40Z, California library row29984
in the same frozen source JSON, SHA
`62162b862ed338e4e1dcac9e1476220b8d0550c0cbb42f981bc5c404129b1e56`.
Actual Hm0=1.6574679484080534m, Tm02=7.23656714864799s, variance
0.17170000000000013m2, peak-density-band period11.4285714286s. All48 measured
usable bands .030-.495Hz and four paired int16 directional moments are complete.
Observed first resultant0.4576612545 and mean propagation3.5289535435deg; actual
realized heading/default bow must be determined from the unchanged seed42 export.

Selection used metadata only: prospective Hm0~1.7m/Tm02~7.42693s/resultant~0.454
within Hm01.5-1.8m/Tm026-8s/resultant.35-.6. No RAO/response computation, phases,
seed fitting, amplitude scaling or motion prediction. The record's variance is
about9.6times the mild observation's; this does not predict boat motion.

## Minimal reuse and exact readiness

New control `reference/sea_states/stronger-buoy-export.json`, SHA
`75755fbc66ec7f4ea7f0d88b045bbf5c6a09f9f896af100fa9adeb0498001f90`.
EXACT unchanged numerical exporter `reference/sea_states/export_calm_buoy.py`, SHA
`22bcf2c0e95f6f021eb7ffb9c1d3b5b34ffd0e1b954bbcd4e1d21ddde38e9621`.
The selected record is in2023, so its existing five2023 raw-file provenance remains
correct. Existing public ObservedSeaState/realize API, seed42/8perband, measured
energies, synthetic phase draws and MEM directions are unchanged. Default bow
faces actual realized arriving mean, not an independently rotated display.

New nonnumerical `prepare_buoy_export.py` accepts explicit fresh config/output
paths, prevents collisions and enforces the unchanged exporter/2023 provenance.
Preparer SHA `9dab68b8a9d65843b31ff3ecfe536a902419c980a020d6126dba22d482f303da`.
It hashed/copied11inputs plus numpy/scipy implementation identities; no generated
waves or solver invocation. Python syntax and whitespace checks pass.

Proposal `artifacts/reference/stronger-buoy-20261024/export-attempt1-source/proposal.json`,
SHA `43faa35a603eed152c289c9475430c1db20841056460e81dc91c62985ef3eed9`.
Independent review and coordinator serial allocation are required before:

```
.venv/Scripts/python.exe reference/sea_states/export_calm_buoy.py --config reference/sea_states/stronger-buoy-export.json --output artifacts/reference/stronger-buoy-20261024/export-attempt1
```

OneCPU/noMATLAB, expected1-10s, internal30s/external60s. Outputs/checks reuse the
accepted mild exporter: observation.json, realization-global.json, components.npz,
case.json, summary.json; all48variances/exact replay/frequency support/body-global
phase closure. Original mild input/train/body/video stay immutable.

## Acceptance limits and stopping point

Stationary source carriers are not themselves locally calm incoming waves.
WEC/researcher own fresh incoming-train assessment and derived finite-train
statistics, native force/body/video. The existing full-heading fixed17 excitation
data can be reused by those owners; no BEM or mesh study is implied here.
No amplitude/phase adjustment, no spectrum rescaling or promise of5-10deg response.
Stop source-ready pending review/grant. No export/field/body/video job has run.

## Executed export receipt

Coordinator subsequently allocated exactly the independently reviewed proposal
43faa35a. One generator export completed1.093s, unchanged producer22bcf2c0/config.
All48band variances/replay/envelope/global-body phase checks PASS. Realized
variance0.17170000000000013m2, source Hm0/Tm02 unchanged,384component IDs0..383.
Actual realized mean propagation2.0671942254deg, first resultant0.4711207531;
default incoming-facing bow182.0671942254deg, yaw3.1776719991015314rad. Actual
carrier frequencies .0302418612652-.494810289810Hz; source edges .030-.495Hz.
Stored carrier frequencies, global synthetic phases and original IDs are exactly
equal to the mild seed42 outputs; actual measured amplitudes and MEM directions
change. There is no uniform amplitude rescaling or phase selection.

Output `artifacts/reference/stronger-buoy-20261024/export-attempt1`:

- realization-global.json SHA
  `ef0df4cc506d1e34620f9bf8535f40868d64539aaac722bb9aca61a6607d7a2e`.
- components.npz SHA
  `4ddb214adb7de476372d0e8b4e5ee38ff044bb84c51f9040ab82046b2b77cbf6`.
- case.json SHA
  `ea2794b3b8bd11c45c2474a51dfab0be155126c3560754c2634a4f037b865a39`.
- observation.json SHA
  `d5ddc42080b10b3ef67f81fc37b2a3305a08d9273a7eee95ecb5434fe51caa22`.

Preflight found two transient read-only Python handles25952/33204; authoritative
same-handle recheck88c58c proved both/children absent before claim. Prior recovered
encoder17992/3976 absent. Claim864953 -> SEPARATE verifyb87fc8 -> runb67e1d,
exactproducer38880/guard32504 terminal0. ElevatedCIM1754c2 proves both/children
absent. All11frozen inputs/dependency/output hashes verified, idle UNDERlease450501,
SEPARATE idle readback1f807f, FINAL release0b0fc7. No ledger writes after release.

WEC/coordinator received actual default yaw and carrier artifacts. No physical
train, calm assessment, native motion or video was run here; each needs a fresh
separate accepted consumer binding/allocation. Existing full-heading excitation
coverage may be reused by its owners. Stop reached; no automatic numerical work.

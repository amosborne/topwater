# Minimum-mesh end-to-end reference: completion audit

2026-10-24 local date. Latest user scope: actual buoy-derived stronger sea,
calm incoming-wave startup without a temporal ramp, full60s, appreciable6DOF
with roughly5-10degroll, video, and radiation-work diagnostics. This closes the
minimum-mesh workflow and requested deliverables, not a precision physical
benchmark or nonlinear-model qualification.

| Requirement | Authoritative current evidence | Finding |
| --- | --- | --- |
|Native producer/official BEMIO/native WEC-Sim|Original minimum native-gl2/consumer H5; required fullheading excitation-attempt1; actual stronger body frozen147 inputs and raw output|Complete; same17-panel boat, original H5/CIC10/dt.02, native Simscape integration|
|Actual collected-buoy generator, no scaled old sea|stronger-buoy-20261024/export-attempt1/observation.json, components.npz and summary.json; independent rawsource row29984 audit|Hm01.657467948m/Tm027.236567149s, seed42/384carriers; measured a/MEM change, omega/phases/IDs unchanged seed draws|
|Bow toward dominant arrival|Actual export meanprop2.067194deg/bow182.067194deg; actual summary R/t and all751 transformed hull poses|Embedding applied once to field and boat; native physical axes unchanged|
|Locally calm/rest/zero radiation/no ramp|Actual train and assessment-attempt1 accepted report; native preflight/raw initial rows|Complete under declared eps.001 engineering policy: eta<=.83572mm,u<=3.686mm/s,Fnorm6.104N,Mnorm12.016Nm; bodyq/v/memory0, ramp0|
|No extra viscous/PTO/mooring damping or prescribed trajectory|Actual native force channels and frozen input/model checks|All6free, viscous0; radiation remains native. No trajectory prescription/history seed|
|Full60s, no prerun/cropping|Native raw3001x60, summary duration60/dt.02 and exact0:.02:60 time grid|Complete; all startup and late motions retained|
|Strong roll and appreciable other5DOF|Actual full3001 native q|Peaks xyz[.562349,.836958,.758899]m,rpy[6.618565,2.677850,1.325248]deg|
|Faithful incident forcing and native consumption|All3001x6 prepared/applied excitation, actual CIC three-row/PCHIP, Solid Newton-Euler and Ar checks|Fe errorEXACT0; CIC<=9.095e-12,NE<=1.4552e-11,Ar<=3.34694e-10; independent PASS|
|Radiation-work check|Native all6world-CGvelocity/worldomega dot resistingFrad, saved-grid integral; independently reduced|Net32971.03955J,mincum0,returned4074.89555J; incident46582.19287J. Not complete fluid energy closure|
|Close centered full-duration video, correct scale|Exact accepted cache+ROOT/independent literal previews and decoded0/30/60 frames; full751poses reconstructed independently|12m orthographic/equalmetres/no amplification; no clipping; field is prescribed incident input, independent of camera crop|
|Actual deliverable decodes fully|render-encode-attempt1/final-review/verification.json; ffprobe full-frame count; actual MP4 SHA42c8ca1c...|H2641280x900/751frames/12.5fps/60.08s, native60 endpoint held one frame; all hashes and literal checks PASS|
|Serial execution and terminal recovery|Shared execution ledger and exact process/idle/readback/release receipts|All physical/render jobs terminal, no lease. Stronger encoder exit0/FFmpeg0; earlier mild owner-loss exit status honestlyNULL with complete output verified|
|A-priori linear-use assumption|Existing exact geometry-envelope curves and stronger-relative-water-proxy note|Approximate fixed-mean linear use is stated, not certified: proxyrelativeheel3.10deg/trim2.88deg/relheave.090m; existing3deg staticrighting~7.3%/tangent~15.3% deviations. Curvature/wetting/combined geometry accuracy unqualified|

Main artifact root: `artifacts/reference/stronger-buoy-20261024`.
Final video: `render-encode-attempt1/reference-isometric.mp4`, SHA256
`42c8ca1c80e961784200d076faff7e9a34a67c782120019e606bbf553ec004e2`.
Motion SHA256 `3eae21f9efa9987d77225b0ebe3f2cac4babec0a482ea01bbbf7801370e6b072`.

Reproducible stage commands/controls are in
[stronger-buoy-pipeline.json](../../reference/wecsim/stronger-buoy-pipeline.json).
Source Python/MATLAB checks and actual stage outputs pass; independent reviewers
were distinct from consequential source authors. No Topwater solver work or CI
was run. Active branch baseline/fluid-body only; final commit/push/issue30 and
handoff record the integration checkpoint.

Physical limitations are part of this result: coarse spatial/interpolation,
finite radiation memory and time-step accuracy have not been converged. A small
geometry-error or nonlinear dynamic validation claim is not established. User
asked for an a-priori assumption rather than proof from an unavailable nonlinear
model; static/proxy context is reported transparently. Tighter applicability or
accuracy requirements need a separately scoped follow-up. No further studies
are automatically authorized by completion of this workflow.

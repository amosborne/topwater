# Checkpoint handoff

Updated 2026-10-27 UTC. Latest user steering controls scope. Historical evidence
and queued assignments do not authorize execution. Follow AGENTS.md/coordination.

## Current priority: packet timing and 1D demonstration

Latest user steering supersedes the former next-speed priority. Step back before
further radiation reconstruction or optimization to investigate wave-packet
timing and simulation duration. Retain quintic smooth-join constraints. The CLI
requested time means DWELL (full-amplitude interval); investigate ramps tied to
each original carrier period, initially2cycles per ramp, with2–3s calm startup
and measured settling afterward. Thresholds/alignment/runtime changes remain
under investigation; no production behavior has changed.

The requested standalone1D demonstration is complete:
[time/frequency construction](wave-packet-timing.md) and
docs/images/wave-packet-dwell-time-frequency.png (2430x1530,
SHA48420d5c703abdff089c26deedb946bda4a8b74ab4222e86be3a61756893f093).
Two original periods4s/6s have ramps8s/12s and identicalDWELL14.5–24.5s.
Four panels show envelopes, time traces, envelope spectra and packet spectra.
Complex FFT/IFFT reconstruction error <=4.811e-16m; Parseval error4.441e-14m2s.
This is an intended1D profile demonstration, NOT propagated wave physics,
startup/calm certification or a WEC-Sim run. No previews/images were inspected.

Root job wave-packet-timing-plots-attempt1 terminal0/frozen/guard0/driver0/
Jobempty, exact producer32352/owner9992 and all guarded host identities absent.
finish_check terminal0 then idle/null readback/release-last completed. Current
ledger idle/null, lease absent. NMR authored illustrative source; TOP reviewed
source; root alone executed. Actual plot work0.844s. Report/complex spectra are
under that attempt's results/; tracked numerical evidence is
docs/experts/numerical_methods_researcher/wave-packet-timing-demo-2026-10-27.json.
All canonical output and production sources remain unchanged.

[Issue36](https://github.com/amosborne/topwater/issues/36) is the new focused
timing investigation, recording DWELL and cycle-based ramps. Remaining work is
to establish physical packet alignment/propagation, threshold-defined startup
and ending, sustained measured settling/cap behavior, input/CLI implementation
and appropriate new-case qualification. No further numerical job is allocated.

## Completed goal and acceptance

Issue #31 is resolved and closed (2026-10-27T18:45:47Z). Accepted checkpoint
`2fe5e96` is pushed to origin/baseline/fluid-body, including both canonical LFS
objects. Issue30/31 authoritative bodies are updated and remotely verified;
issues30,32,33,34 remain OPEN. No upstream or WEC-fork change was made.

Requirement: maximum numerical water-elevation difference PLUS reference
uncertainty <=0.010 m against the chosen LINEARIZED model over the full declared
5929 targets × 751 frames / 60 s, including continuation. This does not measure
error against nonlinear physics or experiments. Preserve accepted WEC-Sim
motion/forcing, boat/sea, mean hull, incident/diffraction/radiation, frame/clock
identities, and all physical opening/contact/startup/export gates.

Full-region qualification and independent artifact review PASS. Maximum actual
runtime/reference difference 0.341183 mm; maximum pointwise difference plus
reference uncertainty 5.784425 mm. No motion or force recalculation occurred.
The production readiness flag is now TRUE; its promotion changes only the flag
and comment. All production function/class ASTs match the qualified source.
The gate regression explicitly exercises both False rejection and True admission.

## Execution state and final evidence

Final integration regression job issue31-final-tests-attempt1 PASSED50 tests in
24.20s with no skips. Exact root session82774 terminal0/frozen, guard0/driver0/
Jobempty and all owned host identities absent. finish47428 terminal0 verified
idle/null readback then released the lease last. Independent receipt review
issue31-final-tests-review1/acceptance.json SHA
b394c85da947f3c4549e02f7c9fc88da6556102a5d7617bef70b699a96f61396 accepts
actual coverage and ENTIRE module AST unchanged apart from readiness Boolean.
No jobs remain allocated. Final saved-array metadata reporting also passed:
issue31-max-location-attempt1,1.188s, terminal0/frozen/guard0/driver0/Jobempty,
owned host identities absent, finish59599 terminal0/idle/readback/released.
Both maxima are at44.48s,target3283; actualdifference0.341183mm plus pointwise
uncertainty5.443242mm gives5.784425mm. No qualification or data changed.

Documentation and warning triage are complete. The tracked evidence bundle is
docs/experts/numerical_methods_researcher/issue31-acceptance-evidence-2026-10-27.json,
SHA10a1e9489eb9d32104a784758d6ec6e819b288c49b3b9ec492a51bb2383919a1.
Source, artifacts, tests, review, documentation, warnings, LFS publication and
issue closure are complete. Root's requirement-by-requirement audit found no
remaining #31 work. No jobs remain live/allocated; ledger idle/null, lease absent.
The packet-timing investigation above now precedes further reconstruction/speed
work. Future jobs need their own bounded serial allocation.

## Accepted method and scope

Runtime preserves intact complex H_eta,v with FFT(v), all six saved native
velocities and their exact knots, zero prehistory, linear interpolation,
1024 s period / half-native timestep (.01 s), 20 s C1 Hermite continuation to
rest, 48 rad/s cutoff / .025 spacing / 128-point chunks. Raw startup remains;
no extra acceleration term, phase fit, causal-only truncation or forced t0 zero.

Independent reference separates the actual finite-low discrete operator's
acceleration endpoint A, applies regular R=H+i*omega*A to C1-Hermite(v,a),
and adds exact saved A*a for validation only. Pointwise uncertainty covers
padding, timestep, interpolation, frequency spacing/cutoff, future/prefix
controls, startup and the omitted regular continuum. The latter assumes
|R| <= twice the same-location maximum calibrated from dense36–48 and64/96,
with known32/prior-probe checks. This is an explicit selected-case engineering
envelope, not a theorem excluding unseen resonances or a model-resolution bound.
Exact-input Fourier enclosure uses8192 s quadrature period through512 rad/s
plus derivative-jump remainder to infinity. Original failed policies/studies stay
preserved. See the immutable Hermite/tail derivations and final acceptance note
under docs/experts/numerical_methods_researcher/.

## Evidence and published artifacts

- Final independent acceptance: artifacts/workflow/issue31-final-method-review1/acceptance.json,
  SHA778e72871f57f90b01ba7940809cfb2561e1dd2350aaa5c16f4ec4552933c318.
- Independent every-array checker: issue31-independent-check-attempt1/results/report.json,
  SHA55575ea13936296500202d9234e75e88f45c5c6e03fe74a94ab9a618ffe3eb50.
  Every packed water vertex/connectivity matches the qualified field exactly;
  common incident/diffraction/TOTAL identities agree to about1e-15 m.
- All751 geometry gates: issue31-recovery-attempt1/results/report.json,
  SHAb1e8f511d5d3833ab816bb3ac054d47c2f4f84b61da13cf42c2c6fd5d3071ae3.
  Replay540.829s, including wet-set355.109s and disturbance metrics157.192s.
- Published examples/simple-boat/result.npz SHAd3438810cfd0aa899600c90d3b885b8e02abf05c354741aacaf8cea113f387c5;
  result.png SHA275ac4c519ff0c2ef0ca557b798e7d3a15cecb3123c1756f19a4351f3485dfe8
  is regenerated and byte-identical because all3001 motion samples are unchanged.
  Publication report: issue31-publication-attempt1/results/report.json. Backups
  of all8 previous canonical files are retained under its results/prior-canonical/.
- Actual public render: issue31-render-attempt1/terminal.json, terminal0/frozen;
  guard0/driver0/Jobempty, owned process absence independently checked, finish
  completed idle/readback/release-last. video.log confirms751 decoded800square
  frames/12.5fps/60.08encoded seconds. Command142.641s, renderer128.688s.
  No previews or screenshots were generated or inspected.
- Readiness promotion proof: issue31-readiness-promotion1/report.json;
  source6d8edb3e6c21e21878715380a0fdb2ca2480f9c88e344dc6896c54e7e11ed575,
  gate-testb91ed02440277fda9e82ef405e8597e6be053f6bd6846a3f85e021121a36e357.
  Previous source bytes archived and all production definition ASTs unchanged.
- Focused full-complex regressions:42 PASS17.81s at d3630fe (unchanged functions),
  artifacts/workflow/issue31-offline-tests-attempt3/terminal.json.

Canonical sea/video logs have0warnings, result.log7:6 inherited calm-restoring
records and1 current wet-set warning. Wet-set departure peaks18.27% at44.64s,
exceeds10% from29.04–53.12s. This is fixed-mean applicability, not nonlinear
motion-error percentage. Actual reconstructed opening/contact gates pass.
Total-minus-incident fullpatch RMS5.08038mm/1.90312%; onebeam5.99504mm/2.22013%;
resolved peak39.1241mm at43.36s. These are disturbance metrics, not error.

## Priorities and ownership

31 is complete for its original case. Current priority is packet timing/start-stop
design and the requested1D time/frequency plots, BEFORE further reconstruction
or measured-speed work;32(model resolution),33(field/load consistency),34(linear
applicability) remain separate later work. No speed/model campaign is allocated.
CAP owns physical field implementation, NMR independent method review, TOP
execution receipts/integration helpers; root owns serial execution and acceptance.
Native Topwater dynamics remain a nonworking prototype. WEC reference6b774228
on officialv7.1.0 with ordinary public coefficient input remains unchanged;
fork cleanup is complete and upstream contributions deferred.

Preserve all old failed studies and other worktrees. Historical evidence and
WEC adoption records are linked from the expert notes and
[archived handoff](experts/history/checkpoint-handoff-2026-10-27-before-full-qualification.md).
Commit/push substantive checkpoints only to baseline/fluid-body; update issue30
main text and handoff. No main push, assistant merge, bookkeeping PR, upstream
issue/PR or outside-message authorization. User forbids render preview self-review.

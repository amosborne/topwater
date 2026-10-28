# Checkpoint handoff

Updated 2026-10-27 local date. Latest user steering controls scope; historical
evidence and queued assignments do not authorize execution. Follow AGENTS.md.

## Current outcome: optimizer-free quintic examples, issue #36

The user's latest requested example is COMPLETE: three base frequencies
0.01/0.1/1Hz, phases0/45/90degrees DIRECTLY at dwell MIDPOINT65s, each shown
for N10/20/40. Exact profile:10s rest,5s ramp,100s dwell,5s ramp,10s rest.
Latest user preference is NO OPTIMIZER: fixed quintic envelope, direct Fourier
coefficients, largest-N magnitude selection. No LP/greedy method in this run.

Plots: docs/images/wave-quintic-n-time.png (2720x1760),
wave-quintic-n-spectrum.png and wave-quintic-n-phase.png (2720x1600).
Rows are base frequencies; columns are N10/20/40. Images were NOT inspected;
source/metrics/header/dimensions/hash verified for user review. Narrative and
reproducible exact source: docs/research/pulse40-three-examples.md and
docs/research/quintic_packet_example.py. Full coefficients/metrics/identities:
docs/experts/numerical_methods_researcher/quintic-packets-2026-10-27.json.

Actual helper time6.8867303s, guard8.5s. Each case normalizes its midpoint sum
to unit amplitude, applies requested phase directly, uses j/130Hz positive
harmonics and keeps nested10/20/40 largest magnitudes. Symmetric envelope
magnitudes hold to2.2e-14; anchors to3.4e-16; quadrature128/256 differences
<=3.1e-16. Tailbound separates omitted positive coefficients from top40.

IMPORTANT limitation: this is positive projection of the COMPLEX template
a(t)*exp(i(2pi*f0*(t-65)+phi)), NOT ordinary real-wave Fourier truncation.
DC/negative analytic content is discarded, never folded; normalization changes
the rawL2 projection. N40 worstcomplexscore bounds29.1142/7.6348/5.8205percent;
real-target maxerrors24.4837/9.9826/5.9692cm forunit1m. The0.01Hz case loses
2.895% analyticenergy toDC/negative content and has+3.8236cm meanbias. This is
a method-specific limitation, not a physical impossibility or accuracy pass.
Do not call these minimax optima, a production sea or qualified boat physics.

Raw evidence: artifacts/workflow/pulse40-three-quintic-attempt1/results.
Sourcea957e80224959bec5c3bc1b694a498568188bc9d7348b9b1bb45614e6a69244c;
source/actual independent reviews accepted consistency, not quality. Root alone
ran the guarded job; all9cases complete. Producer32800/owner32496/guard33420/
driver20416/conhost7680/python37344 hostabsent, terminal0/frozen/guard0/driver0/
Jobempty, finishterminal0; separateidle/null/nojob/nolease readback confirmed.
No next numerical job is allocated. No production/CLI/canonical/WEC-fork change.

User now has a fully explicit per-component interface:
f0,phi(midpoint),rest,ramp,dwell,N. Bandwidth/frequency placement/quality outputs.
Any next study must follow latest user steering; do not resume optimization,
surface reconstruction, dynamics or performance campaigns automatically.
At checkpoint commit/push baseline/fluid-body and update authoritative #36/#30
main bodies; preserve all historical failed-policy studies under originalpolicy.

## Preserved superseded studies

pulse40-three-attempt1 and pulse40-three-midpoint-attempt1 were CLAIMED ONLY,
then cancelled beforelaunch when user changed the requested formulation. Their
cancelled-unlaunched.json records preserve verified manifests; no producer was
started. The earlier pulse40-three-symmetric-attempt1 ran3.010s with start-phase
anchor under its old policy and was fully cleaned; do not publish it ascurrent.
wave-grid100-attempt1 ran8.922s with an obsolete common-envelope/final-band
restriction (55s/12s ramps), not the user's current per-component definition.
All old source/config/results stay unchanged; no old publisher scripts should
be run against the current issue descriptions or canonical examples.

## Completed current evidence

- Latest illustrative plot: docs/images/wave-packet-dwell-time-frequency.png,
  2380x2040, SHAaacb6df90f8af35c42b987abfe918203969c5d6d0906c0ff634587c243225c2e.
  Periods4/6s, ramps8/12s, commonDWELL14.5–24.5s, specified phases both0 AT DWELL.
  Both carriers are1m crests there. Discrete cosine amplitude/phase bars; phase
  panels use tau=t-14.5 and x0. Fourier bins can have nonzero envelope phases.
  New waveforms/coefficients intentionally differ from former globalt0 examples;
  old inputs/envelopes remain preserved. Global/dwell inverse errors4.441e-16m,
  selected independent direct sum1.998e-15m. This is an illustrative time trace,
  NOT a propagated field, calm certification or WEC-Sim run.
- Phase plot source/run: wave-packet-dwell-phase1/ and attempt1/results;
  new complex spectrum SHAea1e42032205bc04747962d4f255784c94afd71790a40185b57062df5db4f1c8.
- Actualsea statistics: all384carriers,0.03024–0.49481Hz, periods2.021–33.067s,
  variance-weighted median9.795s, central90% variance4.308–15.971s. Two-cycle
  ramp median19.59s, maximum66.13s PER ramp. No carriers changed/dropped.
- Cheap COMMONseconds quintic baseline: referencex0 only,60sDWELL and existing
  cutoff0.53345874Hz. R1/2/4/8/16s maximum target-to-bandlimited differences
  98.074/13.615/6.967/2.363/0.391mm. Explicit DC removal;1024→2048padding
  R1/2s changes dwell by0.321/0.199mm. These are temporal projection screens,
  NOT shortest physical ramp acceptance or whole-hull1cm qualification.
  Existing requested native positiveomega grid has671nodes independent ofR;
  diagnostic FFT bins are a different grid, not a complete cost measurement.

Tracked reports and explanation:
[construction/timing note](wave-packet-timing.md),
docs/experts/numerical_methods_researcher/wave-packet-timing-demo-2026-10-27.json,
wave-ramp-period-stats-2026-10-27.json and wave-ramp-bandwidth-screen-2026-10-27.json.
Raw results/helpers remain under artifacts/workflow/ using the corresponding
wave-packet-dwell-phase*, wave-ramp-period-stats* and wave-ramp-bandwidth-screen*
names. Prior figures preserved in Git5aa0170/84fbc8f and private backups.

All THREE latest root jobs are terminal0/frozen/guard0/driver0/Jobempty with exact
guarded host process absence independently checked, finishterminal0/idle/null
readback then release-last: phase22888/34540;stats3440/16496;screen12460/3292.
Latest ledger idle/null, lease absent. No images/previews were inspected.
NMR authored bounded helpers; TOP reviewed source; root alone executed.
Separate /root/minimal_pulse_research authored only the research note; original
NMR reviewed it independently. No numerical job is currently allocated.

## Next scope and ownership

Issue36 is the authoritative investigation; update its main body and issue30
at checkpoints. Deliver the three options/tradeoffs to the user before selecting
and implementing a new synthesis policy. No automatic optimizer/BEM/WEC/video
campaign. Keep further reconstruction/speed work behind this design decision;
issues32–34 remain separate. Shared contract owns physical requirements.

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

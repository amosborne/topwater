# Checkpoint handoff

Updated 2026-10-27 local. Follow AGENTS.md; the latest user request controls scope.

## Current request: re-review surface reconstruction padding (#31)

The user accepted the appearance of the new simulation and asked to revisit #31,
specifically whether padding is still needed now that the packet has calm ends.
They also asked whether WEC-Sim used10ms: YES,13001 samples through130s. The
separate reconstruction FFT currently uses5ms; surface/video sampling is80ms.

Issue31 is reopened with a current-review narrative and preserved historical60s
acceptance. NMR read-only review separates the invented20s motion continuation,
1024s periodic FFT record and half-native-step interpolation. No setting or code
has changed and no new numerical job is allocated. Quiet waves alone do not prove
the boat settled. Proposed cheapest next check: saved endpoint v/a and a few
all-six-mode height probes, changing native10ms FFT, zero inventedfuture and a
shorter record one at a time. Probe agreement is not whole-region qualification;
the target remains1cm numerical linear-model elevation error including reference
uncertainty. Do not automatically launch a large study or adopt a padding policy.

## Completed goal: packet integration and simple-boat regeneration

Implementation and numerical work are complete and published at e0d6ba8 on
baseline/fluid-body, including LFS outputs and updated issues30/36. No numerical
job remains active; finish60701 and independent idle/null/no-lease readback passed.
Do not resume a physics, correctness or performance campaign without new scope.

The confirmed case is rest10/ramp5/dwell100/ramp5/rest10:130s total, fixed10ms
packet/native clock, phases directly at dwell midpoint65s. Canonical boat.npz,
sea.npz and sea.log are byte-identical. New result.npz/result.png/result.log and
video.mp4/video.log are published. The case has13001 motion samples and1626
surface/video frames at80ms. Video fully decoded:800x800,12.5fps,130.08s encoded;
PNG2160x1170. The example directory contains exactly the eight intended files.

Source integration checkpoint e416ca8 is pushed on baseline/fluid-body. The last
source delta changes only the CLI render budget to max(240s,water-frame-count),
its tests and documentation.170 integration tests passed56.93s;25 focused CLI
tests passed3.99s. Render API default remains240s. Physical solver/rendering
algorithms and the maintained WEC-Sim6b774228 reference were unchanged by that fix.

The actual130s native calculation passed all production phase/startup/fullmass/
inertia/Newton-Euler/applied-force/CIC gates. Its original solve invocation then
failed during root-owned MATLAB cache restoration because TMP was inside the
initial CacheFolder. Corrected sibling folders passed a real restoration control.
The saved native data were audited (phase error1.7583157152500917e-14m) and genuine
rebuilt operators matched all native numeric/netCDF values before fresh total
incident+diffraction+radiation field construction. Completed field hash is
9f3870fb9f3e91deb9b82ea83b25f62ed24e6acb43ad11dfd0d1acf3157660fb.

Existing geometry replay passed all1626 opening/contact gates and preserved
cached motion exactly. The final ordinary render succeeded after the original
240s video ceiling was corrected. Logs identify these continuations honestly;
do not claim one successful original solve invocation or a repeated native solve.
All failed attempts, complete fields and original eight-file rollback backup
remain under artifacts/workflow. No previews or image inspections were performed.

Final successful job: session82186, packet-render-finish-attempt1/manifest.json,
packet-render-finish-guard1/guard-receipt.json, helper packet-render-finish1.
Terminal0/frozen-match/guard0/driver0/Job-empty and all owned host process absence
are verified. Field replay proof is packet-field-replay1/results/report.json;
native/operator/control proof is packet-native-recovery2/results/report.json.
Tracked run evidence: docs/experts/packet-integration-run-2026-10-27.json.

Current result is a provisional linear reference. Old selected60s5.784mm surface
qualification does not qualify this changed packet. Mesh/irregular-frequency,
field/load closure and fixed-mean applicability remain separate #32–34 work.
Actual wet-set difference peaks25.13%; disturbance is6.836mm RMS/full patch and
8.036mm RMS/near band, with51.664mm peak. These are disturbances, not errors.
One nonblocking log label uses computational peakZ=.4378107785m; physical world
Z adds preservedDC=.0121672297m, giving.4499780082m. Published geometry is correct.
Timing/settling policy remains fixed for this case, without automatic extension.

Coordinator owns final acceptance/checkpoints and issues30/36. Experts have
completed their bounded source/doc work. Keep baseline/fluid-body; preserve other
worktrees. No upstream submission, merge or WEC-Sim-fork change is authorized.

## Historical checkpoints (not queued authorization)

## Current outcome: fixed10ms minimum-cutoff utility implemented

Latest authorized implementation COMPLETE: src/topwater/wave_packets.py exposes
quintic_packet/minimum_cutoff_packet. Fixed TIME_STEP_S=.01; durations align10ms;
odd/evenFFTcounts supported. Solepasscriterion ABSOLUTE max(0,abs(wave)-A*a)
<=.01m includingrest/ramps/dwell; prior dwellmatchrequirement superseded.
Eligiblecutoffsstartfirstbin>=carrier, all lowerbinsstillretained. FFTonce,
fullarrayincrementalharmonicsperbin, onefinalIFFT, no optimizer/IFFTinscan.
Between-samplecurvature+roundoffguard gives firstCERTIFIEDcutoff, notexact
continuousminimum. No sample_count/finergridAPI. WEC-Sim currently.02s unchanged;
user wants eventualworkflowaligned10ms, futureintegrationseparate.

21tests PASS0.10s in rootserialguard; completeAPI14.05/26.33/61.82ms for
f0.01/.1/1Hz, phi0/45/90atmid65,rest10R5D100,A1,tol.01m. Cutoffs.276923/
.353846/1.253846Hz,36/46/163positivebins. Certifiedbounds9.9173/8.1251/9.0022mm;
guardonlyrejection0/0/1. DwellerrorsdiagnosticONLY. Plotwave-minimum-cutoff-10ms.png,
docs/wave-packet-search.md and packet-cutoff-search-2026-10-27.json current.

Source/actual reviewed independently. Test+benchmark guardterminal0/frozen/
driver0/Jobempty, allownedprocesseshostabsent; finishsession48349terminal0;
separateidle/null/noleaseconfirmed. Noactive/nextnumericaljob. No CLI/seaarchive/
WEC/native/canonicalintegration in this change. Rootcoordinatorownscheckpoint.
Previousruntime/source/failedstudies staypreserved; userlateststeeringcontrols.

## Earlier examples and superseded implementation proposals

LATESTuserauthorizesIMPLEMENTminimumcutoffsearch. LoadedTopwaterexpert
receipt_reliability owns src/topwater/wave_packets.py, tests/test_wave_packets.py,
docs/wave-packet-search.md. IndependentNMRissue31_method_review owns read-only
review; minimal_pulse_research owns NEWpacket-cutoff-search-example1 driver/plan.
Rootonlyrunsreviewedtestsandbenchmarkseriallyafterclaim. Scopeonlynewpacket
utility, no productionsea/CLI/native/WEC/canonicalintegration or physicsjobs.

LATEST timepolicy USERexplicit FIXED10ms forconstruction/search/output, not
32768samplesdefault and nofiner-gridsecondIFFTbydefault. T130=>13000FFTpoints.
Each timingconstant mustalign.01s; oddtotalstepcountsneedproperrfftHermitian
handling. TOPupdatesmodule/tests/docs, NMRbenchmarkdriverupdatede01aa496...
plan7f24953a... no sample_count/finerargs. AwaitNEW10msreviewreceiptbeforelaunch.
Core1cmalltimeenvelopegate/incrementalalgorithmunchanged. Userquestions about
peakchecks/timestepareconceptual, notfurtherimplementationrequirements.

FINAL criterion: userSUPERSEDES prior dwell-accuracycondition. Only all-time
envelopeEXCEEDANCE max(0,abs(realwave)-A*quinticenvelope)<=0.01m absolute.
Rootstatedeligiblecutoffsstartatfirstbin>=f0(carriermustremain), not DC-only.
InputsphaseDIRECTmidpoint, amplitudeexplicit(default1),rest/ramp/dwell/tol_m.
FFTonce; accumulateDC+onepositivecosine/conjugatepairperstep onFULLgridplus
exacttimingboundaries; noIFFTinsidecandidate loop, onefinalIFFTvalidation.
Firstpassingeligiblebin; no binarysearchpeakmonotonicassumption. Dwellerror
diagnosticONLY, notgate. Reportminimumonvalidationgrid, no continuumproof.
Userapprovedthisdesignwith'ok sounds good do it'. Measureoverallruntime.
Oldfmaxplotsalreadyshown4930cba; no rerunpending. Lastjobterminal/leasefinished.

LATEST authorization replaces N input with explicit maximumPHYSICAL frequency.
Buildrealquinticwave time→FFT→retainDCand|f|<=fmax→IFFT; allbinsbelowcutoff,
no largest-coefficient ranking/optimizer/rescaling. User requested3cutoffplots.
Rootselected1.05/1.2/1.5Hz acrosssamef0.01/.1/1Hz,midpointphi0/45/90,
rest10/ramp5/dwell100/ramp5/rest10. Npositive becomesoutput. Newsourcepulse-fmax-
design1, reviewpulse-fmax-review1, guardpulse-fmax-guard1. Rootsole60s1GiBjob
pendingreview, noagentexec/noimageinspection/no productioncanonicalchanges.
Userobservesnextstepcouldsearchminimumcutoff foracceptablewaveformerror;
rootexplaineddiscretecutoffscanwithtolerancenotyetselected. FinishCURRENTthree
cutoffplots first; no automatic error-tolerance study or extra optimizer job.
Show generated plots before bookkeeping. Lastjobidle/no lease; PIDsreconciled.

User identified the prior error: retain DC and negative-frequency conjugates
FREE; N counts positive components only. Corrected9plots are GENERATED AND
SHOWN TO USER. They use rest10/ramp5/dwell100/ramp5/rest10, f0.01/.1/1Hz,
phases0/45/90degrees DIRECTLY at midpoint65s, and N10/20/40 (columns).
No optimizer or rescaling. Real target a_even*cos(2pi*f0*(t-65)+phi).
Positive cosine C(f)=e^(iφ)G(f-f0)+e^(-iφ)G(f+f0), DC=G(f0)cosφ.
Retain largestNpositive plusNnegativeconjugates and1DC (2N+1coeffs).

Current plots docs/images/wave-real-quintic-n-time.png, -spectrum.png, -phase.png.
Immediate plot checkpoint8f2bf13 pushed after user objected to publicationdelay.
Time rows are frequencies, columnsN10/20/40. Phasesreferdirectly65s. Images
were NOT inspected. Show completed plots immediately; do not delay them for
documentation/GitHub bookkeeping. Earlier complex-projection plots are WRONG
for this reconstruction task, not a physical lower bound. Preserve ashistory.

Corrected source docs/research/real_quintic_packet_example.py and narrative
docs/research/pulse40-three-examples.md. Fullcoeffs/metrics/identities:
docs/experts/numerical_methods_researcher/real-quintic-packets-2026-10-27.json.
N40 peak REAL errors0.0044874833/0.0127837423/0.0568887706m perunit1m;
previouswrongerrors0.2448370989/0.0998263302/0.0596923459m. LargestNrealFourier
selection is L2-optimal on declaredgrid before truncation, not peak-optimal.
Fullfinitecandidate reference errors0.410/0.242/0.206mm; sampled fullFFT
roundtrip1.22e-15m; independentquadrature-vsFFT selectedcoefs3.34e-14max;
inverse/directagreement9.26e-14max; imaginaryresidue<=5.62e-16.
The remaining error is ordinary finite-N truncation. Do notplot an analytic
signal's magnitude as the prescribed envelope at lowf0; show actualrealwave
versusrealquintictarget and target envelope guides. PreserveDCmean exactly.

Root alone ran correctedpulse-real-quintic-attempt1 in6.7850583s (guard8.156s).
Independent source/actual reviews accepted numerical consistency, not physical
quality. Producer12880/owner11776/guard32120/driver25300/conhost8048/python31884
hostabsent; terminal0/frozen/guard0/driver0/Jobempty; finishterminal0/idle/null/
lease release verified. No active or next numerical job. No canonical/production
CLI/solver/WEC-fork changes. No upstream activity or merge.

Accepted inputs:f0,phi(midpoint),rest,ramp,dwell,Npositive. Fixed quintic profile,
direct realFourier projection, largestcoefficients. Bandwidth/frequency placement
and achievedquality are outputs. No additionaloptimizer or quality target.
Any further work follows latest user request, not queuedhistoricalstudies.

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

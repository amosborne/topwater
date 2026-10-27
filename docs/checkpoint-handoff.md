# Checkpoint handoff

Updated 2026-10-27 UTC. Latest user steering controls scope. Historical studies
and queued assignments do not authorize execution. Follow AGENTS.md/coordination.

## Active goal and acceptance

Resolve #31: maximum numerical water-elevation difference PLUS reference
uncertainty <=0.010m against the chosen LINEARIZED model, full60s and declared
valid exterior output region (5929 targets/751 frames, including continuation).
This does not measure error against nonlinear physics or experiments. Preserve
accepted motion/forcing, mean hull, incident/diffraction/radiation, phase and
coordinate/time conventions. No fitted motion, sourceK/mesh correction or WEC
rerun. Native dynamics remain an incomplete prototype. Priority #31, then measured
end-to-end speed, then #32-34. Upstream contributions are deferred.

Active branch baseline/fluid-body; implementation checkpoint d3630fe is pushed.
Full-H offline surface
implementation is independently reviewed and its42 focused tests PASS17.81s.
ONE production readiness gate remains FALSE pending whole-region qualification.
Canonical result/video/log files are unchanged; no corrected output is accepted.
CLI solve documentation describes the pending policy and closed readiness gate.

## Execution state and measured pilot

CURRENT EXECUTION: reviewed FULL RETRY2 is LIVE, root tool session62749.
Manifest artifacts/workflow/issue31-full-attempt2/manifest.json; owner/root,
producerPID30120, ownerPID24388, started2026-10-27T17:01:13.403434Z.
Guard1200s/3GiB/one thread, chunk128/internal1140. Do not launch another job or
edit any consumed sources/helpers/proof notes/data/runtime until exact terminal,
frozen-input/host identity checks, idle readback and release-last complete.
Poll SAME tool session; a timeout is not terminal. No readiness/output promotion.
Variant2 source review48c06635 and actual tiny review718a6b46 passed. New root
policy2 SHA555e05d3 is immutable, with same C/physics/10mm; oldpolicy1 preserved.

The cheap cached128 check PASSED3.719s: combined12.3233→0.731152mm solely by
tightening the input-tail bound11.7158→0.106155mm. Actual candidate/reference
max0.009218mm, C and all other U unchanged. Independent review718a6b46 recomputed
every751x128 value. Original/folding/offgrid/scaling controls pass. This is not
full-region qualification. Cleanup52668 terminal0/frozen/CIM/idle/released before
the current retry. New helpers/full runner are under issue31-whole-surface2/;
runtime producer remains source235989db, readiness FALSE.

MOST RECENT: full-attempt1 stopped after128/5899 queries on its prospective
combined-budget gate, NOT an execution failure. Report at
artifacts/workflow/issue31-full-attempt1/results/report.json and chunks/.
Actual candidate/reference max ONLY0.009218mm, but combined12.323328mm because
the analytic triangle tail bound alone reached11.715798mm. Do not relax10mm.
No common diffraction/geometry/full-region acceptance. Child158.921s includes
140.436s unprofiled prelude before18.485s first chunk; retry must time that.
Independent actual reviewb32f412cbd17189499106a28a2fb198dafdcc9421d68bfd0b1e61c5bf4394dac.
Terminal0/frozen/guard0/driver0/Jobempty; hostCIM owned absent, finish92199
terminal0 separateidle/readback/release-last. Ledger idle/null, lease absent.

The accepted NEW variant preserves old sources/policy1. NMR verified a tighter
INPUT-ONLY tail enclosure, same C/model/10mm/runtime:
evaluate exact Hermite derivative-jump Fourier transform on aligned bins from
floor(48/h)*h to ceil(512/h)*h, h=2pi/8192 (prospective FIXED resolution, not a
runtime-period change). Numerators use periodic FFT of native-knot jump arrays;
denominators use actual unwrapped omega. Add magnitude-trapezoid plus rigorous
complex-interpolation curvature enclosure (B-A)*h^2*sup|V''|/12, off-grid phase
allowance, and original analytic jump remainder beyond the actual upper bin.
No new BEM frequencies or physics/settings changes. Independent pulse/direct-bin/
folding controls, exact source and actual cached-result reviews PASS. Root
created policy2 bound to final new helper; policy1 stays immutable. The current
full retry uses those accepted inputs as stated above. No runs delegated.
Full goal remains open; all consumed dependencies remain frozen while live.

Earlier tiny Hermite polynomial controls PASS1.203s, pulse6.31e-16 and multi7.47e-15
transform errors. Evidence issue31-hermite-attempt1/results/report.json;
finish62855 terminal0/frozen/CIM/idle/released. They verify the original jump
formula, not the forthcoming tighter numerical enclosure or hydrodynamic C.

No numerical job allocated. Earlier pilot
terminal0/frozen/guard0/driver0/Jobempty; exact hostCIM owned identities absent;
finish session92978 terminal0, separate idle readback and release-last completed.
Revalidate liveness before another allocation; recorded state alone is not proof.

Pilot evidence: artifacts/workflow/issue31-whole-surface1/pilot-results/
(report.json, source-bank.npz, pilot-evidence.npz). Actual child62.641s:
source warmup50.484s, cached128-query field5.813s, productionFFT1.093s, six controls
0.547--0.906s each.5899 unique exterior queries map5929 targets;384 are continued.
Bank1920x17x6 retains11520 RHS. Actual peakJobcommit543100928B, sampled aggregate
WS491053056B, complete sampling/no gaps; reportPID2080 matches guard identity.
Independent inverse-sum error8.67e-16m. Source-bank and exact mesh identities pass.
Actual review: issue31-whole-pilot-review1/acceptance.json, SHA0757ef55e8713520eab7bfa7c4e8ec70e04290068ceeeae0ad619920f77fe9c7.

Pilot control maxima(mm): period.00520,dt.02065,final-future.95667,band36.02172,
spacing.05=.14029/.1=.70794, rawstartup.04970/pre.05264. No cheap rejection,
but pilot is not accuracy qualification. Field-only extrapolation267.9s omits
all reference controls/common terms/geometry. NMR recommends prospective1200s,
3GiB,chunk128/internal1140, with first-chunk timing/accuracy rejection and saved
checkpoints. Root accepts these as prospective caps, NOT a launch allocation.

## Immediate authorized work and stopping points

CAP is authoring fresh private qualify.py and hermite_tail.py under
artifacts/workflow/issue31-whole-surface1/, plus concise expert notes. Source-only;
no run, production changes or readiness promotion. Reuse the exact pilot bank.
Explicit bounded extra work is permitted in the eventual reviewed full runner:
SAME-operator64/96rad/s at all actual queries (12 RHS total plus2 field frequencies
per chunk), to calibrate the pointwise response envelope. This is not a band,
mesh or source-physics campaign. Stop at source ready for independent review.

NMR owns independent method/source/results review; never correction authorship.
Its current review task is complete, so use followup_task (not send_message alone)
when CAP source arrives. TOP/receipt_reliability is idle after pilot preparation;
no expert has numerical execution authority. Root owns all serial allocations.

Whole qualification requirements:

- Actual runtime candidate stays intact full-H FFT(v), raw startup retained,
  linear interpolation with every native velocity knot exact, zero prehistory,
  explicit20s Hermite future,1024s period/.01s FFT step/48rad/s/.025 spacing.
- VALIDATION reference may separate R=H+i*omega*A on C1-Hermite(v,a), plus exact
  saved A*a at output knots. Derive A from the chosen finite-low operator;
  native symbolic infinity switches operator and is not its endpoint. Never fit A
  or add it to runtime. Measure candidate/reference difference pointwise.
- Conditional remainder envelope C(point,mode)=2max_known|R| must use dense36--48
  AND64/96, not discard known growth. The independent Hermite jump bound covers
  ALL omitted regular continuum above48, including above Nyquist, under that
  explicit envelope assumption. Method note below. No fixed0.5mm allowance.
- Budget pointwise padding,FFT step,frequency spacing(.025/.05/.1),band,future10/20
  and prefix-preserving20/40 cuts, input interpolation,startup and tail uncertainty.
  Group correlated controls. Independent Fourier sums test implementation.
- Verify common incident/diffraction from unchanged train/operator against saved
  incident and saved TOTAL using old cosine24 only as an identity control.
  Do not merely assume cancellation. Preserve trace mapping and output geometry.
- Cover every declared valid target/time, retain worst location/time/component
  evidence and enforce max difference+uncertainty<=10mm. No near-waterline exclusions
  or seven-probe extrapolation. Existing physical calm/startup/contact/rim/export
  gates still apply. Stop early on a measured failure or infeasible resource cap;
  preserve checkpoints, never tune physics/motion to pass.
- Independent actual-result review, regressions/log triage/docs and appropriate
  canonical result/PNG/video publication remain before completion/gate promotion.

Independent tail derivation/controls:
docs/experts/numerical_methods_researcher/issue31-hermite-tail-2026-10-27.md,
SHA39973a4e6165a50128dcd4faf972b523f957b5c92993ac763938848d34e5e1d0.
For compact C1 piecewise-cubic v, integrate derivative-jump Fourier bounds under
|R|<=C: C/pi*(|a0|/Omega+sum|jump v''|/(2Omega^2)+sum|jump v'''|/(3Omega^3)).
Include outside rest endpoints, native knots and future joins. Exact saved direct
term and actual runtime/reference semantics difference are accounted separately.
CAP implements; NMR reviews exact source/independent cubic controls before use.

## Exact inputs and preserved evidence

- NEW motion artifacts/workflow/reference-adoption-case3/native/motion.mat,
  SHAaf127598fd6bcf3b8bd9bbe727e10f8109df8f9a73f7d87411aa18f69f9c31df.
- Completed field artifacts/workflow/reference-adoption-private-evidence1/
  field-work/surface-recovery/completed-field.npz,
  SHA2000203da3645e3861a667cfd0fe0ee287bc352a3cc8b8a2f3367ba30bc63ec7.
  Native q/v/a payloads match accepted snapshot;3001x6 in fixed mean native basis,
  not Euler rates/world-velocity substitution. No source bank originally retained.
- Actual geometry/native.mat come from native-fullband-producer-attempt1/
  coefficients/operator-manifest.json, SHA f1f788f9; exact mesh/context identity
  checked by pilot against reference-adoption-private-evidence1/surface-context.json.
- Production source235989db/tests544ec0fa+19595990, independent reviewa75ccbde;
  test evidence issue31-offline-tests-attempt3/terminal.json and allocation3 logs.
- Seven-probe actual old/full-complex difference42.660mm (directnear31.927mm);
 24-to48 extension changes height.05916mm. Causal-only cut REJECTED194.933mm,
  despite passing code controls. Do not resume kernel tuning.
- Input-tail screen actual0.188s, seven necessary conditional tests pass,
  max assumed full-Htail+known U=4.457671mm. Actual review02c29fa6. No48-to96
  band extension justified by its rule. It is not whole-field qualification.

Historical policies, all failed studies/setup attempts, full hashes and WEC fork
adoption/cleanup records are preserved in
[the archived pre-qualification handoff](experts/history/checkpoint-handoff-2026-10-27-before-full-qualification.md)
and the Capytaine diagnostic note. WEC runtime remains6b774228, standard public
interface on v7.1.0; fork cleanup complete. No MATLAB window changes or body job
needed for this phase. Never overwrite old failed evidence or canonical artifacts
before acceptance. No preview self-review during actual rendering.

Commit/push substantive checkpoints to baseline/fluid-body, update issue30 main
description and handoff. Preserve user wording and others' changes. No main push,
assistant merge, bookkeeping PR, upstream issue/PR or outside message authorization.
Keep #31/goal open until full acceptance evidence exists.

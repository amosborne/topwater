# Seeded actual-wave geometry audit: incomplete first attempt

One allocated job ended at the unchanged300s compute budget, terminal0 as a valid
INCOMPLETE diagnostic. It does not accept full-case wetting or linearity. Inputs
are the new seeded actual trajectory only; rejected held-release motion excluded.

## Identity and execution

Exact reviewed runner5ea25c1422e5529f724bbb7802b7dae55167be66068dbc7ed7fd50931c567af4,
core7dda063b57328c0e72b6f0fe390d9a8ffd92724cfc37d0051fc39c3af68cb104.
Only READY status changed in seeded config061e685998fd9c02448aae02380b0638f32a48f6b7ee094e48309489402a402e.
Frozen proposal6a42089ed2d0068712325677598400a0b48517ffc892eb20110e36dd094c602b,
source commit e5dc02a1181913c3f1a421d4fb034afeb9591c97. Actual motiond69652e3 unchanged.
Output `artifacts/reference/actual-wave-geometry-seeded-20261024/attempt1`:
manifest5abedd837b01b368586ae9405fd070ea33e1fdbb37051c5b63fb3092070f6527,
geometry6bc8bea487932835538a1c5dfbfe1203c1f8325ec53ab0f1357296a67c9b0196.
postcheck.json includes exact input/output checks and prefix ranges.

Claimc85f71 exit0, separate Verify9daa8c exit0, Run session4488, producer33716/
guard3436, UTC launch19:51:20.979274. Exact sole child waited terminal0 and tool
terminal6d1c72 exit0; both identities absent after terminal0e81a6 exit0. Declared
compute budget300s; observed terminal-log lastwrite minus launch301.5357s (includes
setup/final output, not measured CPU time). Idle48f0f0 under lease, SEPARATE idle
proof39a20b, final removal6c740b, all exit0; ledger not written after release.

## What is known and unknown

All prepared numerical fixtures executed PASS: exact flat-plane triangle clip,
independent mean wet area13.979548845462237m2, canonical reflection parity and
global/native analytic-field closure1.6654e-16m. All frozen input identities and
exact processed trajectory rows verified after terminal.

Only123/3001 saved times0..2.44s were processed. Samples123..3000 (2.46..60s),2878
samples, remain UNKNOWN.122 samples met0.1%mean-area and1mm-rim reporting targets;
the final sample retains valid coarse wet bounds with rim UNKNOWN when budget
expired. All122 evaluated rims were dry at those saved times. No continuous-time
proof or inference about the unprocessed trajectory follows.

Across the processed prefix only:

| Diagnostic |Reported envelope |
| --- | --- |
| Physical wet-area bounds |12.85372..14.93444m2 |
| Material wet-set turnover fraction bounds |3.0184..9.8047% |
| Evaluated rim-clearance bounds |0.46776..0.59529m |
| Max absolute body/local-water roll,pitch |1.2738,1.1816deg |
| Max absolute CG-relative heave |0.033538m |
| Maximum footprint-curvature lower/upper |0.23573/0.48052m |

These are envelopes over prefix bounds, not single-pose values or complete-case
extrema. Footprint-curvature uncertainty has no adopted target and remains broad.
The roundoff allowance is engineering conservatism, not interval arithmetic.
No forces, nonlinear trajectory, incident wave changes, BEM or extra job occurred.

Independent executed prefix/integrity review requested. Full3001-time geometry
applicability remains unresolved even though the seeded native motion is finite
and independently checked. No automatic precision relaxation or rerun is authorized.

## Measured cost and ONE bounded recovery proposal

Measured740489 triangle evaluations and2238572 analytic wave-point evaluations
within the300s compute cap: approximately2468 triangle evaluations/s and0.41 saved
poses/s, roughly6000 triangle evaluations per processed pose. Straight extrapolation
of this implementation would take about2h for3001 samples; that is not proposed.
There are no per-kernel timing counters, so a measured dominant kernel cannot be
claimed. Source inspection points to repeated scalar polygon clipping/area sums,
leaf-list scans/sorts and repeated spectral-constant reductions as likely overhead;
wave evaluation is already batched, though only32-triangle refinement groups are
normally passed to it. This attribution is inference, not profiling evidence.

Cheapest ONE recovery proposal, requiring a new bounded source/review grant:

1. Prospectively use1%mean-wet-area interval width (0.1397955m2) and10mm rim width
   as numerical reporting precision. Original0.1%/1mm choices are not user physical
   requirements. Physical wet-set/righting assumptions or acceptance gates do not
   change, and the old narrower-prefix results remain preserved.
2. Batch the same affine lower/upper/XOR polygon clipping with bounded padded
   vertex arrays; cache spectral constants and vectorize per-triangle diameter/
   roundoff calculations. Preserve scalar clip/area routines as an independent
   oracle and compare deterministic actual-input subtriangles/flat-plane fixtures.
   Keep analytic Hessian/set/rim inequalities, physical faces and frame unchanged.
3. Add low-overhead field/clipping/adaptation timing counters. Use one later
   allocated300s/1CPU job for fixtures and every saved sample in order, unchanged
   depth8/2Mtriangle/8192active/2048point caps. Explicit unknowns if still incomplete;
   no benchmark campaign, decimation, attractive-frame selection or automatic retry.

Reaching all3001 poses in300s needs roughly24x effective throughput versus this
attempt. Coarser reporting bounds plus batched clipping may supply that gain, but
it is unmeasured and not guaranteed. Source math/oracle review must precede launch.
No recovery implementation or numerical job has been performed or allocated.

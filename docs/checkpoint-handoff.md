# Checkpoint handoff

## Latest: cold2772 fields complete, output integration fails three couplings

Read docs/fixed-field-2772-001.md and prospective protocol first. Branch
baseline/fluid-body; producer858bba6. Overall goal ACTIVE/INCOMPLETE: no qualified
source-mesh/finite-wave reference, coupled Topwater motion or usable runtime.
Previous goal turn is concrete progress: actual saved fields now expose the
output error that must be resolved before larger source jobs.

- Exactsession94492 TERMINAL exit1 expected gate failure, no exception. Complete
  true, native_calibration_passed true, integration_passed FALSE, passed false.
  ZERO BEM solves, three full both-half Gauss12/16/20 fields retained.
- Native S*sigma/load calibration1e-10 passed5.94284s. Actual fields took
  102.8695/174.0592/265.4823s, total558.1784s<900. Peakworking234,070,016B,
  peakcommit1,647,284,224B<2GiB; disk186,309,340B<2GiB. No native job live.
- BOTH integration increments fail exactly[1,5],[3,5],[5,1] (zero-based modes:
  sway-yaw,roll-yaw,yaw-sway). Ratios3.175457/1.646947; trend passes. All36
  diagonal/coupling comparisons retained; no symmetry fitting/zeroing/relaxation.
- Uq remains empirical max of both increments. One charged1232->2772 source
  interval fails15 entries,max421.0975, DIAGNOSTIC ONLY; unresolved integration
  prevents source-error attribution. No two-level/source/full qualification.
- All294tests pass8.63s (session5629 exit0), JUnit archived. Verifier reconstructs
  all native/order loads, links field loads to gate arrays, checks imaginary
  loads, all36 gates/Uq/masks and interval metadata. Fresh restoration verifies
  complete actual graph, exact arrays and external inputs. Independent postrun
  audit passed; preserved FAILED evidence, no physics promotion. Evidence fixed-field-2772-001.zip
  SHA5c24ecbfba6a8ff839c77e5fb33085b9f729c72b9c108b329a950050e4523924;
  resultdoc lists restore command/external inputs. Failed result preserved.
- Source-export002 still provides24 columns across four six-mode bases, all
  native A exactly reproduce old raw. Producer242c839,18new+6historical,
  66.4498s, source-only8GiB cap. This is not the field cap or source convergence.
  Old1232 G12/16/20 conditional output gates pass (producercde8f2c); all original
  native/finite failures remain failed. See prior resultdocs/source packet.

**Exact next bounded goal:** FIXED2772 OUTPUT CONTINUATION ONLY, G24/28 plus
verified reuse of priorG20. First implement executable controls/protocol and
prelaunch adversarial review. Scope: same saved002 level2 six-mode sources,
same308seed physical output surface, zero BEM solves. Native calibration repeated
in a fresh cold process at1e-10 scaled. Regenerate/verify priorG20 points, owners,
normals, weights and fields' exact producer/input/payload identity before reuse.
Independently evaluate both halves for newG24/G28 targets177408/241472, total
1,161,135,360 new field pairs,100-target chunks, one native thread.

Prospective1500s TOTAL walltime (imports/inputs/calibration/verification/packaging
included),2GiB peak working AND commit,2GiB accounted disk. Max actual measured
rate8.36720097e-7s/pair: (rate*newpairs +5.94284s native) *1.25 +60s =1281.8602s
<1500; this screen is
not an accuracy guarantee. Model retained prior/new fields and input/archive
copies before launch; keep the prior62.6MB packet as an explicit SHA-bound
external dependency if needed rather than recursively duplicating payloads.
Fresh host margin/identities, minimal selected input loading, archive actual
producer/physical graph BEFORE native work; immutable until exact process terminal.
Use existing shared evaluator/consumer/streamed checkpoints/payload verification.
Do not alter source reconstruction or run a larger source mesh at this stage.

Acceptance: SAME all36 relative0.001 +1e-6 normalized entrywise policy for BOTH
20->24 and24->28 increments, latest<=previous+1e-6; explicit Uq=max both. Preserve
failed12/16/20 evidence unchanged; higher resolution resolves a failed gate,
not relaxation. Require native calibration, finite complete actual fields,
resource/identity/ownership checks. No Richardson or continuum-certification
claim. Single1232->2772 comparison with BOTH Uq remains diagnostic only even
if this integration control passes; no source/full reference promotion.

Deliver actual new fields/matrices/uncertainty/masks/resources, source and reused
field identities, independent postrun audit, durable evidence with fresh full
restore including bound external packets, all tests, checkpoint commit/push,
issue30 update and fresh handoff; THEN STOP. Failure requires review before more
work; no automatic4928/7700 job. Their field cost/quadrature controls must be
separately assessed only after this output control passes.

Exclusions: finite solves, remesh/formulation/major speed campaigns, source/full
promotion, coupled Topwater simulation. All finite A/B/complex excitation/phase
and24 missing diffraction exports remain open. Broader time-domain target still
>=1s matched mild-sea histories within1cm/0.1deg including uncertainty, then
<=200 wallseconds/s including normal preparation. Reject5000x routes early;
major runtime optimization waits physical validation. Offline Capytaine cost
is separate and establishes no Topwater simulation throughput.

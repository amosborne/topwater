# Cold2772 field evaluation completes; integration gate fails

Producer `858bba6ef23e7017094fbc7758887783afc7bdcb` ran the
[prospective protocol](fixed-field-2772-protocol.md), using saved infinite-frequency
source-export002 data and zero BEM solves. Exactsession94492 is terminal exit1
because integration failed; execution itself completed without exception.
Complete=true, passed=false, native_calibration_passed=true, integration_passed=false.
Source/full-reference qualification and full_run_authorized remain false.

Native calibration took5.94284s, reproducing saved potentials and all36 native
loads at1e-10 scaled. Potential/load tolerance ratios are2.14655e-5/2.31320e-6
(fractions OF that gate, not physical relative errors). All target coordinates,
owners, normals and weights equal the admitted common308seed rule. Both halves
were evaluated independently; every actual complex field and pressure matrix
is retained. Native7,683,984 plus field683,020,800 source-target pairs ran.

Gauss12/16/20 took102.8695/174.0592/265.4823s. Total558.1784s<900s; peak working
234,070,016B and committed1,647,284,224B<2GiB. Accounted disk186,309,340B<2GiB.
These are offline reference costs, not Topwater simulation timings.

Both required integration increments fail the unchanged all36 entrywise policy:
maximum ratios3.175457/1.646947, trend passes. Exactly the same three couplings
fail BOTH increments; diagonal entries pass. Physical values and final allowances:

| Response, acceleration | G12 | G16 | G20 | Allowance |
|---|---:|---:|---:|---:|
| Sway, yaw (kg m) | -0.52848340 | -0.54656847 | -0.53718867 | 0.00569527 |
| Roll, yaw (kg m²) | 31.28561872 | 31.36211373 | 31.32027930 | 0.03511747 |
| Yaw, sway (kg m) | -0.89840137 | -0.91288534 | -0.90563111 | 0.00606371 |

The output uncertainty remains the empirical maximum of the two increments,
not a certified continuum bound. The charged1232-to2772 single-source-interval
diagnostic fails15 entries, max421.0975. Because2772 integration is unresolved,
this does not separate source error from remaining output error. Even agreement
would not establish source convergence. All earlier native/finite failures stay.

The evidence verifier independently reconstructs all36 native and each order's
loads from retained fields, links their real stack to quadrature_A, checks
imaginary loads, recomputes gates/Uq/masks and charged source-interval metadata.
Fresh extraction/restoration verifies exact arrays, complete producer graph and
all physical input archives. All294 local tests passed, zero skips/errors/failures
(session5629 exit0, JUnit in packet). Independent postrun review verified
producer/input/dependency/source/coarse identity, actual fields and all pressure/
gate/Uq/mask/interval calculations, payload and resources. No material blocker;
approved as preserved FAILED integration evidence, with no physics promotion.

[Durable packet](../validation/evidence/fixed-field-2772-001.zip),62,585,567B,
SHA256 `5c24ecbfba6a8ff839c77e5fb33085b9f729c72b9c108b329a950050e4523924`.
Adjacent summary retains all36 Uq and failure masks/diagonals/resource records.
Actual field NPZ/source are in immutable payload; source-export packet (including
its old pilot), preflight and physical JSON are embedded. After extraction run
`python restore_inputs.py --frozen-bundle PATH --raw-bundle PATH --sources-bundle PATH`
using the tracked archives and exact SHAs in external-inputs.json:

- frozen-fixed-001.zip: `fcaf28963be08e4443410b8c19d680f1fb81cff3c44a291c34aa69c6e1a9847d`
- reference-study-003.zip: `e18bf7366eb651766edc5487addd84c798c1ae53183f58a0e88b3b623ed6729d`
- periodic-sources-003.zip: `4a317cbb91cc51348b0059b40b4a6437315e3656624ad2903743dec4721214e5`

Next proposed bounded correction: fixed2772 output continuation G24/28, reusing
EXACT priorG20 field/targets after full identity verification; final sequence
20/24/28 under unchanged gates. Preserve failed12/16/20 separately. New targets
177408/241472,1,161,135,360 pairs, zero solves and repeated native calibration.
Maximum measured seconds/pair8.36720097e-7 gives projected1281.8602s including
repeated5.94284s native calibration, factor1.25 and60s preparation reserve,
below proposed1500s total; retain cold
2GiB process/disk caps and one thread/100-target chunks. This is a prospective
cost screen, not an accuracy prediction or runtime qualification. Source and
full reference remain false even if it passes. No larger mesh is authorized
automatically. Full scope/acceptance/stopping checkpoint is in the handoff.

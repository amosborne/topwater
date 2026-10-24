# A-priori geometry operating-envelope study, 2026-10-24

Scope: static geometry/buoyancy/restoring curves on the unchanged canonical full
hull/mass/CG. No BEM solve, motion trajectory, native nonlinear adapter, Topwater
solver change, bespoke beam-sea optimization or multi-case dynamic validation.
The user now wants ONE future buoy-generator sea/initial-heading case with roll
5-10deg and appreciable other motions. This static curve family precedes that
choice and does not alter the trusted minimum baseline.

Implementation: reference/capytaine/geometry_operating_envelope.py; controls in
reference/capytaine/geometry-envelope-controls.json. Native Capytaine3.0.0 Mesh
GL2 integrates the physical hull after exact planar clipping. Artificial deck
closure faces are removed. Opening-rim contact/submergence invalidates physical
closed-hull buoyancy BEFORE any cap volume can be counted.

Pose: R=Ry(trim)*Rx(heel) about mean WORLD CG, then vertical shift h. Waterplane
is z=0; h is relative immersion displacement against this test plane. Buoyancy
is rhogV; net Fz=rhogV-mg and world torque=(CB-CG)cross(rhogV*ez). Compare the
full6load vector with -Kmean*[0,0,h,heel,trim,0], retaining physical coupling.
Also report torque projected onto the heel axis. Free-heave/trim righting curve
solves V=m/rho and CBx=CGx at each fixed heel. No symmetry is imposed.

Frozen proposed coverage:

-49heel points -12:0.5:12deg, each fixed-heave/trim and free-heave/trim.
-175additional static poses: heel[-10,-5,0,5,10]deg,
  h[-.25,-.1,-.05,0,.05,.1,.25]m, trim[-5,-2.5,0,2.5,5]deg.
-Free equilibrium bounds h+/-0.5m, trim+/-15deg; maximum40least-squares nfev,
  plus finite-difference Jacobian calls (all count in recorded geometry evaluations).
-One thread, expected under60s, declared300s soft budget, zero BEM/dynamic jobs.

Metrics: displacement volume/CB/buoyancy/full loads, equivalent arm (GZ only on
mass-balanced free equilibria), physical wet area, waterplane area/firstXY/full
secondXY moments, free-heave/trim offsets, rim clearance, waterline vertex and
geometric knuckle-edge contacts. Material wet-set symmetric difference tracks
port exit AND starboard entry on original physical faces; aggregate net wet
area alone can hide the turnover. Ordinary chine crossing is diagnostic, not
an automatic rejection. Actual opening immersion, invalid wet body, or a
nonrestoring branch rejects that physical state. Failed solver/invalid trial
leaves equilibrium existence unknown; it does not prove physical nonexistence.

Prospective restoring-shape labels use effective free-heave/trim roll stiffness
Keff=Krr-Krh*Khh^-1*Khr, with h=[heave,pitch]. Esec=abs(R-Keff*phi)/(Keff*abs(phi));
Etan=abs(dR/dphi-Keff)/Keff. Meanphi0 uses the independently checked derivative
limit. Report GZ/(GM0*sinphi) separately to explain sine-versus-angle effects.
The 5%Esec/10%Etan labels are explicit engineering assumptions requiring
coordinator adoption for an operating-range decision; no hard wet-area or
wet-set percentage is inferred. Report the connected positive/negative range;
only the contiguous converged component containing zero receives tangent-based
labels. An outer failed root cannot erase valid inner coverage; no disconnected
island bridges a gap. Synthetic49point checks independently verify these rules.
These criteria do NOT establish finite-wave hydrodynamic or trajectory accuracy.

Future buoy cases must interpret heave relative to local water elevation and
body orientation versus local wave slopes. Large global heave alone is not a
failure/flooding criterion if boat and water rise together. Calm-water heel
curves diagnose fixed-mean assumptions; no curved free-surface nonlinear model
is introduced here. The actual sea's entire spectral support remains governed
by the separate shared contract; this geometry study neither filters nor fits it.

Independent expected checks (part of the one proposed job):

- Wall-sided box at three vertical shifts: exact V, CB, wet area, waterplane.
- Every valid native GL2 pose independently checked against clipped-polygon
  triangle-fan analytic area/first/second moments (separate clipping/integration).
- Mean physical full hull V/CB matches accepted native-gl2; centered numerical
  force derivatives agree with full accepted heave/roll/pitch K block.
- Free equilibrium residuals explicitly bounded; all +/-heel parity evaluated:
  h/trim/area/Ixx/Iyy even, CB_y/righting/Ixy odd, wet-set turnover even.
- No output replaces baseline hydrodynamics or motion. Invalid/flooded cases and
  failed roots remain explicit rows, never accepted/artificially symmetrized.

Expected artifacts: manifest/source snapshots, geometry-results.json (fixed roll,
free roll and pose-grid tables), geometry-curves.png with six static plots. A
per-pose numerical consistency check is separate from engineering applicability.

Command from C:/Users/alex/topwater (fresh output directory under allocated lease):

.venv/Scripts/python.exe -u reference/capytaine/geometry_operating_envelope.py --case reference/capytaine/canonical-boat.json --native artifacts/reference/minimum-20261024/native-gl2 --controls reference/capytaine/geometry-envelope-controls.json --output artifacts/reference/geometry-envelope-20261024/attempt1

Original proposal is retained above; executed evidence follows below.
Independent researcher math review is
../numerical_methods_researcher/geometry-operating-envelope-review-2026-10-24.md;
source signoff passes after the connected-coverage interpretation correction;
independent source/protocol and executed review now PASS within the static scope.

## Execution and protocol

Attempt1 terminal exit1 in0.26347s: abs() received a Python-list residual; no
envelope was accepted. Its manifest/log/source/proposal remain preserved. Protocol
deviation: atomic lease was held, but an existing UTF8 BOM made claimed-ledger
read fail and the launcher improperly continued without that state. It was not
retroactively labeled claimed/running. Exact terminal evidence was written idle
under lease then released; post-release idle/no-lease readback passed.

The numerical fix is np.abs on the vector; diagnostics now preserve mean checks
and tracebacks. Independently reviewed geometry_execution_guard.py uses UTF8-sig,
separate record-claim/verify/run invocations, a matching owner/timestamp/proposal/
lease/frozen-identity receipt, and exact Popen child wait even after a ledger-write
failure. No execution policy or generalized framework change. Exact inline/direct
commands: artifacts/reference/geometry-envelope-20261024/proposal-attempt2.json.

Retry proposal SHA2562f24e09a1c6b8d022ca279dff307f09991eba57caafddc585b5a4f5892d6647b;
physics sourcef5604fa1e8195d436f61b1e5758a88c98411226d7c2157c4ca01933ccb5fd4be;
guard1d85dd2a94391a48fd48664b66e6b23429267e135fb04e99adc0d026b804534d;
controlsf78729aa8a74abe3e6da639f366e3f1fa18d9aaf077a8b781a975ce04eb0be3a.
Separate Claim243d14exit0 and Verifya11457exit0 preceded Run; receipt is
attempt2-claim-readback.json. ProducerPID40832/guard26852/session33660 terminal0
chunk593ac4. Elapsed15.2183183s/951geometryevals/one thread/zeroBEM/zero trajectories.
Idle written under lease, SEPARATE under-lease readback602ed0exit0, final removal
2ed4b7exit0. See attempt1-terminal-record.json, attempt2-terminal-record.json and
attempt2-idle-readback.json. No additional numerical work is allocated.

## Executed curves and physical interpretation

Output: artifacts/reference/geometry-envelope-20261024/attempt2. All9postcheck
groups pass; all49free-heave/trim roots converge, max normalized residual5.535e-14.
Every262retained valid pose passes independent analytic moment checks; all signed
parities/box/meanK/connected-coverage checks pass. Eleven of175additional poses
have actual opening-rim immersion against the tested flat plane and are flagged
invalid before accepting buoyancy. Free branch remains restoring and unflooded
through+/-12deg. That stable branch does not establish fixed-mean linearity.
Independent researcher executed integrity/interpretation review PASS, including
input/output hashes and visual plot; no reviewer numerical rerun.

| Relative calm-water heel | Exact righting below mean-K prediction | Tangent deviation vs mean K | Material wet-set turnover | Free heave | Rim clearance |
| --- | --- | --- | --- | --- | --- |
| 1.5deg | 3.14% | 7.40% | 5.38% | 0.648mm | 0.562m |
| 2.5deg | 5.96% | 12.83% | 8.53% | 1.655mm | 0.539m |
| 5deg | 12.33% | 24.03% | 15.87% | 5.953mm | 0.482m |
| 7.5deg | 17.73% | 32.63% | 22.58% | 12.537mm | 0.427m |
| 10deg | 22.33% | 39.32% | 28.82% | 21.154mm | 0.373m |

Mean free-roll stiffness41055.66852695Nm/rad. Under UNADOPTED exploratory
5%secant/10%tangent labels, connected scanned range is+/-1.5deg. Relative flat-water
5-10deg states have substantial geometry/restoring departure; no threshold was
retuned to rescue them. At5/10deg net wetted area only changes-4.98%/-8.33%, while
material entry/exit turnover is15.87%/28.82% and waterplane area-7.05%/-11.95%.
Sharp-knuckle/waterline changes explain why aggregate wet area alone was not a
sufficient diagnostic. The effect is not merely sin(phi) versus phi.

This is NOT an automatic GLOBALroll ceiling. Body following a local wave slope
may roll5-10deg globally with smaller relative immersion/orientation. However,
that requires actual generated-case geometry evidence; it cannot be assumed as
a blanket rescue. Future interpretation must use bodyheave-localeta and
body-versus-wave slopes/curvature, not global heave or angle alone. These still-water
curves are conservative a-priori diagnostics, not a nonlinear wave-surface or
trajectory model. No final case gate or dynamic accuracy claim was adopted.

Artifacts: geometry-results.json (full fixed/free/175pose rows), geometry-curves.png
(visual QA passed), manifest.json, postcheck.json (source/output hashes). Next:
coordinator and WEC owner bind any candidate to faithful buoy-bank/heading and
relative-water geometry before authorizing new native hydrodynamics/dynamics.

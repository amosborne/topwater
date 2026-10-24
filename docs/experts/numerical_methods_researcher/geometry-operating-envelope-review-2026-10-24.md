# A priori geometry operating-envelope review, 2026-10-24

Independent read-only assignment; Capytaine owner implements/runs the geometry
study later under coordinator allocation. No implementation, numerical job,
weakly nonlinear trajectory, native nonlinear adapter or Topwater work here.
Current objective is ONE buoy-generated canonical case with significant 5-10
degree roll and appreciable other five motions. The linear model's physical
applicability is an a priori engineering assumption, supported or contradicted
by geometric buoyancy/wetting evidence, not proved by a nonlinear simulation.

## Independent math review

Cap's proposed pose convention is internally consistent with the mean-world
reference: rotate the full mean hull about its mean CG with
`R=Ry(trim)*Rx(heel)`, then translate CG vertically by heave. Keep this explicit;
do not silently exchange world/mean-body axes or rotation order.

For each flat-water clipped pose, exact volume `V` and buoyancy center `B` give

`Fz=rho*g*V-m*g`, `tau=(B-G) cross (rho*g*V*ez)`.

At fixed heel, free-heave/trim equilibrium therefore solves
`V=m/rho` and `B_x=G_x`. Define signed righting moment `Rheel=-tau_x`, positive
for a positive stable heel, and equilibrium righting arm `GZ=Rheel/(m*g)`.
Outside mass equilibrium, label the latter an equivalent moment arm rather than
an intrinsic geometric GZ. Preserve the fixed physical CG/mass; do not trim away
a residual by changing the boat.

Compare the full exact wrench with `-Kmean*q`, where
`q=[0,0,heave,heel,trim,0]`. For the free-heave/trim heel curve, the linear
prediction uses condensed stiffness

`Keff=Krr-Krh*inverse(Khh)*Khr`, `h=[heave,pitch]`, `r=roll`.

The actual lateral symmetry makes the relevant coupling near zero, but compute
it rather than hardcode it. Record equilibrium heave/trim drift: even when
righting is nearly linear, a nonzero second-order equilibrium shift is useful
physical evidence of the approximation.

## Metrics that answer the user's assumption

- Full signed righting/GZ curves over both signs through 12 degrees, explicitly
  marking 5 and 10 degrees, secant departure
  `abs(Rheel-Keff*heel)/(Keff*abs(heel))`, and tangent departure
  `abs(dRheel/dheel-Keff)/Keff`. At zero use the derivative limit. Also plotting
  `Rheel/(Keff*sin(heel))` explains geometry beyond trivial sine-versus-angle
  kinematics; the fixed-K comparison itself uses angle, not sine.
- Wetted area, waterplane area and its full second-moment matrix, CB and
  equilibrium heave/trim, with signed/absolute changes from mean geometry.
- Material-surface wetted-set symmetric difference
  `area(Wheel XOR Wmean)/area(Wmean)`. Net wet area alone can hide substantial
  port entry/starboard exit. Compare partial face polygons in common material
  coordinates; whole-face/vertex classification is inadequate.
- First knuckle/chine/rim events and remaining opening clearance. Ordinary
  waterline/chine crossings are diagnostics, not automatic failure. Rim
  submergence/downflood, loss of a valid watertight immersed hull, or loss of the
  stable restoring branch reject the assumed supported operating region.
  Nonphysical opening caps must never add accepted buoyancy after contact.

The proposed fixed-pose heave/trim coverage is an envelope diagnostic, not a
matrix of simulated sea cases. Its heave is RELATIVE IMMERSION:
`h_relative = body_global_heave - eta_at_body_reference`. An ocean-following boat
may have large global heave while its wetted geometry remains nearly unchanged;
do not reject it because global heave exceeds calm-water draft. Similarly compare
body heel/trim against the local incident-wave slope, and inspect curvature and
waterline variation over the footprint. A local tangent plane has normal
proportional to `[-eta_x,-eta_y,1]`; it does not remove curvature effects.
Flat-water curves do not capture spatially varying incident crests, and no
uniform-wave-level assumption follows merely from wavelength exceeding twice
beam. Preserve generator output and the selected heading convention rather than
fitting another sea.

## Cheap independent checks and tolerance ownership

Use analytic clipped-triangle moments already reviewed for GL2, and independent
tetrahedral volume/first moments on the closed immersed polyhedron. Upright must
recover `V=2500/1025`, physical equilibrium and reviewed mean K. Finite-difference
exact wrench slopes converge to `-Kmean`; the free-equilibrium righting slope
converges to Keff. A rectangular-box analytic upright volume/CB/waterplane and
restoring slope is a small decisive fixture. These are geometry identities,
not new hydrodynamic cases.

Mirror symmetry predicts righting/CB_y odd in heel and equilibrium heave/trim,
wet area and waterplane area even. Do not impose fore/aft symmetry. Check both
heel continuations and root residuals; identify branch/contact events rather
than silently interpolate across them. Tight moment/root tolerances establish
geometry computation, not physical linearity.

Engineering tolerance proposal ONLY for coordinator/user adoption before judging
the curve: label a region with at most 5 percent secant and 10 percent tangent
righting departure as a candidate small-departure region; show the actual
numbers at 5 and 10 degrees. These are explicit engineering screening choices,
not universal linear-regime thresholds or motion-accuracy claims. Wet-area,
symmetric-difference and waterplane changes remain visible supporting evidence;
there is no defensible universal conversion from a percentage of wet area to a
six-DOF error tolerance. Do not invent an unowned dynamic 2-percent criterion.

A candidate case is accepted only when the owner/user explicitly states the
assumed envelope and its limitations. Even a favorable righting curve does not
prove that mean-wetted-surface radiation/diffraction coefficients remain accurate
at 5-10-degree roll. That uncertainty is acknowledged, not replaced by a demand
for a nonlinear model the project does not have.

Stopping point: proposed math/assumption review and subsequent read-only review
of Cap's source/results. No new case or operating envelope is certified yet.

## Initial source review

Read `reference/capytaine/geometry_operating_envelope.py` and its controls without
executing them. Pose/plane transforms, material wet-set intersections, analytic
polygon moments, native GL2 comparison, exact buoyancy/torque, free-heave/trim
residuals and Schur stiffness agree with the independent derivation. The three
analytic box-heave fixtures and upright finite-difference K checks are decisive
inexpensive checks; all valid static poses compare native clipping/integration
against the independent polygon moments. Reflection parity includes wet-set
turnover and correct waterplane Ixy sign.

One required result-validity correction was returned before execution: a failed
outer free-equilibrium point must not erase valid inner 5-10-degree coverage.
Tangents/engineering labels should use the contiguous valid component containing
zero, with unresolved outer coverage explicit. An invalid solver trial or failed
root is not proof that no physical equilibrium exists. Also retain the distinction
between world-x equivalent righting arm and torque about the actual heel axis.
Final source signoff awaits this correction; no numerical result exists here.

## Final independent source signoff

PASS for the bounded static-geometry study source before coordinator allocation.
Reviewed current source SHA256
`d2c395d26f478bb7e5b42bb3b179585ce529d477d66bce1c54c6213c646651c7`
and draft controls
`59eb2dc920f6dceb8c4942bafba9aa205723e5dfb7473fc0b749fc8937fe7219`
(controls status still awaiting owner readiness flip/refreeze).

The required connected-coverage correction is resolved: the component containing
zero retains valid inner coverage despite failed outer points; failed inner
points do not bridge disconnected islands; unresolved zero yields no connected
branch. Synthetic expectations correctly retain +/-11.5 degrees when only
outer +/-12-degree points fail. Tangents and labels use that component, with
one-sided boundary differences explicit and invalid plot gaps left unbridged.

Root/trial failures explicitly leave physical equilibrium existence unknown.
Native-versus-independent moment inconsistency now raises a distinct runtime
error that cannot be hidden as a root failure. World-x and projected heel-axis
arms are separately labelled. Only a residual report wording cleanup was
recommended to keep its short physical-limit description consistent with these
statuses; there is no remaining blocking mathematical issue.

Cap owns implementation and the reviewer authored none of this correction.
No numerical study was executed by this reviewer. Source PASS supports the
coordinator's allocation of the already authorized bounded study; it does not
grant scope or accept an operating envelope or a 5-10-degree generated-sea result.

## First execution failure and bounded correction review

Owner attempt one terminated exit one before producing an accepted curve:
`abs(residual(solution.x))` applied Python's scalar `abs` to a list. Independent
diff review of current source
`f5604fa1e8195d436f61b1e5758a88c98411226d7c2157c4ca01933ccb5fd4be`
passes: `np.abs` performs the intended vector residual absolute value; remaining
changes only persist earlier successful checks and record tracebacks/evaluation
counts on failure. No force, geometry, root or engineering-label criterion changes.

The ledger honestly records a separate orchestration deviation: an atomic lease
was held, but the claimed-ledger read failed on a UTF8 BOM and the orchestrator
continued despite the preparation failure. No concurrency occurred. Future
launch preparation must read JSON as UTF8-sig, stop explicitly on every nested
nonzero exit, and obtain a separate successful claimed-ledger readback before
any numerical launch. Protocol implementation/proposal review is pending;
source math PASS alone does not justify a retry.

Final independent retry source/protocol review PASS: geometry source remains
`f5604fa1...`; reviewed stdlib guard is
`1d85dd2a94391a48fd48664b66e6b23429267e135fb04e99adc0d026b804534d`.
Its UTF8-sig reads and separate Claim/Verify/Run invocations stop preparation
failure from launching numerical work. A successful separate claimed-ledger
readback creates a receipt bound to owner, claim timestamp, proposal hash and
frozen input/dependency checks; Run requires that matching receipt before its
only spawn site. Inline PowerShell commands explicitly return native exit codes,
and each prior tool result must be inspected as exit zero before advancing.
No execution-policy modification or blocked `-File` invocation is used.

One additional source-review finding was corrected before PASS: if writing the
running ledger fails after Popen, the guard now preserves the claim/freeze,
waits for that exact child to terminate, prints its PID/actual exit, and only
then propagates the bookkeeping failure. Completion still belongs to the owner:
terminal evidence, idle write under the lease, a separate successful idle
readback under that lease, and lease removal as the final operation.

At review, draft proposal/control status still blocks execution pending a fresh
coordinator retry allocation and owner refreeze. No retry was executed by this
reviewer and no geometry result is accepted by this source/protocol PASS.

## Final independent executed geometry review

PASS for the integrity and interpretation of the bounded static geometry study;
this does NOT accept a 5-10-degree relative-heel linearity assumption.

Inspected authoritative attempt-two `manifest.json`, `geometry-results.json`,
`postcheck.json`, plot and frozen proposal under
`artifacts/reference/geometry-envelope-20261024`. Recomputed every declared
input and output hash, including archived source/control/physical input: all
match. Claim readback's proposal hash
`2f24e09a1c6b8d022ca279dff307f09991eba57caafddc585b5a4f5892d6647b`
matches the actual frozen proposal. Separate claimed readback preceded launch;
guard waited exact producer PID 40832, terminal exit zero; separate idle readback
was performed while the lease remained held, then the lease was released. The
lease is absent at review. Attempt one and its protocol deviation are preserved.

All 49 free-heave/trim roots converge; maximum normalized root residual is
`5.5351e-14`. All 175 additional static rows are retained, with 11 opening-contact
invalid states explicitly excluded. All 262 valid retained fixed/free/grid poses
pass independent analytic-versus-native geometry checks. Upright K, analytic box,
reflection parity and connected-coverage checks pass. Independently recomputed
secant and centered-grid tangent metrics at 5 and 10 degrees from the stored
righting moments; differences from reported values are below `4e-16`. The
rendered curves agree with the stored results and label their limited scope.

| Relative flat-water heel | Mean-K righting departure | Tangent-stiffness departure | Material wet-set turnover |
|---|---:|---:|---:|
| 5 degrees | 12.33% | 24.03% | 15.87% |
| 10 degrees | 22.33% | 39.32% | 28.82% |

The mean K is too stiff relative to these exact righting curves. Net wet area
changes are smaller (-4.98% and -8.33%), illustrating why wet-set turnover matters.
Free-heave equilibria move upward by 5.95 and 21.15 mm. The free branch remains
restoring and avoids opening immersion through +/-12 degrees; that is distinct
from small-departure applicability. The connected +/-1.5-degree 5%/10% label is
an exploratory engineering label, not a final user-adopted operating limit.

These results contradict describing 5-10 degrees of RELATIVE flat-water heel as
a small geometry/restoring perturbation for this hull. A large GLOBAL wave-following
roll need not equal relative heel, just as global heave need not equal relative
immersion. That possibility is not evidence that the generated sea rescues the
assumption. The next evidence must examine the one selected generated sea's
actual relative immersion, body-versus-local-wave slopes and footprint curvature/
wet-set changes, with frame conventions explicit. Flat-water curves remain a
conditional diagnostic; a tilted-water geometric comparison is not an exact
nonlinear wave-load calculation. Do not fit thresholds, invent another sea,
require a cancelled nonlinear model, or claim dynamic motion accuracy here.

The implementation author and this reviewer are distinct. Reviewer performed
only stored-artifact identity/algebra/visual checks, no geometry/BEM/dynamics job.

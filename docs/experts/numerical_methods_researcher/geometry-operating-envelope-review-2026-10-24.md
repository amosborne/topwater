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

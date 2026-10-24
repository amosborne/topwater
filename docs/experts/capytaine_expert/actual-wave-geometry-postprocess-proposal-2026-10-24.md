# Actual native trajectory: incident-water geometry audit design

SOURCE IMPLEMENTED AND INDEPENDENTLY REVIEWED; no geometry allocation or execution. This follows the isolated lid
control and uses a future actual native January90 trajectory, not particular RAOs
or a fabricated trajectory. It changes no body/sea/force/operator and computes no
nonlinear hydrodynamic trajectory. Goal: applicability evidence about material
wetting and opening clearance under the prescribed analytic incident surface.

## Inputs and frame

Freeze canonical full hull/mass/CG, all native timestamps/CG positions/orientations,
and the unchanged January90/global bank plus explicit heading R/t. Transform mean
world physical hull vertices once to CG-local coordinates, then each native pose
as `v=Rz(yaw)Ry(pitch)Rx(roll)*(v_mean-CG_mean)+CG_native`. Use stored native radians
and raw sensor CG, not displacement plus another q0 rotation. Check t0 matches
the declared actual initial pose. Do not bake initial heading twice.

Evaluate the exact384-component analytic incident elevation/gradient/Hessian at
actual hull/world positions and times. Use either native-frame bank and native
positions or global bank and its exact R/t; cross-check both at test points.
No phase shift, scaling, tail removal, moving origin, wave ramp or elevation
interpolation from a renderer. Label it prescribed incident geometry; the model
does not supply a computed scattered/radiated free surface.

## Bounded geometry metrics

At every saved native output sample, report:

- CG height relative to eta at actual CG, actual body-versus-local tangent-plane
  orientation from relative rotation matrices (not Euler-angle subtraction), and
  full-footprint departure from that plane. Report full surge/sway/yaw drift.
- Physical material wet area, each material face's wet fraction and symmetric
  difference from mean material wet set. Net area cannot hide entry/exit turnover.
- Opening-rim clearance minimum and actual wet/contact/unknown flags. Artificial
  closure faces never contribute physical wet area or conceal an immersed opening.
- Mean-relative pose/restoring-curve comparisons, clearly labeled flat-plane
  geometry diagnostics. Do not turn curved incident geometry into a claimed
  nonlinear buoyancy-force model. Existing independent free-heel curves are only
  one slice; combined pitch/heave cannot inherit their equilibrium assumptions.

The January90 particular proxy's~2.1deg relative roll and~2.5deg relative pitch
are smaller than global~5.2deg roll, but still not negligible: stored free-heel
2/2.5deg gives4.57/5.96% secant,10.19/12.83% tangent departure and6.97/8.53% wet-set
turnover. Those equilibria require only1-1.7mm heave, not the proxy's~30mm relative
heave. Its sampled0.307m footprint curvature departure cannot be dismissed by the
small CG-relative height. No final physical acceptance percentage is adopted here.

## Spatial convergence with explicit unknowns

For each transformed planar material triangle, set f=z-eta(x,y,t); wet means f<0.
Evaluate f at vertices, clip the affine nodal interpolant exactly on the material
triangle, and adaptively split only ambiguous/contact triangles. Do not classify
an entire large face from its centroid.

For the exact analytic bank, horizontal Hessian norm is bounded by
`C2=sum(abs(a_j)*k_j^2)`. If h is projected triangle diameter, the interpolation
error satisfies `abs(f-f_PL)<=C2*h^2/2`. Thus clipping f_PL at-minus/plus epsilon
gives certified lower/upper material wet areas. Apply the same nested set logic
to mean-wet/current-wet intersections to bound wet-set symmetric difference.
Track bounds and successive-level changes per face and summed over the hull;
refine until both explain the reported numerical precision. These are geometry
integration tolerances, not motion-accuracy or linearity thresholds.

For each physical rim segment use exact analytic surface samples and the
one-dimensional second-derivative interpolation bound. Positive lower clearance
proves dry at that saved time; an actual negative sample proves immersion;
straddling bounds stay UNKNOWN until refined. No contact sign is forced to pass.

Prospective numerical reporting target: wet-area interval width<=1e-3 of mean
wet area and rim-clearance interval width<=1mm, with actual achieved bounds
reported. Compare at least the last two subdivisions and retain failures. A
prospective one-CPU/300s cap, maximum subdivision depth8 and2million triangle-field
evaluations prevent an unbounded geometry campaign; uncompleted samples/metrics
remain explicitly unknown. Those resource/tolerance choices require source review
before a later grant. No automatic second job, relaxed tolerance or decimation.

Evaluate all saved samples if the budget permits; record exact coverage and
unprocessed indices if not. No assertion about continuous time between samples
follows without a separately justified temporal bound. Merely checking extrema
of global q or selecting attractive frames does not certify the trajectory.

## Output and stopping point

Proposed output: immutable input/source manifest, per-sample/per-face metrics and
uncertainty intervals, rim flags, relative-pose histories, worst-frame geometry
images with numerical bounds, and explicit incomplete coverage. Parent owns
engineering interpretation; metrics establish neither nonlinear wave-force
accuracy nor a trusted linear operating envelope by themselves.

## Prepared implementation and stopping point

Coordinator authorized source implementation only. Files:
`reference/capytaine/actual_wave_geometry.py`, `actual_wave_geometry_core.py` and
`actual_wave_geometry_config.json`. Runner consumes actual `motion.mat` fields
time, position_world, q and mean_world_cg; it cross-checks displacement/absolute
CG and Euler columns. All saved samples are processed in order or their exact
unprocessed indices remain unknown. Current January90 ACTUAL trajectory is
REJECTED (roll~107.5deg), so no geometry run is allocated and no actual file hash
is bound. This tool must not be used to guess that rejected run's cause.

```powershell
.venv/Scripts/python.exe -u reference/capytaine/actual_wave_geometry.py --config reference/capytaine/actual_wave_geometry_config.json --output artifacts/reference/actual-wave-geometry-20261024/attempt1
```

The command is prospectively specified, currently blocked by draft status and
unbound actual-motion SHA. Future allocation must freeze the two source modules,
config, actual motion/runtime provenance, canonical hull and both frame banks,
plus numerical dependency/Python identities; independent source review precedes
the numerical fixtures and geometry processing. Original trajectory is read-only.

Prepared fixture logic checks an exact tilted-triangle flat-plane area, independent
mean wet area against original native wetted triangles, canonical reflection
parity and global/native analytic field closure. Only AST/syntax inspection has
run; fixtures and all geometry evaluations are UNEXECUTED. There is no synthetic
motion trajectory or alternate body integrator in the fixture path.

Implementation retains one CPU/300s/depth8/2million triangle evaluations and uses
wave batches of at most2048 points. An8192-active-triangle per-sample memory cap
also returns UNKNOWN with existing valid coarse bounds. No tolerance is relaxed.
The analytic interpolation bound includes a reported conservative double-precision
roundoff allowance; this is not a formal interval-arithmetic implementation.
Full-footprint curvature has explicit lower/upper bounds and no adopted numerical
width target; the area/rim targets do not certify its adequacy automatically.

Outputs are checkpointed geometry-results.json and manifest.json, including every
processed sample's per-face intervals, numerical statuses, rim status, relative
pose/eta/gradient, curvature intervals and exact missing indices. Resource/depth
termination can be a valid incomplete diagnostic exit0, never a physical PASS.
Source errors exit1 and preserve checkpoints. No source or phase modification,
physical load, free-equilibrium root, nonlinear dynamics or BEM is included.

Stopping point: source review/readiness only. No further implementation campaign,
geometry launch or application to the rejected January trajectory is authorized.

Independent SOURCE math/frame/coverage review PASS at runner
5ea25c1422e5529f724bbb7802b7dae55167be66068dbc7ed7fd50931c567af4,
core7dda063b57328c0e72b6f0fe390d9a8ffd92724cfc37d0051fc39c3af68cb104,
config45da2e2050f6a47454d25290abea5398b5455a008abed56ebafd1dae6e608a30.
Source-only provenance manifest666cfcaa51200463aa4a0bf820bb69d9b98846417018376a630ac4e4b2eeebb9
at `artifacts/reference/actual-wave-geometry-source-20261024/proposal-source-only.json`.
Review required two validity fixes before this PASS: exact3001/full0:.02:60 grid
and raw native-body time identity when available; two-sided rim-clearance
roundoff allowance and subtraction from reported curvature lower bound. These
prevent a decimated history being called complete or near-zero floating-point
contact being falsely proved immersed. No geometric fixture or numerical job
was run to obtain source readiness. The floating-point margin remains an
engineering allowance, not interval-arithmetic machine certification.

# Independent first native motion review, 2026-10-24

Reviewer capytaine_expert, read-only output/source review. Consumer authored and
executed the WEC-Sim case; reviewer authored no consumer physical correction and
launched no simulation. This review inspects stored native outputs with algebraic
Newton-Euler and block-timing identities; no independent motion solver/reference
load is substituted for WEC-Sim physics.

## Outcome and scope

PASS as the minimum17-panel native WEC-Sim fixedstep reference result, with the
identified temporal discretization limitation below. No mesh-convergence or
temporal-convergence accuracy claim. Final encoded video audit passes below.

Artifacts under artifacts/reference/minimum-20261024/case-attempt1:
summary.json, preflight.json, motion.mat and native-output.mat. Frozen input/
source manifest is ../case-attempt1-source/manifest.json; maintained fork commit
f962103a348d63101ce669c2611f7d7c109755c3. Actual execution session44817 terminal
exit0/chunk860bab; native model closed and process termination recorded before
lease release. MATLAB R2025b, ode4 dt0.02, one thread, all10s retained.

Reviewed output SHA256:
summary.json2c46f4e54e24f3f8a275dae35251f2a0ddc174e377088b8bad029cea31f2763b;
preflight.jsone43fda1a84865edf0fc9e1b0952209dd14aa459b89b3e7dfad7c803f7102d513;
motion.mat78c9a354e81c87e91a1294db15ea59645f52504bd6e6bfafb33d7e3f1c879b30;
native-output.mat80a16af1bec95a09aef43377e4390c1b5e221e871772bb481dc500b9169d1ca8.

The result has501 native samples on[0,10]s, finite six motions, q0 at most
8.67e-19m, all initial velocities exactly0, initial radiation force exactly0.
Peak absolute motion[m,rad]:
[0.00522222086609,0.00496991746600,0.00593544998690,
0.000333667120062,0.000533874486830,0.000356447704650].
Actual H5 CG/CB/volume/fullK/nativeAinf readback assertions passed. Selected
compiled FileSolid plus residualA reconstructs fullM+A to1.5916157281e-12.
Prepared/applied excitation agrees with independently accepted501x6 retained
history to at most1.72e-12/1.61e-12, including correct signed three-direction sum.
No artificial damping, PTO/mooring, force replacement or startup prerun.

## Apparent initial full-mass imbalance: resolved native numerical mechanism

At initial rest, instantaneous physical reference equation would use M+Ainf.
Stored initial acceleration multiplied by M+Ainf differs from initial netFe;
this discrepancy is not excused by the small incident tail or force magnitude.

Source bodyClass.adjustMassMatrix transfers a symmetric subset of Ainf into the
physical FileSolid tensor. Let S denote that complete adjusted Solid6x6 matrix,
and Ar=M+Ainf-S the remaining, possibly nonsymmetric, acceleration-force matrix.
Selected native HydroBody TransportDelay SID8:4071:6251 delays From2/acc1 by
10e-8s=1e-7s before the B2B->WaveRadiation Acceleration input. Product1 in
system_8_4071_6402 multiplies Ar by that delayed/extrapolated acceleration; its
F_AddedMass diagnostic is subtracted in the native '+-----' total-force sum.
The delay InitialOutput is not explicit, so native default0 applies.

The [official R2025b Transport Delay documentation](https://www.mathworks.com/help/releases/R2025b/simulink/slref/transportdelay.html)
specifies zero initial output by default and forward extrapolation from stored
points when delay is smaller than integration step. It does not simply delay
the signal by one whole step. Native initializeWecSim suppresses the small-delay
warning, but retains that documented numerical behavior.

Independent read-only identities on the stored actual output establish the cause:

- At t0, F_AddedMass=0 and S*a0 matches stored Ftotal0 to1.39e-17. Thus the
  omitted initial residual-acceleration term completely explains the initial
  fullM+A imbalance; CG/frame or products-order mistakes are not needed.
- Across all501samples, physical Solid Newton-Euler balance holds with world
  inertia Iw=R*Srot*R' and torque Iw*alpha+omega cross(Iw*omega).
  Max absolute residuals by axis[N,Nm] are
  [1.7053e-13,1.1369e-13,1.7053e-13,1.6209e-14,2.1316e-14,1.0658e-14].
  This simultaneously verifies the native RzRyRx/worldomega/worldalpha mapping.
- Full native diagnostic force sum
  Ftotal=Fexc-Frad-Fadded-Frestoring-FmorisonViscous-FlinearDamping
  matches exactly in all501samples.
- At the first positive sample, Fadded1=Ar*a0 exactly. For every n>=2,
  Fadded[n]=Ar*(a[n-1]+(1-1e-7/0.02)*(a[n-1]-a[n-2])).
  Maximum absolute deviations by axis are
  [1.1795e-12,1.0090e-12,9.2371e-13,3.6859e-14,4.6629e-15,2.7534e-14].
  This identifies the precise nominal-delay forward-extrapolation rule across
  the entire trajectory, rather than guessing a one-sample delay from t0 alone.

Physical consequence: result is the native WEC-Sim mass-transfer plus explicit
acceleration-feedback discretization at dt0.02, including its initial zero
feedback and subsequent linear extrapolation. The nominal1e-7 delay does not
make finite-step effects automatically negligible. Maximum diagnostic
Fadded-Ar*a_current discrepancy is0.3733N among translations,0.01133Nm among
rotations; these are force diagnostics, not a motion-error bound. Radiation
memory is separately initialized zero as required. No frame/inertia/force-sign
correction is supported by this evidence; no physical code was changed.

## Acceptance limits and cheapest follow-up if needed

The first minimum-case histories are internally consistent with the actual
selected native algorithm and all declared rest/forcing/frame requirements.
Their absolute motion accuracy is not certified by the algebraic checks, energy
claims or small residuals. Both17-panel hydrodynamic discretization and dt0.02
feedback discretization contribute unquantified reference uncertainty.

No new experiment is needed to identify the observed initial mismatch: stored
outputs satisfy the exact source/documented relation to roundoff. If coordinator
requires a temporal-sensitivity bound before relying on finer accuracy, cheapest
decisive experiment is ONE same17panel full10s run at dt0.01, keeping exact
hydrodynamics, packet coefficients and rest state. Compare all six whole native
histories on common timestamps, reporting metres/degrees and startup, without
clipping, warmup, coefficient gates or automatic mesh campaign. Consumer owns
implementation/execution if separately allocated.

Coordinator accepted this limited first minimum fixedstep reference following
independent review. No source correction or temporal campaign is allocated.
Final artifact rendering and independent audit are now complete below. Later
motion variation and absolute accuracy remain outside this minimum-case acceptance.

## Final centered encoded video/frame audit: PASS

Accepted artifact: case-attempt1/render-attempt3/reference-isometric.mp4, SHA256
a79a75cd88fa03c8b8d46d833cd445645e58a906995a826be9e5674ef345664f.
Consumer exact encoder session6890 terminal0, all FFmpeg decodes/ffprobe terminal0,
no encoder/Python child remains, lease released. Runtime physics/motion unchanged.

Reviewer visually inspected literal encoded frames, not pre-encode previews:
render-attempt3/qa-encoded/frame-000.png, frame-125.png, frame-250.png. Actual
0/5/10s states show full centered hull without occlusion or header/telemetry
overlap; fixed orthographic isometric camera, equal actual metres, no motion/
height exaggeration. Prescribed incident-only surface label is explicit.
Native t0 telemetry is zero; t10 telemetry[5.222,4.970,0.396]mm and
[0.0191,-0.0306,-0.0122]deg agrees stored native pose and reviewed CG/RzRyRx
transforms. Small apparent hull motion reflects actual weak baseline response.
Six-motions.png displays correct ordered native histories/units over full0..10s.

qa-encoded/ffprobe.json confirms H2641280x900,251frames at25fps,10.04s duration:
native states cover0..10s inclusive with final10s sample held one encoded frame.
qa-encoded/source-invariance.json confirms all5 frozen renderinputs unchanged,
including exact packet coefficients, native motion, canonical hull and renderer;
the original draft video remains unchanged. Literal decoded frame SHA256:
0=37bcf483d03ad52b7d7ed11aa59cd2d21ecf0603fe13c30095e8bd3bed737a8e;
125=fb0b50d5dc23c471326b156d07b068452f410836ed1d1b068f7769552f15af10;
250=a2c07e35eae5d029fc83fe6f14703e65481567ec3ae9a808e540493897d11ac6.

This completes independent first minimum native motion/video signoff. Temporal
and mesh motion accuracy remain unestablished. Preserve the trusted baseline and
all prior iterations. Latest user requires future larger-roll case from project
sea-state generator; canceled bespoke RAO period/heading optimization was never
implemented or launched by this role. Await next bounded assignment.

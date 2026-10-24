# Consumer interface review, 2026-10-24

Reviewer: capytaine_expert; read-only consumer ownership review. No consumer
implementation changes or numerical calls by reviewer. Applied MATLAB skill
for source review. Producer hydrostatic correction was authored here and was
independently reviewed by startup_review, not self-approved.

## BEMIO mapping source review

Reviewed reference/wecsim/prepare_bemio.m at SHA256
1d42b6c569d63db9206e2b8ebf51abbb649a4b76642ed0ea6f0591afbe13a8c0.

- Official readCAPYTAINE remains importer. Exact native-gl2 CG/CB/volume replaces
  rounded native legacy text. K is normalized by rho*g once, as expected by
  official H5 writer/body preprocessing.
- Native direct Ainf divided by rho overrides radiationIRF's finite-band estimate
  after kernel construction. Official writer stores normalized Ainf, body class
  restores dimensional added mass with rho.
- FK+diffraction excitation compares named native dimensions and canonical DOF
  indices against imported real data and negated imported imaginary data. One
  imaginary sign conversion matches native exp(-i*omega*t) interface requirement.
- Native M is present in consumer MAT workspace but needs explicit case/body
  mass, diagonal inertia and products mapping; H5 hydrodynamic import alone does
  not establish full dry inertia. Case mapping/compiled evidence pending.

Source compatibility body_name alias was added in producer's separately reviewed
native-gl2 output; raw native3 body coordinate retained.

## Executed hydroimport joint signoff: PASS

Official BEMIO execution terminal exit0, exec chunk547243, engine9168 evaluated
synchronously and detached. Log artifacts/reference-bemio.log reports completed
official import, radiation/excitation IRFs and H5 write. Frozen source/input
manifest artifacts/reference-bemio-source/manifest.json matches all current
recorded source/input hashes, including accepted exact native-gl2 output.

Read-only artifact inspection of stored consumer/bemio.mat (no new simulation,
MATLAB session, solver or physical numerical experiment) confirms:

- CG/CB/volume exactly equal accepted native MAT. Density1025kg/m3,
  gravity9.81m/s2, infinite depth, sixDOFs and propagation headings0/45/90deg.
- Stored normalizedK equals nativeK/(rho*g) exactly; normalizedAinf equals
  direct nativeAinf/rho exactly. Saved full dryM equals nativeM exactly.
- Imported A/B/Ainf/Khs/ex_re/ex_im and radiation/excitation IRF arrays finite.
- Executed interface-check.json reports max real and single-imaginary-sign
  residuals1.4551915228366852e-11 each against native FK+diffraction.
- Official H5 output exists,1107552bytes,
  SHA256d54d5f5474b2d71028e5953cb6fb59dcfdeaa950439955b86056ed7e5583033e.
  bemio.mat SHA256b175794970defd2faa7571d05f5077d0ccc2ecb0fc7f6aa8ab15d2190fa2114c;
  interface-check.json SHA256200a13712564913b0cec6c832c6d0b20df67a508a3191823ed03ee63cae8e553.

Acceptance scope: hydrodynamic import and documented official writer units/DOFs/
phase/metadata. Reviewer did not independently reload H5 storage; eventual
bodyClass readback, full physical body-frame/dry-plus-added inertia tensor,
actual-rest forcing/startup and motion remain separate acceptance steps.
No motion uncertainty/accuracy result follows from this signoff. Corrected
hydrostatics' physical heave/pitch coupling remains preserved.

## Products-order defect: independently reviewed reproduction

Selected maintained fork HEAD a23d2993807410fe78cc63546c1ba9d81912024f;
bodyClass.m SHA256697552816c164c34207a7cdb32002d43a139871203864571764ef9b08f0a59bd;
WECSim_Lib_Body_Elements.slx SHA256
0ef996e4ef79f1eef3bcf451719a39e397820937f3e4586a3e99bcf8c61beaeb.

bodyClass.adjustMassMatrix constructs rotational added products from
[fAddedMass(4,5), fAddedMass(4,6), fAddedMass(5,6)]. Library ZIP member
simulink/systems/system_8_4071_6208.xml sends adjustedInertiaProducts directly
to ProductsOfInertia. Official [MathWorks R2025b Inertia block documentation](https://www.mathworks.com/help/releases/R2025b/sm/ref/inertia.html)
specifies tensor off-diagonal entries ordered [Iyz Izx Ixy], retaining the
negated-integral convention already represented in the inertia matrix.
Canonical dry matrix only has nonzero Ixz, the middle entry, so that part alone
cannot expose this order problem. Native added-mass upper products are
[855.7901890552873,46.24564776255444,-215.8532092739416]kgm2, and the opposite
entries differ: raw coarse collocation data remain preserved without symmetry
projection.

Consumer actual-path reproduction is now independently reviewed PASS:
artifacts/reference-inertia-stock-attempt2/result.json,
artifacts/reference-inertia-stock-source-attempt2/manifest.json, and
artifacts/reference-inertia-stock.log.attempt2. Source script SHA256
3e55fd483ba9179e183b9e90e22dd6538d964db7e27a8a6e6ef25371b58999ef;
maintained fork2f43f547301b223cd9a07e84dc88c9375249af0e on official v7.1 base;
bodyClass/library hashes above unchanged. MATLAB selected R2025b engine9168.
Actual selected block MaskType='File Solid', actual ProductsOfInertia expression
evaluated after official H5 read/bodyClass hydroForcePre/adjustMassMatrix.
No compiled body simulation was run. Ledger records synchronous eval/command
terminal exit0, chunk136dd2, library closed without changes.

Expected versus actual checks independently match:

- Canonical dry-only middleIxz: tensor error0, the immune control.
- Distinct dry offdiagonals[11,13,17]: block reads yz=11, zx=13, xy=17, yielding
  xy/yz swap and max error17-11=6kgm2, exactly reported.
- Canonical fullM+Ainf: upperxy855.7901890552873 versusyz-215.8532092739416 swap
  produces max error1071.6433983292288kgm2, exactly reported. The full original
  nonsymmetric residual added-mass matrix is included in reconstructed total;
  physical coupling error is not attributable to rejected coefficient asymmetry.

Independent actual-block reference additionally checked:
[MathWorks R2025b File Solid](https://www.mathworks.com/help/releases/R2025b/sm/ref/filesolid.html)
uses [Iyz Izx Ixy], directly the off-diagonal tensor entries. No new minus sign
belongs in the mapping. Physical consequence: stock selected path applies the
wrong roll/pitch versus pitch/yaw inertial coupling for this actual boat.

Consumer proposes reordering only three actual library ProductsOfInertia
expressions by[3 2 1] using Simulink API, retaining public/internal[Ixy,Ixz,Iyz],
all raw coefficients and asymmetric residuals. That bounded boundary correction
is mathematically appropriate and minimal. Consumer owns implementation and
regression; reviewer will not author the correction it reviews. Source edit and
executed regression are reviewed below; compiled body/frame evidence is pending.

Read-only correction source review PASS: reference/wecsim/fix_inertia_products.m
targets exactly the3identified native FileSolid fields in Flex, Drag and Hydro
body variants via load_system/set_param/save_system. It appends only the boundary
permutation, preserves public/internal semantics, and rejects repeated application.
check_inertia_mapping.m repeats all3immune/distinct/native cases through actual
H5/body adjustment plus actual block expression, requiring each reconstructed
M+A error below roundoff and verifying all3expressions. No consumer code edits
authored by this reviewer.

## Executed inertia correction and R2025b library save review: PASS

Evidence: artifacts/reference-inertia-fixed/result.json,
library-semantic-diff.json, artifacts/reference-inertia-fixed-source/manifest.json,
artifacts/reference-inertia-fixed.log. Exact execution session26461 terminal
exit0/chunkca26f4, library closed; recorded engine9168 R2025b. All3expression
checks passed; tensor errors[0,0,1.8189894035458565e-12]kgm2 (last is roundoff,
not exactly zero). Original nonsymmetric residualA remains unchanged, and
reconstructed totals recover original fullM+A. Correction source hashes:
fix_inertia_products.m8c3c9e1c09a65dc7ec5ef3de0d1dab5b1b4ba0b49a02b423724fdc603a6b479d;
check_inertia_mapping.mf5efbd8e45bc4360e1bda16432b6979bfb4546aecd6f207672025532162a7434.
Original library snapshot retains SHA0ef996e4ef79f1eef3bcf451719a39e397820937f3e4586a3e99bcf8c61beaeb;
corrected saved library SHAe15104ac2310d2e13125fdecdb1b84f4a515f2b2b7812c69afca34b35694f3cd.

Proper save migrated R2020b serialization to R2025b; it did not change only3
serialized fields. Independent read-only ZIP/XML inspection confirms:

- All531 old Ports vectors exactly match new PortCounts attributes when keyed
  by system/block name. Those omissions are serialization migration.
- All146system graphs have identical source/destination connectivity after
  generated SID renumbering is normalized to block names/types. All39embedded
  MATLAB/Stateflow scripts unchanged. No blocks added or removed.
- Exactly3ProductsOfInertia boundary expressions contain intended permutation;
  all other physical parameters unchanged on selected rigid hydrodynamic path.
-117ZOrder changes and8LibraryVersion updates are save metadata.
- RigidBody AllowZeroVariantControls='on' becomes EmptyChoice='on'. Official
  [Simulink R2025b rename documentation](https://www.mathworks.com/help/simulink/slref/simulink.variantutils.converttovariantsubsystem.html)
  explicitly retains previous zero-output/no-active-choice functionality.
- FlexBody MotionSensor reference gains AllowFlexibleInterface='on', the
  [R2025b Variant Subsystem default](https://www.mathworks.com/help/releases/R2025b/simulink/slref/variantsubsystem.html).
  This reference is in the inactive FlexBody variant; selected RigidBody motion
  sensor subsystem parameters and signal connections are unchanged.

Disposition: boundary correction plus R2025b save migration accepted for selected
rigid hydrodynamic case; no extra correction needed. No blanket all-library
behavioral certification. Body compilation/initial frame/outputpose and actual
motion remain separate; renderer must transform mean-world hull about worldCG
exactly once, without reapplying equilibrium orientation.

## First case physical/frame source review: PASS, runtime gates pending

Read-only reviewed reference/wecsim/wecSimInputFile.m, run_reference_case.m,
build_reference_model.m against actual selected native library. Source permits
a first bounded body launch after independent full physical-hull XY startup
diagnostic passes; no actual motion/readback/compiled result claimed here.

Reviewed source SHA256 identities:
wecSimInputFile.m603a0ce666f7fa88d1d65fdf8abcdc65c990ea6240321574133a4df8dfabc48b;
run_reference_case.m9c90d6c5d4ed6a7fa8ac0777cfe28a79241fb64fd85f128aa6b5979142acb439;
build_reference_model.me003aa055a971742f1dd4b0cef4f5cc45da8ec9a8ae448bebc8f9fe7a4a8b555;
native quaternion2EulXYZ.m6d00d27278f9e892f69883b9d11c90b0eee34b7813e0ce831349f4243917e5da;
native Constraints library117a5a1688b84c29664ad551a186401501517170939d500d332d3105a1904da8.

Mean-world producer boat.stl already includes equilibrium orientation. Runner
subtracts meanworldCG exactly once for local FileSolid geometry. Native selected
Structure setupCG rigid transform then adds body.centerGravity once; initial
transform is zero displacement/angle; FileSolid localCOM=[0,0,0]. Thus native
worldM and meanorientation remain consistent, with no second equilibrium rotation.
Runtime asserts meanCG/CB/volume and fullK/directAinf against exact accepted MAT.
DryM uses explicit mass, diagonal inertia and raw publicWEC products with the
separately reviewed boundary permutation.

Native Floating(6DOF) is the standard BushingJoint: all6 motion primitives are
ComputedMotion, actuationNoTorque, springs/damping0 and limitsOff. No PTO/mooring
or empirical damping is added. initialq/angle0, actual native initialvelocity
checked0 after run; fixedode4 dt0.02 over[0,10]. initializeWecSim clears
convolutionIntegralInterp/Surface persistent memory, and convolution history
initializes with zeros. Startup forcing remains the independently accepted
three-direction packet with rampTime0; no force replacement or prerun.

Selected world MotionSensor spans actual CG against WorldFrame, so translation
is nativeCGworld minus meanworldCG. Native quaternion2EulXYZ normalizesQ and
returns roll/pitch/yaw radians in the standard Rz(yaw)*Ry(pitch)*Rx(roll)
representation. Velocity angular entries are WORLD angularvelocity omega;
acceleration angular entries are WORLD angularacceleration, not Euler-angle
derivatives. Directbus columns1:6/7:12/13:18 match stock responseClass extraction.
Raw bus/time retention bypasses presentation postprocessing, not solver forces
or hydrodynamic physics. Postrun assertions require finite6histories, actual
q0/v0rest, complete10s fixednativegrid and nontrivial motion.

Eventual renderer formula for original mean-world hull vertices is
v(t)=CGmean+qtranslation(t)+R(t)*(vmean-CGmean). Equivalent localCG geometry
uses CGmean+qtranslation+R*vlocal. Apply neither extra equilibrium quaternion
nor extraCG shift. Renderer source and compiled/runtime case evidence pending.

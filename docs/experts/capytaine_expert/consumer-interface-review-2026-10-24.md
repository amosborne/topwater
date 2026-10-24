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

## Focused products-order lead (not yet reproduced)

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
Therefore the selected block's precise mapping needs independent consumer
reproduction before reliance; source suggests xy/yz swapped in the added-mass
adjustment. Canonical dry matrix only has nonzero Ixz, the middle entry, so that
part alone cannot expose this order problem. Added-mass xy/yz can be nonzero on
the asymmetric coarse hydrodynamic discretization. This is an exercised-path
lead, not a confirmed numerical defect or authorization for broad library audit.

Consumer and coordinator informed with exact source paths and R2025b primary
reference. Consumer owns smallest decisive reproduction and any confirmed fix;
reviewer will not author the consequential consumer correction it reviews.

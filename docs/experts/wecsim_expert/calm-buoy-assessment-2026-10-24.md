# Explicit REF002 numericalcalm reassessment: same immutable train2

Coordinator declared revised engineering acceptance, informed user and updated
sharedcontract. Original REF001 .1N/Nm caps remain attached to earlier rejected
train1/train2 reports; do not rewrite those results as old-policy PASS.
No sea/phase/window/C/force/H regeneration, rescaling, retiming or history seed.

Current train2 initial heave12.7276N is about.052% of bodyweight. New epsilon=.001
policy uses rotation-invariant ||Fxyz||2<=epsilon*m*g and ||Mxyz||2<=epsilon*m*g*Lchar,
PLUS each |F_i|<=epsilon*S_i in the SAME fixednative meanboat basis as H/F.
Positive carrierS is not a vector to rotate. Lchar=maxmaterialhull L/B=7.3152m,
rotation-invariant geometric property. Bodymass2500/g9.81 gives force24.525N,
moment179.40528Nm caps. These are explicitly stated engineering numericalcalm
assumptions, not guarantees of physical accuracy/equilibrium/linear applicability.

Existing unchanged eta<=.001sum(a) and fluidspeed<=.001sum(a*omega) fullhull-region
criteria still REQUIRED. q/v/radiationmemory0 remain exact native body initial
requirements after eventual field acceptance; no body has run for this new sea.
Arrivalgate/physicalparent/train2 packet/force are unchanged.

## Minimal source and rigorous combined-field bounds

NEW `reference/wecsim/assess_calm_buoy_train.py` reads immutable
`artifacts/reference/calm-buoy-20261024/train-attempt2/{train.npz,body-force.mat,diagnostic.json}`.
All six existing train2 files are hashed before/after; no new train/coefficient
file is written. Force0 crosschecked against original stored report; original
rejections reported openly. Only output is a NEW policy/field assessment.

Full meanworld hull planform bbox x wetdepth:0, originalneutralgeometry unchanged.
Combined ALLsignedmodes are summed as COMPLEX derivatives before abs or vectornorm,
so unrelated parent-component maxima do not manufacture conservative failures.
Taylororder6, initialmaxbox1mx1mx.2m, batch8centres. Forbox halfwidth h,
r_l=|kappa_x|h_x+|kappa_y|h_y+|kappa|h_z; surfaceeta usesxyonly. Remainder
sum|C_l| r_l^7/7!, velocity sum sqrt(2)*Omega_l*|C_l| r_l^7/7!.
Exp integral remainder is valid because completebox upperz<=0, so vertical
growth times center decay is bounded by1. Allfieldterms are exact spectral sums,
not profileinterpolation or regenerated forcing.

Adaptive subdivision is a computational proof refinement ONLY: samefield/order/
gates, up to3levels within480s. An actualcenter/sample violation is a definite
regionviolation; an unresolved conservative upperbound or runtime-limited pending
box is explicitly INCONCLUSIVE, not a claim of a bad sea. Leaf/pendingboxes,
upperbounds, center values and levels saved. No additional accuracythreshold or
physical-parameter search. Assessment can stop with honest partialcoverage flags.

## Exact readiness, resources and stopping point

ProducerSHA1f0e8f32affc5b38cfe83176d57a6460047a37f585d12a92d5eae239f6ac2a1b;
controlsSHA946a257879808b7b1372940786d2294ac0b9be49b3d89032bcf7ac9079167ee1.
Prospective manifest:
`artifacts/reference/calm-buoy-20261024/assessment-attempt1-source/proposal.json`,
SHAb6dd95cc6dfb60e1f69349b5620cf8ef65092b24f3367be3fc6cd57a529eeda9.
Actualassessment/module/physicalinput+sixoldtrainfiles and installed numpy/scipy/
threadpoolctl files frozen separately. No field calculation at readiness stage.
Independent source/math/freeze review requested before any numerical grant.

```
python reference/wecsim/assess_calm_buoy_train.py --output artifacts/reference/calm-buoy-20261024/assessment-attempt1
```

OneCPU/noMATLAB, internal480s/external600s, expected30-180s. Same guard Claim ->
SEPARATEverify -> Run -> exactterminal -> idleUNDERlease -> SEPARATEreadback ->
FINALrelease. Stop after actual assessment; no automatic body/encoder/fieldchange.
If accepted under this declared policy AND field proof, later nativecase/renderer
must reference this new assessment report with exact oldtrain2 C/force identities,
not pretend oldpreflight accepted. Nativebody remains originalH5/CIC10/.02/free6/
rest0/memory0/viscous0/ramp0, with one-time meanboat/globalframe mapping.

## Executed assessment PASS and next body readiness

Independent SOURCE/FREEZE PASS then coordinator allocated exactly abovejob.
Actual20.984s accepted under NEWexplicit engineeringpolicy. Fullhullregion
eta rigorousbound21.033616um <2444.851712um; vectorfluidspeedbound2.290736mm/s
<2.489977mm/s. All24eta/72velocity leafboxes level0 cover fullarea/volume without
overlap/pending/unresolved/actualviolations; no subdivision needed. NormF12.728635N
<24.525N, normM1.162293Nm<179.40528Nm andallfixed-basis epsilonS/arrival pass.
No C/force/physicalparent change; oldtrain2 report remains old-policyREJECT.
Independent actual result/hash/coverage review PASS after terminal/release.

Output `artifacts/reference/calm-buoy-20261024/assessment-attempt1/diagnostic.json`,
SHA367a08f33af6fec4aca2cf410d2d8492a294dcaa7ba2e02d6ccefd595113f1e2;
boxproofSHAe90a8b82d8fb31cfc8d5def595cb8f2f6327c7d6438fb6df226712796441a18b.
Exact20280/30460 terminal0/session13794/chunk1e0731/CIMa562d9 absent;
all12unchanged/idle5c5e64 -> SEPARATE98df90 -> FINALrelease0aee90.

Body minimalSOURCE/FREEZE independentlyPASS0849bcdd29d5aa1c37964406eb6e288c38ac2c871f0bc5c5a0571d9b458d470a,
`artifacts/reference/calm-buoy-20261024/case-attempt1-source/proposal.json`.
Runner2af37237/input9d7d67aa/preparebc0f7bac useoriginalH5/nativeGL2M/K/Ainf/STL,
CIC10/dt.02, q/v/memory0/ramp0/viscous0/free6. Same immutabletrain2 Fe applied
once via reviewedcasehook; ALLsixassessment input hashes crosschecked beforefreeze.
Globalstation bow201.308deg R/t is embeddingmetadataONLY; native M/Htheta/forces/
modelworld unchanged. NativeabsoluteCG andoriginalCG-relativeSTL mapping happen
once, then Rglobal/t once for display/globalcoordinates. No heading doublecount.
Fullraw3001/CIC/finite/Fe/H5mass-inertia/NE/Ar tests retained.

```
python reference/wecsim/run_matlab.py reference/wecsim/run_calm_buoy_case.m --log artifacts/reference/calm-buoy-20261024/case-attempt1.log
```

Oneexistingengine/oneCPU/max600s/expected30-180; no actual body/encoder yet at
this readiness checkpoint. Source146inputs+archives/deps verified, no remaining
sourceblocker. Bodyexecution still requires separate coordinator allocation.

## Actual mild60s native BODY implementation PASS

One body run completed27.776s after independent source/freeze review and exact
coordinator allocation. Fresh output `artifacts/reference/calm-buoy-20261024/case-attempt1`.
Full3001samples0:.02:60 preserved. Knowninitialq<=8.7e-19m/rad, allinitialv and
nativeCICforce exactly0, ramp0/viscous0/free6. ALLsix prepared and actual applied
excitation histories match immutable train2 forces with EXACT0 error. Native
3rowRT/PCHIP/CIC check<=2.274e-12; rotatingSolid NewtonEuler<=4.548e-12N/Nm and
nativeAr extrapolation<=6.336e-11N/Nm. OriginalH5/nativeM/K/Ainf/CG/fullinertia
readbacks pass. These are implementation/consumption checks, not a mesh/time
accuracy or nonlinear validation. Independent stored-body review requested.

Peak absolute native motions: surge.2380058m/sway.2744882m/heave.2402980m;
roll.02585130rad (1.481deg), pitch.01156304rad (.663deg), yaw.01076037rad
(.617deg). Whole60 response bounded, no cropped startup/detrend/retiming or rerun.
Physicalglobalstation embedding is recorded in ACTUALsummary: yaw3.51348990965954,
R=[[-.931639590951216,.363383643787459,0],[-.363383643787459,-.931639590951216,0],
[0,0,1]], t=[0,0,0]. Rawposition_world/q/v/omega remain NATIVEmodelworld/meanboat;
apply globalR/t exactlyonce for globalCG/geometry/vectors/video. Originalequilibrium
trim/CG handled once in nativeSolid; no modelyaw or second CG/phase rotation.

Hashes:

- native-output.mat `6b8a9adc5224d6bbb92541af3c6350771b99ee638af98e7c30d4999d31739e57`
- motion.mat `f822fde31bae1f5fd4e0691d6fd9467143786af6631cc981c591df63cd8daaf4`
- summary.json `aab25bc9a17af2900269555debdc0732baf6a76effe87090af4b835895b539bf`
- native-cic-output.mat `a366294ce6e125361e3868b456e45367734f226627450306920221a520a2b998`

Prelaunch approval review rejected explicit engine.exit in readiness query as
possible sharedservice shutdown. Safe read-only query without exit/quit passed;
root authorized minimal clientcleanup compatibility change from explicit exit()
to normal attached-client reference cleanup. Independentreview clarified destructor
internallycloses attachment; sharedconnect semantics disconnect attachedclient,
not a confirmed enginebug. Old0849source snapshot preserved. Exact finalbodymanifest
`case-attempt1-source-ready2/proposal.json` SHA
5e7f4836173c8bebe0dcf5c4252f090372f2129f81b37425f6f87acb64209349,
client91b5d4ec/bodyphysics2af37237; independentall146hashes/archives PASS then
coordinator granted this exactnewidentity. No unreviewed numerical retry.

Exactclient29068/guard39568 terminal0/session67588/chunk1fa9ae; authoritative
CIM671424 shows both/children absent, onlyknownnonnumericalengine docschild33628.
All146frozeninputs unchanged; idleUNDERleaseec580f -> SEPARATEaf7ffc -> FINAL
releasedd4b47. Noledgerwriteafterrelease. Log warnings(modelname shadow onpath,
disabledlinks for copiedcase diagnostics) retained; localexplicitmodelpath,
compiledreadbacks and actualforces checked. No physical/source correction claimed.
No further body/energy/mesh jobs or encoder by this owner. STOP afterreceipt;
independentstoredreview and separatevideo binding/preview/export owned by root/video.

Independent final MILD60 storedreview nowPASS: all146inputs/output/assessmentlineage,
native3001restIC/visc0/grid, Fe/fullwrench exact0, independentlyreconstructed CIC
and NE/Ar/mass/frame checks match. Reviewer confirms globalR/t/yaw exactexport and
unrotatednativebasis, boundedfinalangles. Ready for root-controlled video binding;
no additional numerical work by bodyowner.

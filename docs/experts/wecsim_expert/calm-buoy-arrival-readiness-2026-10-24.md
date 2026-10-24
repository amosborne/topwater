# Mild buoy-derived incoming train: fixed physical preflight source ready

Latest user regroup supersedes PTO-zero/bug-chasing/parameter-search paths.
PTO-zero source proposal is UNRUN and cancelled; preserve it as history. Target
is a genuinely milder observed buoy sea, initially locally calm water and body at
rest/zero radiation history, physical incoming waves, full60s and centeredvideo.
No temporal elevation/force ramp, amplitude/phase fitting or bodyhistoryseed.

## Selected observation and fixed physical construction

Input owner exported actual NDBC46086 record2023-03-19T08:40Z, Hm0.5351635264m,
Tm027.42692891s,48bands .03-.495Hz, seed42/384carriers. Files in
`artifacts/reference/calm-buoy-20261024/export-attempt1`; componentSHA
814644c993a957528432a6b5aed92d16d1ee94d52119e65e2e7dc26916b7eaa3;
caseSHAb63d82d280b5fd3a70ef5dab93e707e300d1ca7cdbf156efcf66931fc15b26ee.
Default bow201.308143185deg/yaw3.51348990965954rad follows actual generated
dominant arriving energy; no heading fit. Preserve all generated carrier a/omega/
direction/phase arrays. No changing the generator/catalog or source observation.

Root FIXED construction: percarrier propagationcoordinate s, initial complex
profile a_j exp(i*k_j*s+i*phi_j) W(s/cg_j+30), cg_j=g/(2*omega_j).
W=1 for|tau|<=15s;0 for|tau|>=25s; transition1-10q^3+15q^4-6q^5 with
q=(|tau|-15)/10. This is a SPATIAL INITIAL train, never W(t)*globalwave/force.
Nominal upstream front is5cg_j, center30cg_j, rear55cg_j; exact subsequent
dispersion changes that profile. Waves approach from each observed carrier's
upstream direction. Wider physical extents are independent of the close camera.

G(nu)=(1/2pi) integral W(tau)exp(-i*nu*tau)dtau, even real; constantplateau term
analytic, polynomialtransition integrated with fixed256Gauss nodes. Signed
kappa=+/-Omega^2/g is retained. nu=cg_j*(kappa-k_j), coefficients
C_j,+/-,n = a_j G(nu) exp(i*(phi_j+30nu))*(2cg_j*Omega/g)*trapezoid_dOmega.
Common Omega grid.1:.005:physicalcap with exactlast endpoint; cap is min(.6Hz,
deepwater wavelength>=2beam),9ftbeam. Omega<.1 excluded for original native
coverage; abovephysicalcap excluded explicitly. NO coefficient renormalization.
Negativekappa sidebands propagate at theta_j+pi and are retained/reported, not
silently projected away. Same coefficients define surface, fluidvelocity and Hforce.

Finite spatial trains necessarily redistribute the stationary spectrum/MEM.
Report actual60s meanpoint periodogram/variance/bandpower versus observed target,
directional spectralnorm moments, reverse spectral fraction and removed support.
Component Parseval norms percrestlength are not global2D fluidenergy or unchanged
stationaryvariance. No claiming the derived train has the original stationaryHm0.
Underlying unpruned spatialwindow is locallyzero beforearrival; supportpruning
introduces tails that MUST satisfy the actual numericalcalm gates.

## H coverage and cheap-first physical acceptance

Cap exported required excitationONLY, SAME fixed17body, no radiation/H5/mass
change: `excitation-attempt1/excitation.npz` SHA
31c0b956f81f8819c01dad7dcf6b6d42f3526881ff6569d41671b3be70894b8f.
34omega .1:.1:3.4 and24periodic headings plus exact2pi seam, full6 FK+diff,
nativeexp(-iwt)/no secondconjugation/no extraCGphase. Consumer real/imag bilinear
interpolation modulo2pi, boundedOmega, with includedseam (not appendedagain).
Interpolation/mesh accuracy UNQUALIFIED; originalradiationH5/CIC10 unchanged.

Order is fixed: construct/save shared train/support/actual60PSD, compute ALLsix
initial forces and carrier scales. Require |F0_i|<=min(.1N/Nm,.001S_i), explicit
degenerate axes. If forcefails, save fullscreen/report and STOP before costly
fieldbounds. No automatic retiming/front change/search.

Existing25%arrival intent uses combined R(t)=sqrt(sum_j|z_j(CG,t)|^2/sum_j a_j^2),
where z_j includes both signedbranches/quadrature and parentamplitude. Require
max0..60 R>=.25 and report percarrierpeaks. This is an arrival-envelope diagnostic,
not stationaryvariance/Hm0 proof; phases are not fitted.

If F0 passes: use fullhull meanworldplanform circumscribed projectionradius and
wetdepth:0. Evaluate each component on shared s-interval/z-rectangle, spacing.05m;
bound maxima by localgradient radii+global spectralHessian remainder, sum component
bounds by triangleinequality. Four shared spectral matrix products perz exploit
u_h,s=-u_z,z and u_h,z=u_z,s; no587kpoint3D brute-force spectrum evaluation.
Require eta_bound<=.001sum(a), speed_bound<=.001sum(a*omega). All values/depth/
projection/support and conservative bound construction recorded. No fieldclipping.

## Fixed source/API, freeze and prospective job

`reference/wecsim/calm_buoy_train.py` SHA
226b6ca75dfc0d6ecf6f86abb1931f3b27f939e33b24cec4371620e119f18c82;
`prepare_calm_buoy_train.py` SHA
42ae8147cfc6617aa3931251b8efda2935998f9a81480a30c05e54f304c89ef6;
controls SHAb235ae43571005889c06ace0a58833223d061ba87de12e5cc6184ea9b9551d0c.
Ready archive `artifacts/reference/calm-buoy-20261024/train-attempt1-source-ready2/proposal.json`,
SHA01e96e720d32e35f7765b19bdb2079d2cc07b3be1f99d1f9ce645b1939502e7f.
11actualinputs and numpy/scipy/threadpoolctl implementation identity frozen;
source syntax checks pass. Independent source/freeze review requested, no field
numerical job granted/executed at this source checkpoint.

Independent final ready2 PHYSICALSOURCE/FREEZE PASS received at exact01e96e72
proposal: all11current+archivedinputs/installeddependency files match. Reviewer
accepted signedkappa/Jacobian/window/velocity/Hseam/frame, honestfinitePSD/norm
labels, combinedarrival gate and efficientbounds. Existing Topcarrier/Houtputs
independently checked against declared dependency readbacks/FK/hash evidence.
No extraaccuracygate/campaign required. This remains source readiness, not field
execution or body acceptance.

```
python reference/wecsim/prepare_calm_buoy_train.py --output artifacts/reference/calm-buoy-20261024/train-attempt1
```

OneCPU/noMATLAB, internal480s/external600s, expected15-180s; cheapforcefail avoids
allfieldbound work. Guard Claim -> SEPARATEverify -> Run -> exactterminal ->
idleUNDERlease -> SEPARATEreadback -> FINALrelease. No automatic body/encoder.

Shared `CalmBuoyTrain(train.npz)` API: C[384,2,Nomega] is complex and ALREADY
quadratureweighted; omega_rad_s/k_abs_m_inv shared, directions_body_rad/carriers/
IDs/Rglobal/t retained. Field Re sum C exp(i*sign*k*n(theta).x-i*Omega*t).
`profiles(component_ids,s_grid,times)` returns float64[batch,time,s] using CZT
for the uniformOmega portion plus exactlastendpoint. Cache in smallbatches;
spatial interpolation is VIDEOONLY, report approximatecacheerror without inventing
an extra physical acceptance gate. Forcing/fullcalm bounds use rawcoefficients.
Video role owns only NEW renderer/controls/note, no packet/source/model overlap.

After separately accepted preflight/rootgrant, minimal native consumer will retain
originalM/K/Ainf/STL/H5/CIC10/.02, all6free, q/v/memory0, ramp0/viscous0 and actual
sixforce table once. Sameexistingraw/frame/NE/Ar/CICchecks required. Body and encoder
are not authorized by this source receipt; full60s + centered equalmetrevideo
remain the goal. No density/energy-bug study or root-cause assertion.

## Actual ready2 preflight: initialforce REJECTED

One allocatedjob completed2.343s. FixedC/PSD/support/arrival are preserved at
`artifacts/reference/calm-buoy-20261024/train-attempt1`. Initialforce[N,Nm]
[-.1195462,-.1187317,194.5438296,.04918827,-15.6304824,-.1837571] failsfixed.1
absolute caps in fiveaxes. Relative.001carrier scales alone wouldpass, but both
criteria are mandatory. Combinedarrival1.086208 at30.88s passes. Fullhull eta/u
bounds were SKIPPED exactly as force-first policy requires; no watercalmclaim,
bodyrun or encoder. Removedcomponentnormfraction.00249413/reversefraction.00947017;
actual60point centeredvariance.00660023m2 versus observed.0179m2. No renormalization.
ReportSHA58317383106f3cfadd4b5d5811b7d2d122c9b994b559fdc59cb94f8bf95dffe9;
trainSHA712fe11cea4ae69deb1e82096bc6335c7740dc9419a18034a8b989b534b91359.
Exactproducer36188/guard35000 terminal2/tool329b67exit1; authoritativeCIM24175a
absent. All11frozeninputs unchanged/idleUNDERlease802788 -> SEPARATE6cd78d ->
FINALrelease340535. No furtherjob/retiming/search authorized by that receipt.

Root then requests INPUT-only lowfrequency support reasoning, not solverbugs.
Bounded storedC reduction: retained eta(meanXY,0)=.001424052277m while unpruned
initialW(s/cg+30) is exactlyzero there. K[:,Heave]*eta0 gives180.776910Nheave/
-14.416218Nm pitch versus actual194.543830/-15.630482. Remainingdifference
13.766920N/-1.214265Nm is not proof of physicalcause or promisedbridgePASS.
Hard-support cancellationloss is material; lowcut/highcut/nonlocalpressure
contributions have not all been separated causally.

Independent researcher analyticreview accepts H0=K[:,Heave] real/headingindependent:
uniformwaterlevel incidentpressure rho*g integrates -p[n;r cross n] overwetbody,
fixed-body diffractionpressure vanishes asOmega->0. Originaltriangle constant
pressure equals exactKcolumn within2.37e-11. This is excitationONLY, no radiation
orArDC correction. Proposed lowOmega0.. .1 boundedcomplexlinear bridge to existing
H(.1), samecarrier/window/position/no renorm; Jacobian2cgOmega/g makes C(Omega0)=0,
so no duplicatedDCstep. Bridge remains a lowfrequency approximation, not BEMdata.
Implementation and newjob await coordinatorgrant; currentphysicalattempt immutable.

## Attempt2 minimal lowfrequency INPUT correction: source ready, UNRUN

Coordinator authorized source only, preserving failedattempt1. NewfullOmega grid
starts0, same.005step/upperphysicalcap, regenerated trapezoidweights (old.1 node
is now interior). H0 computed from actual wetpanel constantpressure/moments about
meanCG, compared to native exactK[:,Heave]. Independent researcher refines bridge
to H0+(Omega/.1)^2*(H(.1,theta)-H0), consistent with leadingO(k)=O(Omega^2)
fixed-body longwave correction under bounded zero-frequency response. This is
declared excitation approximation ONLY; no native BEM/radiation/H5 change.
C(Omega0)=0 from Jacobian, asserted, never overwritten/renormalized. NoDCstep or
double +/- branch contribution. Window30/15/25, generateda/phase/seed/directions,
upperphysicalcap and gates remain unchanged. Force-first again; no promisedPASS.

Newproposal `artifacts/reference/calm-buoy-20261024/train-attempt2-source/proposal.json`,
SHAfcea8dd2e01992eed3f0bc9a17dd82ebff0aacb764c311ba019028fbe7879730;
moduleSHA60de70b863b5f94ffc711b96374e93eb579ef412e22894eb5c28406ffbb5b6a0;
producerSHA76b14d1af2b0bb071e3d539c7b3eb99a129d6017a7fe40e6d39279860d5215be;
controlsSHA4e1dc021f4858281ad60b1f8890e3d3bd3244608e6a889acc8032939f7e27461.
Actualinput/dependency snapshots preserved separately; Python syntax checks pass.
Independent finalsource/freeze review requested, no numerical allocation yet.
Futureproducer command (only after separate grant):

```
python reference/wecsim/prepare_calm_buoy_train.py --output artifacts/reference/calm-buoy-20261024/train-attempt2
```

Same oneCPU/noMATLAB/internal480s/external600s/expected15-180. No blindupstream
retiming, phasefit, forceinitialoverwrite, source-rescale, bodyjob orencoder.

Attempt2 subsequently independentlySOURCE/FREEZEPASS and coordinatorallocated.
Actual3.062s preflight REJECT: F0[N,Nm]=[-.12332777,-.10914483,12.72756912,
.03484314,-1.14673347,-.18631368], fixed.1caps still fail in fiveaxes. Restoring
lowOmega reducedheave194.54 to12.73N, not a guaranteed or completecalmfix.
H0pressure/Kheave maxdifference2.3647e-11; combinedarrival1.081971 at30.7sPASS.
Fullhullbounds againSKIPPED(force-first); no surfacecalm/body/video acceptance.
CompleteC/PSD/support/arrival preserved; no parameter changes after result.
ReportSHAeefe021e87af066be0bd3ffcca622ed5ff3c88dc3df779641e7a8d0d35d08cbb;
trainSHA4942d9a77cb9a98ae5a2414ac624d3a9c63a8c21a9348d61e4a423217362b3a6.
Exactproducer35804/guard35372 terminal2/tool34f7ecexit1, authoritativeCIM1ded6f
absent; all11frozen/idleUNDERlease2789d0 -> SEPARATEfa2bde -> FINALrelease88686b.
No furtherallocation/retiming/rerun/body/encoder; coordinator nextdecision required.

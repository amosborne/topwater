# One mild incoming packet: rejected incident preflight

Coordinator assignment2026-10-24: prepare a simpler physical incident diagnostic
AFTER the two official supplied baselines passed. No body or incident computation
has run here. No new BEM/export, final buoy-sea design, pruning choice, model repair,
viscosity/damping fit, seed/history initialization or encoder is authorized.

## Frozen prospective diagnostic

One nonzero packet: amplitude.015m, carrier.25Hz, phase0, propagation45deg CCW
from+x, cos^4 compact spectral weight with halfwidth.5rad/s, focus40s at origin.
These are prospective diagnostic values, not fitted to a desired boat amplitude.
Actual support (.25 +/- .5/(2pi))Hz, approximately.1704-.3296Hz, corresponds to
deepwater wavelengths14.4-53.8m, within the existing project lambda>=2beam envelope.
The centroid-frequency group packet center at0 is approximately125m upstream;
exact field is the full dispersive spectral sum, not a translated rigid envelope.

Reuse unchanged reviewed `IncidentPacket` class and derivative/Hessian calm bound
from REF001. That evaluator has three slots; amplitudes[0,.015,0] make precisely
ONE physical packet. Original canonical JSON, original controls/packet/results
and all previous failed studies remain unchanged. No source bank is modified.
Same spectral coefficients generate the field, initial fluidvelocity and forcing;
future visualization must load this packet archive rather than invent a wavefront.
This is not a collected-buoy final case.

Linear finite-band packets have noncompact tails; literal exact flat water ahead
of a nonzero bandlimited packet is unavailable. Original accepted NUMERICAL calm
policy is retained, without clipping/overwriting: full hull planform x wetdepth:0
grid.025m with conservative local-gradient/global-Hessian remainder. At0 require
eta<=1e-3 sum(a)=15um and fluidspeed<=1e-3 sum(a*omega)=23.562um/s. Every signed
direct-native-H and retained-IRF F0 must satisfy absolute <=.1N/Nm AND <=.001S_i;
zero S_i means exactlyzero force, no denominator floor. Report all values/gates,
initial hullregion and spectralband. At least25% amplitude arrival within60s.
No success asserted until this incident-only preflight actually passes.

Forcing origin/phase convention exactly follows REF001: nativeexp(-iwt), total
FK+diffraction H at45deg, interpolated in existing native frequency coverage,
with no extraCGphase. Actual applied force uses original BEMIO+/-20s excitationIRF,
convolution retaining input[-20,80] to produce output[0,60],3001samples at.02s.
Direct spectralH is a separate tail/discrepancy check, not the body force source.
This temporal padding supports future and past lags without body prerun.

## Native body proposal (conditional on independently accepted preflight)

Full60s retained, q0=0 relative to known neutral meanCG/orientation, allWORLDv/omega0,
native CIC memory0 and no previous body motion. No inherited stationary history
adapter/RT initialization/joint velocity target. Native clear/reset path and all6
free joints retained. Ramp0; linear/quadratic/Morison viscous paths explicitlyzero.
Same original17panel fullM/K/Ainf/Ar/CIC40 H5 as the earlier diagnostic, NOT a newly
qualified operator. Default nominal1e-7 acceleration-feedback delay remains native.
Negative native radiation-power witnesses and mesh sensitivity are unresolved.
Bounded motion would not validate those coefficients or prove passive dynamics.

Copied-case builder and reviewed properAPI force hook are reused without fork
editing. A scalar ZERO elevation file is native initialization compatibilityONLY;
it does not represent the physical field or applied force. The one actual positive
imported branch is replaced by the preflight's six-column retained force history.
Full3001 x 6 prepared/applied match and zero initial body/CIC readbacks required.
Raw60-column native output saved before assertions. Preserve absolute CG, WORLD
omega (not Euler rates), native RzRyRx orientation, trajectories and all startup.
No crop/detrend/rebase/hidden force cancellation. Additional direct indexed
zero-history convolution check uses actualWORLDv and hf.irkb; endpoint halfweights
match native trapezoid. Sink-only properAPI logs actual pairedRateTransition times,
threeforcehistory rows and consumed nativePCHIP output. The unchanged reviewed
`verify_seeded_cic_output` is reused ONLY as a postprocessor with identicallyZERO
pre0velocity array; it compares directindexedconvolution rows->PCHIP->actualbus.
No seeded initializer/history adapter is called. Actualviscous/linear bus columns
must be exactlyzero; fullnativewrench sum checked separately.
RotatingSolid NewtonEuler and nativeAr extrapolation checks
are reused from the independently reviewed calm-control runner.

## Prospective commands/resources; no execution yet

Sources/controls: `reference/wecsim/mild-packet-controls.json`,
`prepare_mild_packet.py`, `prepare_mild_packet_manifest.py`,
`prepare_mild_packet_inputs.m`, `wecSimInputFile_mild_packet.m`,
`run_mild_packet_case.m`, `verify_mild_packet_cic.m`. New MATLAB sources passm_check.
Existing producer/consumer/builder/evaluator/fork sources remain untouched.

First archive prospective incident identity (no computation):
`python reference/wecsim/prepare_mild_packet_manifest.py --stage incident`.
Then ONLY AFTER independent source review + coordinator numerical allocation,
claim/SEPARATEverify/run its proposal under the existing guard. Exact producer:

```
python reference/wecsim/prepare_mild_packet.py --output artifacts/reference/mild-packet-20261024/incident-attempt1
```

OneCPU, noMATLAB, internal480s/external600s budget, expected120-480s dominated by
fullregion bound. Save diagnostic.json/body-force.mat/wave-physical.txt/packet.npz
even if calm gates reject (no automatic retiming or parameter search).

Body archive command requires actual accepted incident report first:
`python reference/wecsim/prepare_mild_packet_manifest.py --stage body`.
Its separately allocated native producer is:

```
python reference/wecsim/run_matlab.py reference/wecsim/run_mild_packet_case.m --log artifacts/reference/mild-packet-20261024/case-attempt1.log
```

OneCPU/one existing MATLABengine, maximum600s, expected30-180s. Exactfreshoutput
`artifacts/reference/mild-packet-20261024/case-attempt1`. No automatic body grant
from source/preflight PASS. Each job freezes actual inputs/dependencies, uses
Claim->SEPARATEverify->Run, exactterminal->idleUNDERlease->SEPARATEreadback->FINALrelease.
No new numerical launch is authorized by this source note. Stop source-ready.

Incident actualdependency freeze prepared:
`artifacts/reference/mild-packet-20261024/incident-attempt1-source-ready2/proposal.json`,
SHAad874263acaa5307319762ad54ca7e5f86cc9fa1aff5204653eddb288af800b7.
Only8actualincidentfiles frozen, plus installednumpy/scipy/xarray/threadpoolctl
implementation identities; bodydraft/fork/generatedH5 excluded from that stage.
Producer SHA99d96750e3091c072a710c7f792be797e16e2d1c933a6996e537dc0344bdfab1;
controls SHA9e99469398a0ed57374ba81ba0db1e8b0e7d97abaaa471a898669754449457b8.
First manifest-preparation partialarchive retained: nonexistentoptionalbackend
metadata failed BEFORE computation; current freeze uses installed netCDF3/scipy
path, no package installation or incident diagnostic occurred.
Independent researcher final SOURCE/FREEZE PASS: all8incidentfiles and listed
installeddependencies match ad874263 proposal. BodyRT/PCHIP/zero-past correction
separately SOURCEPASS runner2c14eb43, verifierfac53889, logger06dc6a2b; no bodyfreeze
or numerical body acceptance yet. Coordinator then allocated INCIDENTONLY at this
exactfreeze. Actualsession65212/client32260/guard30676 running under seriallease;
no MATLAB used. No body execution grant follows automatically from incidentPASS.

## Actual incident-only result: REJECTED, preserved

ONE allocated diagnostic completed309.922s. Fullplanform x wetdepth:0 has587,412
gridpoints plus conservative remainder. eta bound2.068649um <15um andfluid speed
4.113119um/s <23.561945um/s pass. Arrival1.0 and support.170422528-.329577472Hz
pass. Initial retained heave+.129757307N and spectral/interpolated-H heave+.147764405N
both fail fixedabsolute.1N cap; allotherfiveaxes pass. The relative.001S_heave
bound alone would pass (S=1246.955726N), but BOTH criteria are required.
No body simulated, cap relaxation, retiming, source correction or automatic retry.
Full60s direct/retained discrepancy preserved (heave max17.604N), not an accuracy
proof or fitted correction.

Output `artifacts/reference/mild-packet-20261024/incident-attempt1` preserves
physicalfield/wave, diagnostic and preparedbodyforce. DiagnosticSHA
ff4d8e4628c3835f9b0e5883a60f9be8a0b60fd0e549c05af16a646dcb7924b2;
packetSHAe4e1f3d00437fe1f244040a3009ccbe0108f6bb12ab7ebfcebc452eec27e25fd.
Exactproducer32260 exit2/guard30676/toolsession65212 exit1/chunk37dde0;
authoritativeCIMd12006 confirms both and children absent. All8frozeninputs match,
idleUNDERlease343fb9->SEPARATE2ec224->FINALrelease202d10. No ledgerwriteafterrelease.
Bodydraft finalSOURCEPASS runnerb8cfa004950132ee78b95986fae0d4b5cfb210db5e6a9539a5dcaf15f049accd,
but no bodymanifest/execution because incidentacceptance failed.

After terminal/release, coordinator requested manifest pathkey correction:
normalize reporthash keys with backslash->slash before looking up controls.
This is archive/readiness parsing only; actual rejected packet/source snapshot
unchanged. Current helper differs from executed archived preparer by that fix.

## Preferred next source-only alternative: extend accepted REF001 unchanged

Coordinator offered original threepacket reuse; recommend it instead of another
singlepacket position iteration. Freeze original amplitudes[.02,.01,.005]m,
carrierfrequencies[.1,.25,.5]Hz, phases[0,pi/2,pi], headings[0,45,90]deg,
halfbandwidths[.3,1.47,.18]rad/s andfoci[100,11.5,200]s. No scale/retiming/search.
The original physical spectralcoefficients already have independent accepted
fullplanform eta3.323565um/speed25.495529um/s andall6retained/directF0 gates.
Keeping their exactfield identity inherits these initial bounds; new focus or
singlepacket revision would require proving them again.

Minimal extension proposal: load frozen acceptedpacket.npz coefficients; assert
arrayexactw/k/c/directions against the original shared evaluator/canonicalconfig;
preserve ALL original[-20,30] inputrows verbatim andappend onlyt>30 through80 from
that same physicalfield. Evaluate originalretained+/-20IRF toproduce0:.02:60 force.
Independently compare its first501 x 6 forceprefix against immutable old expected
retained-force artifact and actual10scase. F0 inheritance then checked directly,
not inferred from generic bandlimited/localcalm claims. No new587k fieldbound
needed if EXACTfield/region/boundidentity is proved; originalbounds stayattached.
All futurewaves/geometry must use originalthreecomponentfield, not onefakefront.

Prefer originalREF001H5/nativeGL2M/K/Ainf/CIC10 operator for a pure duration
extension, rather than silently substituting the latergeneratedCIC40/H5. All6free,
q/v/memory0, viscous0, ramp0, full60startup. First10trajectory is a diagnostic
repeatability prefix only; both originaloperator and60strajectory remain mesh/time
accuracy-unqualified. Existing negative radiation witness not repaired by reuse.
The selected alternative now has incident-extension source/archive ready for
independent review: [exact REF001 extension](reference-extension-readiness-2026-10-24.md).
No numerical extension/body/encoder has run or been granted at this source checkpoint.

# Shared Matplotlib topdown and exact native telemetry: completed

Updated2026-10-24. ROOT accepted the peakposter and linked the result toUSER.
[Five-second topdown clip](../../../artifacts/reference/stronger-buoy-20261024/matplotlib-clip-attempt1.mp4)
SHA954343c52b6bc5690fcf421bc0ac6804c2a4643b2571c26ffc34682f44e87eed.
H264800x800/63frames12.5fps/5.04s; savedtimes50.00..54.96, peakposter51.44s.
17selectedtestsPASS; full63frame decode and ROOTliteralposterPASS. No further
WECcode/jobs; executionowner owns finalterminal/lease+uvlock reconciliation.

## Exact data and telemetry proof

Originalscene1079ae90ec78ea3171204be2996a101caf86015e44771b509b9af82352a3ef4f
is preserved. New scene-telemetry-attempt1.npz,
SHA2475d3befca4483a9ea4b0cacfc06e05586e7cea05e7b7be09ef5dd8dd9e9b23,
adds only optionaltelemetry roll_deg/pitch_deg/heave_m and metadata. All5 original
time/hullpoints/hullfaces/waterpoints/waterfaces arrays retain identicaldtype,
shape and rawbyteSHA. No geometryresample/waveprep/headings/scale/physics change.

AcceptedmotionSHA3eae21f9efa9987d77225b0ebe3f2cac4babec0a482ea01bbbf7801370e6b072;
matchEXACTnative times at zero-basedrows2500:4:2748. Roll/pitch use q[:,3:5]
radians converted todegrees; heave=q[:,2] knownneutralCGdelta, neverfirstclip-frame
rezero. Genericrenderer onlydisplays suppliedtelemetry, no WECpose inference.
ReferenceQA augmentation.json and verification.json contain perarrayidentity and
native-row proof. Missing/nonfinite requiredtelemetry is an expliciterror.

## Shared reusable workflow and timing

WECadapter -> producer-neutral SceneSequence -> sharedTopwater publicCLI:

```powershell
uv run topwater render examples/scene.npz videos/topdown.mp4
```

INPUT andOUTPUT are BOTHrequiredpositionals, relativepaths preferred. Render the
WHOLEprepared NPZ at its savedcadence, fixed800square, onlynamedMP4defaultoutput.
No --extra/start/end/size/fps/output flags. The adapter selectsinterval upstream.
Recordedisolatedenv may use uv run --no-sync withUV_PROJECT_ENVIRONMENT set;
no synchronization/network inside renderjobs. Matplotlib/encoder are coredeps,
optionalPyVista3D remains separate; no GPUcontext needed for default2D.

The standardMatplotlibAgg path uses blitting/directFFmpeg streaming, continuous
blue-white-red water centred at0, lightgrayhullmesh/blackrim/bluegeometrywaterline
andTime/Roll/Pitch/Heave HUD. Mainrender11.25s (setup.765/validation.094/frameavg
.120s); whole referencejob17.797s includingtests/augmentation/QA. This is a
cached-data rerender with no wavepreparation. InitialwideDISPLAYsampling100.75s
belongs to the earlier preparedscene and is reused, not repeated.

Rawgeometry/incidentwater remains the prior provisionallinearphysicalreference;
a 2Dwaterline/color display is not flooding/scattered-fluid/nonlinear accuracy
qualification. No bespokeWECrenderer; genericdrawing/API/CLI/tests belong to
visualization, dataaugmentation/reference QA to WECowner.

## Execution receipt and preservation

Exactmanifestcc6ba729f6104f635f453cd58dd70a44fc9fc26ae9bc8f3fd1fa0ae890a2b589,
matplotlib-clip-attempt1-source/proposal.json; owner/root/visualization,max90s.
SourcecompilePASS, isolatededitablecoremetadata refreshed6a047c beforefreeze127159.
QA matplotlib-qa-attempt1/verification.json and poster.png understrongerroot.
Ownerfinalizes exactprocess/child/inputs/idle/separate/readback/release inshared
ledger; WECowner doesnotwrite anotherrole's executionstate. No morejobs.
All accepted3Dmovies/source/cache and rejectedattempts remain in place.
Any conditional3Dcleanup requires laterUSERreview/authorization; none performed.

## Subsequent cached style/profile and font-only readiness

Cachedstyle/profile attempt2 ROOTliteralPASS: opaqueoffwhitehull/mutedpalette/
fixedmonospacecolumns. ExternalCLIwall12.205s vs libraryrender10.625s; geometry
3.326s/Aggdraw4.283s/pipewait.110s/encoderflush2.187s. Pipewait is blockingwall,
not encoderCPUtime. No encoderbenchmark/tuning is authorized by these metrics.
Originalscene2475d3be/cadence/physics unchanged.

USERnextrequests largerheader only. Explicitownerrelease2abc98 before any WEC
sourcechange. Genericfontfontsize14.5pt is the onlyvisual change; phasecounters
unchanged. Newcontroller skips focusedtestrepeat forcosmetic-only, uses fresh
attempt3.mp4/QA3. Proposal9a77d47fe82c07b54763bd51b6958854c7187eb2290a5ca80075b1349f1b4e87
matplotlib-clip-attempt3-source/proposal.json,90s owner/root/visualization.
Sourcecompile/freezePASS, no dependencies/augmentation/waves/test repeat; no job
byWECowner. Actualgrant/process/ROOTliteralacceptance belongs tosharedhandoff.

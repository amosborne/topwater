# Full60s shared2D reference: complete

ROOT accepted0/peak and executionowner accepted0/30/51.44/60decodedstills.
USER receivedthe fullvideo link. All751frames decode,H264800square/12.5fps/
60.08s (native0..60 inclusive, finalframe heldoneinterval).

[Watch full60s](../../../artifacts/reference/stronger-buoy-20261024/matplotlib-full60-attempt1.mp4)
SHA15510982f6cef3a246f1129222e0525052215a06cfa2e0073af3c53bf803f23b.
QA `matplotlib-full60-qa-attempt1/verification.json` plusframe0000/0375/0643/0750
understrongerartifactroot. Producer/CLI/frozenidentities passed, executionowner
/root/visualization exactterminal/lease FINALreleasebacd1e. No furtherjob allocated.

## Data preservation and minimalborder

Original accepted751hullposes/time and central65x65waterXYZ retained BITEXACT.
Shared fixedview formula allhullbounds/framing1.45 yields center[-.780124,.019896],
half6.100702m. Originalmovingwaterpatch+/-6m missed fixedcorners. Minimumoriginal
.1875m spacing plusonecellmargin yields75x77nodes, halfwidthx6.9375/y7.125m.
ONLY1550newborder nodes/frame were evaluated fromthe SAMEfrozenincidenttrain;
no96m field, newsea/forcing/headings/rescale/IC/physics change. Allfixedviewcorners
covered. NativeHUD matches acceptedq atrows0:4:3000 exactly, not rezeroed.

Reusablepreparedinput `scene-full60-attempt1.npz`,
SHA884667cb0a43f55ec86303c38fae303a76c58a8869c547b1d3290f70b4f330bc,
~104MBwatergeometry. Futurestyle-only changes reuse thiscache, no waveprep.
Fullsceneproducer `prepare_full_reference_scene.py`, wrapper `run_full_reference.py`,
controller `full-reference-controls.json`; explicitreferenceQA separatesdiagnostics
frompublicCLI. Frozenproposal f03211a6f1fd83d7ce34b76f753cfc5f26a159e19c31ef0a9284e05b1072febe.

## Actual timings and interface

Border/data prep40.657s; externalCLIwall148.792s (includesprocess/interpreterload),
libraryrender145.75s; whole referencejob197.911s includingQA. No benchmark or
physics rerun. Current2D Matplotlibblit/directencoderstream, same acceptedlarge
monospaceHUD/opaqueoffwhitehull/graymesh/blackrim/darkbluewaterline/mutedBWRmean0.

```powershell
uv run topwater render examples/scene.npz videos/full60.mp4
```

Bothpaths REQUIREDPOSITIONALS, wholeinput/fixed800/savedcadence/MP4only. No--extra
or crop/size/output flags. Forrecordedpreparedinput use itsrelativeartifactpath
and a FRESHoutputfilename. StandardCLI doesnotknow WECphysics or hull dimensions.

## Approved cleanup and limits

USER explicitly authorized removal of active3D rendering code. WECremoved20
obsolete3Ddisplay scripts/configs/manifests after allreference/library/tests/
physicalproducer importaudit; dependencies were onlyinside thatoldrender cluster.
Newfullproducer and active2D import no deletedhelper. Genericowner removedPyVista
backend/3Ddependencies. NOphysical3Dhydrodynamicgeometry/solver/WECnativefork
visualization was removed; calm_buoy_train and incidentphysics are retained.
Historicalmovies, failedstudies and frozen sourcearchives remain preserved.

The underlying WEC-Sim trajectory is the same provisional fixedmean linear
reference: implementation/forces/startup passed, but mesh/time/nonlinear/curved
waterwetting/flooding accuracy remains unqualified. Displaywater is prescribed
incident field; visualwaterline is not a floodingphysics result.

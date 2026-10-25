# Stronger presentation movie: complete

Updated2026-10-24. USER requested clearer existing motion. No physical input,
body trajectory, force, wave coefficients or accepted751-frame cache changed.
ROOT accepted all4 actualdecodedframes and linked the final movie to the user.
No numerical or rendering job remains allocated.

Final movie: `artifacts/reference/stronger-buoy-20261024/presentation-encode-attempt2/reference-presentation.mp4`.
SHA7a9debb31a440b6967fd3340f2fcd6f89072699dd7ecc2e36136d3bca6079aa8.
H2641280x900/751frames/12.5fps/60.08s. Native0..60 inclusive; endpoint held oneframe.
Full751-frame decodePASS; actualencoded0/30/51.44peak/60stills literally PASS
owner and ROOT. QA `qa-encoded/verification.json` under sameoutput,
SHAdf7a5f7cd944b1d795983d0edc43fe1f8ef9120eb7ba4b24d379c92112648e9e.
Mean absolute RGB differences from original PNGs1.105/1.058/1.325/1.116 levels.

## Accepted display and limits

Fixed24deg camera, meanbow-110deg azimuth,10m span, equal actual metres and XY
following only. No axes/grid/wiremesh/technicalfooter; compactTime/Roll/Pitch/
Heave HUD. Warm inside, ivory outside, flat face lighting, blue submerged tint;
filled fixedrangeblue water, dark2px actualcached-contact outline. No amplitude
or trajectory amplification. All original movies/previews remain unchanged.

The distantwater ray/color pass is320x225, upsampled for1280x900 presentation;
hull, actualcontact/HUD and water/hull-neighborhood intersections remain native
pixel resolution. Finite cached patch edges remain visible. These are display
limits, not changed waves. Contact uses the SAME diagonal piecewise-linear
cachedwatertriangles; no field regeneration or implied flooding/scattered-fluid
solution. A convex display-volume closure masks water inside the boat, with
actualphysicalsurface boundary checked. All751 directional slope ratios<=.3324
permit unique selectedview ray hits. Maximum rayheight residual4.44e-16m,
contactendpoint residual2.22e-16m. Native linear-model applicability remains
provisional; this renderer does not qualify nonlinear/mesh/time accuracy.

## Execution evidence

Finalsource d6ae724b340c4af9bd111ca7ce003b947a8efe16800b123d8865790aefea43e8;
accepted config84554b883bc775d3c2fdfac8e72c4fd674ae559f3579483c4b06ff86060340cb.
Frozenencodeproposalcbf276209e7af879a17177b89d0296f9b803cffb751ff8a3601f67f736498467,
`presentation-encode-attempt2-source/proposal.json`,24inputs plus dependency IDs.
ROOT grant included full751decode/encodedstills under same serial lease.
Toolsession14472 was recorded UNDERlease withguard12232/producer38888,
renderer34364/ffmpeg35848. Encoding1064s (17m44s); both renderer/ffmpeg exit0
fc41ce, authoritativeprocess/childabsencec01bcc. Frozen24/deps/outputhashesPASS;
idleUNDERleasef92057 -> SEPARATE86c8e3 -> FINALrelease8bbf9c. No ledgerwrites
after release. ffprobe/fulldecode/stills da69a4. No further job or style change.

## Preserved presentation iterations

- Attempt1 proposal1a437400:12views at18/24/30deg and0/30/51.44/60s. Completed
  207.594s; exactdisplaycontact/canvasPASS, ROOT selected24deg. Pythontriangle
  water raster8.86-15.14s/frame, full751forecast4.86h rejected as impractical.
  Source/archive/output preserved. Exact37000/11396terminal; FINALa00a48.
- Attempt2 proposal8e4b908b: fullresolutionvectorizedheightfieldray pass, selected
  24deg. Fourframes23.078s,3.78-6.81s/frame; forecast4242s rejected as too slow.
  Rays/contactroundoffPASS. Exact37168/6508terminal; FINAL645866.
- Attempt3 proposal4a74fa70: distantwater320x225 with nativehull-neighborhood ray
  refinement. Fourframes5.469s,1.187-1.266s/frame; forecast921s. ROOT literal
  selection/receipt, no further source/style changes before fullmovie. Actual
  encode1064s includes allframes+FFmpeg. Exact33504/37064terminal; FINAL9f5804.

Scripts `render_reference_presentation.py`, `prepare_presentation_preview.py`,
`presentation-preview-controls.json` are case/display tools. Oldoutput paths
intentionally fail if reused. Any future render needs explicit bounded scope,
fresharchive/output and ordinary serial claim/verify/terminal/release protocol.

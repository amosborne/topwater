# Minimum native motion and artifact record

User-authorized minimum end-to-end reference. Maintained fork project pin
`f962103a348d63101ce669c2611f7d7c109755c3`; independent selected-path reviews
cover retained excitation history, directional preprocessing sum, File Solid
products-order boundary and actual imported dispatch exactly once. These are
scoped checks, not a whole-library audit or generic accuracy claim.

## Actual mask dispatch and cold start

Updated stock graph had three active imported histories, identical total-force
expressions and positive sum connections at ports1/4/5. The native precomputed
total was therefore tripled. Evidence:
`artifacts/reference-imported-dispatch-stock-attempt4/wiring.json`.

Proper mask API changed one selected rigid-body initialization header to call
`importedExcitationBranches`. The native creation loop for other modes remains
verbatim. Imported groups use one branch; persisted native numbered clones and
their positive Add ports are removed; mixed imported/non-imported groups reject
explicitly. Actual compiled fresh/repeat/ordinary3regular/regular-to-import/
mixed/final-restored tests passed, with accepted total unchanged. Independent
researcher review passed source/execution/hashes and all146 system semantics;
only that authored header changed. Evidence:
`artifacts/reference-imported-dispatch-fixed-attempt2/result.json` and semantic
diff, with exact sources/manifests in its corresponding source directory.

Preserved setup failures distinguish model authoring from physics: an empty
template lacked native physical configuration, and preprocessed body/bus data
must exist before masked body instantiation. Final sequence is input objects ->
Global Reference Frame-only native template -> official preprocessing/bus/variant
setup -> native free joint/body -> explicit reload/update. No simulated warmup
or trajectory was used. A regular-wave test fixture's numeric infinity argument
was corrected to the official H5 `infinite` token; accepted physical inputs did
not change. All failed logs/partial models remain preserved.

## First whole ten-second result

`artifacts/reference/minimum-20261024/case-attempt1` contains501 samples at0.02s
covering0..10s, raw `native-output.mat`, derived `motion.mat`/`motion.csv`, exact
`preflight.json`/`summary.json`, authored model and input file. The source/input/
dependency manifest is `case-attempt1-source/manifest.json`; exact command/log is
recorded there. Whole job terminal exit0; numerical children terminal before
lease release. Total16.74s, preprocessing4.89s/compile3.91s/dynamics5.99s, one
MATLAB engine and one compute thread. This is not an optimization measurement.

Actual H5 CG/CB/volume/K/Ainf readback passed. Compiled dry+added tensor roundtrip
error1.59e-12. Initial q is within8.7e-19m of independently known neutral mean CG;
initial linear/world-angular velocity and radiation force are exactly zero.
No first simulated sample was subtracted to manufacture rest. Prepared and
actual applied501x6 excitation agree with the accepted retained history to
1.72e-12 and1.61e-12 respectively. Raw60-column bus and absolute CG pose remain.
Largest motions are about5.94mm translation and0.0306degree rotation; startup
is retained, with no settling requirement.

The initial apparent `(M+Ainf)*acceleration` mismatch was retained and reviewed,
not hidden. Independent all501 native Solid Newton-Euler residuals are at most
1.71e-13N/2.14e-14Nm, the native force sum is exact, and residual added-mass force
matches stock1e-7s Transport Delay zero initialization/forward extrapolation to
1.18e-12. See the
[independent physical review](../capytaine_expert/native-motion-review-2026-10-24.md).
The coordinator accepted this minimum native dt0.02 result with temporal and
mesh motion accuracy explicitly unestablished. No speculative engine correction
or optional timestep campaign followed.

## Artifact status

First video was completed and preserved at the original case root, but rejected
for water occlusion/header overlap. Iteration2 corrected visibility/telemetry
without physics changes; its video is preserved under `render-attempt2`.
The user requested a closer centered view. Iteration3 uses a fixed12m physical
XY square around full mean hull planform midpoint, true orthographic isometric
camera, equal metres and native mm/degree telemetry. Root inspected its actual
static0/5/10 previews and accepted the centered framing. Final encode3 completed
in71.938s under `render-attempt3-source/manifest.json`. Candidate final video is
`case-attempt1/render-attempt3/reference-isometric.mp4`, H2641280x900,251frames at
25fps,10.04s including the final native t10 sample. SHA256
`a79a75cd88fa03c8b8d46d833cd445645e58a906995a826be9e5674ef345664f`.
Literal encoded frames0/125/250 are decoded separately under its `qa-encoded`
directory, alongside ffprobe.json and source-invariance.json. All five rendering
input hashes and the original draft video hash are unchanged. Encoder/decode
processes are terminal, execution ledger idle and lease released. Independent
capytaine_expert final audit passed literal encoded0/125/250 frames, telemetry,
ffprobe/invariance and all six motion plots: centered hull visibility, labels,
equal metre isometric camera and one CG/rotation mapping are correct. Coordinator
final integration/acceptance remains the final project checkpoint. No accepted
packet, motion, force or body simulation changed for these visual iterations.

Source workflow and exact commands are in
[reference/wecsim/README.md](../../../reference/wecsim/README.md).

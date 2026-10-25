# Mild incoming buoy sea: video source readiness

Scope: new renderer/control files only. No packet, body, native model, fork or
old renderer changes. No numerical field evaluation, preview or encoder has run.
Actual accepted incident/motion identities are pending; source-only controls
fail closed until separately bound and allocated by the coordinator.

Source: `reference/wecsim/render_calm_buoy_reference.py`; controls:
`reference/wecsim/calm-video-controls.json`. Geometry and plot/provenance helpers
reuse the independently exercised `render_generated_reference.py`; forcing/display
consistency uses the consumer-owned `calm_buoy_train.py::CalmBuoyTrain` API.

Python syntax compilation and whitespace checks passed; no rendering/runtime
check is claimed. Accepted-input binding will create separately reviewed controls.

Source review correction: immutable train2 includes an exactly zero-frequency
node with exactly zero quadrature-weighted coefficients. The renderer accepts
that zero DC node only when all its coefficients are zero; every other frequency
remains strictly positive and ordered. No field or rendering evaluation was run.
Source-only controls now point to `train-attempt2/train.npz` and the NEW accepted
`assessment-attempt1/diagnostic.json`. The original train2 diagnostic retains its
earlier-policy rejection and is not used as the visualization acceptance receipt.
Controls remain unbound/pending accepted body motion; no preview/encoder grant.
Current renderer SHA256
`5076954c3c2e957fe8f3faae0e7a58cbb9087669b4bd4dbccea2fc72ea4c363c`;
source-only controls
`94b37cb391980954b1aafa9f7350e32c61341a937bf1dd852d560419056955c7`.
Non-author numerical-methods reviewer subsequently SOURCE PASS at those exact
hashes: zero-DC handling, frame/CG/global mapping, shared-field cache checks,
camera/751native-frame coverage and preview/encoder provenance guards. This is
not actual input freeze, preview acceptance or an execution allocation.

## Physical mapping and view

Original full24ft hull/material-CG/equilibrium quaternion are applied once. Native
position and `Rz(yaw)Ry(pitch)Rx(roll)` radians are in mean-boat axes; the generated
mean global heading rotation/translation are applied once. The wave field is
evaluated in mean-boat coordinates, never moved with instantaneous body rotation.
The source verifies native raw position, motion grid, parent carrier amplitudes,
phases, directions, packet signs and global transformation identities.

Full0..60s native motion is retained. Every fourth0.02s sample gives751frames at
12.5fps; the t60 endpoint is held one frame (encoded duration60.08s). No temporal
interpolation, settling crop, motion exaggeration or wave-height amplification.
The camera follows actual full-hull midpoint XY, with fixed global z and azimuth
mean bow minus135deg. Orthographic35.264deg elevation, equal actual metres,12m XY
footprint and fixed[-2,2]m z range. All displayed hull/surface grids must fit.
The close camera has no effect on open-water hydrodynamics or the upstream train.
Surface labeling explicitly says prescribed incident surface; source Hm0 is
distinguished from the finite train's realized statistics.

## Efficient identical field and checks

Use both signed branches of the EXACT accepted packet's quadrature-weighted C.
Compute per-direction longitudinal profiles with the shared exact spectral/CZT
API in batches of4, then interpolate only the displayed surface on a0.05m grid.
Forcing never uses this cache. The report includes the global interpolation-error
bound `h^2/8 * sum(abs(C)*k^2)` and direct full spectral checks at nine displayed
points at0/30/60s. No new physical precision acceptance criterion is introduced.
Scalar phase reconstruction separately checks world/mean-boat coordinate mapping.
Preview precomputes all751surface frames, writes a hash-bound cache and literal
0/30/60 PNGs plus all3001six-motion plot. Encoding consumes that exact reviewed
cache, avoiding duplicate profile computation. Cache, producer, controls, packet,
heading, raw motion, summary and accepted-review receipt identities stay bound.

Motion acceptance receipt must say `accepted_for_visualization`, reviewer and
observation identity, trajectory status, and exact SHA256 for motion, summary,
native output, heading, packet and incident report. Incident report must be complete
and accepted. Encoder additionally requires a non-author `preview_pass` receipt
binding preview preflight/cache SHA256. No encoder is authorized by source readiness.

## Commands and resources (prospective)

After accepted inputs are bound, source/freeze review and separate allocation:

```
.venv/Scripts/python.exe reference/wecsim/render_calm_buoy_reference.py --config <bound-controls.json> --output <fresh-preview-dir> --mode preview
.venv/Scripts/python.exe reference/wecsim/render_calm_buoy_reference.py --config <same-bound-controls.json> --output <fresh-encode-dir> --mode encode --preview-review <review.json>
```

One CPU thread/affinity, no MATLAB. Preview maximum600s, prospective30-300s;
encoder maximum900s, prospective120-600s, unmeasured for this train. Profile arrays
are batched, all75165x65surface grids are about25MB; design below1GiB, not measured.
FFmpeg single thread/libx264/yuv420p,1280x900. Exact encoder PID/terminal status is
printed; child completion and frozen-input checks precede final report. Failed
or partial outputs remain preserved. Claim/separate verify/run/terminal proof/
idle-under-lease/separate readback/final release apply to EACH allocated job.

Stop: source ready for independent review and later accepted-input binding.
No preview, encoding, motion/field computation or new numerical allocation.

## Actual-input binding preparation

Source/control hashes above retain non-author SOURCE PASS. Immutable packet is
`train-attempt2/train.npz`; accepted field report is the NEW
`assessment-attempt1/diagnostic.json`, SHA
`367a08f33af6fec4aca2cf410d2d8492a294dcaa7ba2e02d6ccefd595113f1e2`.
The earlier train2 report remains rejected under its original policy. Native body
output stays in mean-boat/model-world axes; global201.308deg rotation/translation
are metadata applied once in the display, not applied again inside the model.

After independent body review and exact terminal/release, coordinator writes
`case-attempt1/render-motion-acceptance.json` with these required fields:

```
decision: accepted_for_visualization
reviewed_by: <coordinator and independent reviewer>
observation_id: ndbc-46086-20230319T0840Z
trajectory_status: provisional_model_trajectory
motion_sha256: <case-attempt1/motion.mat>
summary_sha256: <case-attempt1/summary.json>
native_output_sha256: <case-attempt1/native-output.mat>
heading_sha256: <export-attempt1/case.json>
packet_sha256: <train-attempt2/train.npz>
incident_report_sha256: <assessment-attempt1/diagnostic.json>
```

Pure bookkeeping helper is prepared, NOT executed:
`artifacts/reference/calm-buoy-20261024/prepare_render_binding_source_only.py`,
SHA `be515d55db8d9bd75ae7ae699a75fd3593953a156f797853819f33bedf23c7e0`.
It hashes/copies accepted input identities and constructs prospective controls and
one guarded-job proposal; no field/renderer import, lease or process launch.

```
.venv/Scripts/python.exe artifacts/reference/calm-buoy-20261024/prepare_render_binding_source_only.py --stage preview --receipt artifacts/reference/calm-buoy-20261024/case-attempt1/render-motion-acceptance.json --owner /root/calm_video
```

This creates `render-ready-controls.json` and
`render-preview-attempt1-source/proposal.json`, for fresh
`render-preview-attempt1` output (all under `calm-buoy-20261024`). Actual binding
and manifest must be independently checked before any numerical rendering grant.
After literal preview review, the SAME controls/cache identities are used:

```
.venv/Scripts/python.exe artifacts/reference/calm-buoy-20261024/prepare_render_binding_source_only.py --stage encode --preview-review <root-written-preview-review.json> --owner /root/calm_video
```

No accepted-motion receipt, ready-controls, frozen preview proposal or preview
execution is claimed by this source-only preparation record.

## Actual accepted-input binding (no rendering)

Coordinator subsequently wrote the genuine accepted-motion receipt, SHA
`7c52776c9342d3ed69b9089d40d543f357c623586d446f4e396dd846852a6767`, after full
independent3001sample native review and body terminal/release. Receipt scope is
`provisional_model_trajectory`; no mesh/time/nonlinear-accuracy claim.

The bookkeeping helper was promoted to versioned
`reference/wecsim/prepare_calm_video.py`, SHA
`8b4aab8a7e71a90f13c2f5ecf938f2127389bcd0a265e5fe1829d92e067c30d3`.
Only its ROOT parent index changes for its new location. Source controls and
renderer hashes above remain unchanged. Pure accepted-input binding completed
3.407s with no field/renderer import, lease claim or process launch:

```
.venv/Scripts/python.exe reference/wecsim/prepare_calm_video.py --stage preview --receipt artifacts/reference/calm-buoy-20261024/case-attempt1/render-motion-acceptance.json --owner /root/calm_video
```

Bound controls `artifacts/reference/calm-buoy-20261024/render-ready-controls.json`,
SHA `0edd1468c42b5fbcbacec43567e2ea042ac1453e8adbf9350d63c83ae34cad0b`.
Preview manifest `render-preview-attempt1-source/proposal.json` under the same
root, SHA `61d7eaf49dc6f8dfae915fb74177ab1f5704033e07a860a3187f7eea857ac508`.
It freezes15files (14project inputs archived, external FFmpeg binary identified)
and13dependency packages. Fresh output `render-preview-attempt1`; oneCPU/noMATLAB,
maximum600s/expected30-300s. Non-author freeze review and coordinator execution
grant are still required. No preview/cache/PNG/video has been produced yet.
Non-author actual SOURCE/FREEZE PASS subsequently received at proposal61d7eaf4:
all15current inputs,14project archives,13dependency identities and external
FFmpeg hash match; genuine motion receipt verifies all six input identities.
Ready for separate ONE preview allocation; no rendering job performed.

## Allocated preview executed

Coordinator subsequently allocated exactly proposal61d7eaf4. Preview completed
28.609s, producer27940/guard29344/session29828 terminal0/chunk379ef0. Authoritative
CIMec98ba shows both exact clients and other project numerical clients absent.
All15inputs/14archives/13dependencies/output hashes verified; idle UNDERlease
4e1103, SEPARATE readback52cb40, FINAL release1619c4. No post-release ledger write.

Output `artifacts/reference/calm-buoy-20261024/render-preview-attempt1` contains
`frame-0000.png`, `frame-0375.png`, `frame-0750.png`, `six-motions.png`, full751frame
`surface-cache.npz` and final `render-preflight.json`. No encoder was run.

- Preflight SHA `baf5301dc49effa6e63dbb304645c632edc642415eb9be25585cd0b0341495bf`.
- Cache SHA `f6af6151c572f6ddf751b1e7c8eff3d66f4b4f167e7df6f6ccd7a93d0b8b57d9`.
- All751 displayed hull/surface grids fit: hull z[-0.6701,0.8655]m, prescribed
  surface[-0.3042,0.2791]m, hull XY radius3.6636m below6m halfwidth.
- Analytic spatial interpolation error bound56.757um; direct spectral spot checks
  max2.203um; world/meanboat scalar-phase reconstruction2.01e-16m. These are
  display-consistency checks, not body-motion physical accuracy qualification.

ROOT and independent literal/provenance review are required next. Encoder receipt
must contain `decision=preview_pass`, `reviewed_by`, `preflight_json` with the
preflight SHA above, and `surface_cache_npz` with its SHA. Both paths are absolute
or project-relative. Encoding uses unchanged producer, bound controls, input
identities and this exact cache. Pure encoder-manifest preparation and actual
encoder execution require separate coordinator direction. Stop after preview.

## Preview accepted; encoder bound but unrun

ROOT and non-author reviewer personally inspected0/30/60s PNGs and all six-motion
histories; independent all751pose/cache/provenance checks PASS. Genuine ROOT
`render-preview-review.json` SHA
`86e43d7047f32f0c04a3cbd2b7577d066f4a0136c74c59a9659cee92c1381c09`
binds the exact preflight/cache and all four viewed images. Source/code checkpoint
69d2431 was committed/pushed by the coordinator; no physical/source changes.

Coordinator then authorized pure encoder binding only. The versioned helper
completed3.051s, without renderer import/field evaluation/lease/process launch.
Encoder proposal `render-encode-attempt1-source/proposal.json` SHA
`5490d747ebd7638d7932356f1a8f0b21cee859139836b3e8051585317f60104e`
freezes18files/13dependency packages, including unchanged accepted inputs plus
the genuine preview receipt, exact preflight and reviewed surface cache.

Prospective producer:

```
.venv/Scripts/python.exe reference/wecsim/render_calm_buoy_reference.py --config artifacts/reference/calm-buoy-20261024/render-ready-controls.json --output artifacts/reference/calm-buoy-20261024/render-encode-attempt1 --mode encode --preview-review artifacts/reference/calm-buoy-20261024/render-preview-review.json
```

Owner `/root/calm_video`, oneCPU/noMATLAB, maximum900s/expected120-600s. Intended
MP4 `render-encode-attempt1/reference-isometric.mp4`. Final non-author freeze
review and separate coordinator execution grant remain; encoder has not run.

# Generated-buoy-sea renderer: source ready, unexecuted

Bounded ownership: `/root/generated_video` owns only the new Python renderer,
controls and this note. No native case, fork, original renderer, current-state,
handoff or execution ledger was changed. The parent coordinator owns binding,
review acceptance and serial allocation. No render, fixture, encoding or sea/
motion evaluation has run. The January held-release motion remains rejected.

Source: `reference/wecsim/render_generated_reference.py`, SHA256
`6f3f5626d4ec20572faefc78d0491efe3fc589e1c6341b4b586ae02643064572`.
Controls: `reference/wecsim/render_generated_reference_controls.json`, SHA256
`7212083be127f32f40bc2d6c1bf578dd30a292a49aba53d07a5bf89eadbe1ad4`.
AST parsing passes. Runtime/visual/source-independent review remains pending.
Metadata inspected: Python3.12.8, NumPy2.5.4, SciPy1.18.1,
Matplotlib3.11.2, threadpoolctl3.7.0. Eventual execution must freeze actual
dependency/Python/encoder identities under the normal execution manifest.

## Required input binding

The checked-in controls are intentionally unlaunchable: motion, summary, raw
native output and acceptance-receipt paths are null; status is source-only.
Bind a NEW artifact copy of the controls only after an actual case is accepted
for visualization. Set status `ready_for_coordinator_render_allocation` and add
exact hashes for every required input, including the renderer source. Existing
boat, global-bank and heading identities are already pinned. The rejected
Jan90 motion SHA is explicitly excluded, not inferred from its provisional
summary label. No output directory may already exist.

Required schema:

- `motion.mat`: `time`3001native samples0:.02:60; `q`3001x6, measured from known
  neutral meanCG rather than first sample; `position_world`3001x6 contains native
  absoluteCG and Euler radians; `mean_world_cg` is the known three-vector.
- `native-output.mat`: raw `native_body.time` and3001x60signal array; its six pose
  columns must exactly equal the motion file. Summary must report complete60s.
- Unchanged global realization with384amplitudes/angular frequencies/propagation
  directions/phases, headingJSON with native-to-global yawR/t and local bank,
  and canonical full hull. The actual initial pose may be nonzero.
- Acceptance receipt: `decision="accepted_for_visualization"`, nonempty
  `reviewed_by`, matching `observation_id`, `motion_sha256`, `summary_sha256`,
  `native_output_sha256`, `heading_sha256`, and `trajectory_status` equal to
  `provisional_model_trajectory` or `accepted_reference_case`. This records the
  coordinator's existing acceptance; it does not ask the user for permission.

The script neither assumes a trustworthy physical trajectory from finite data
nor independently chooses a canonical case. It retains provisional captions
when that is the status in the receipt.

## Physical display mapping

The physical full hull uses its equilibrium quaternion once, then subtracts
known meanCG once. Each native frame uses `Rz(yaw)Ry(pitch)Rx(roll)` and its raw
CG. HeadingR/t maps the entire moved hull to global east/north/up exactly once.
The sea is evaluated from all384unchanged global components at those actual
global coordinates: `sum a*cos(k.x+phase-omega*t)`. Global/native phase and wave
vector closure, an independent scalar sum, raw native pose identity and known
CG reconstruction are checked in the eventual preflight.

The close view spans12mx12m; orthographic isometric elevation35.264deg and
azimuth meanheading-135deg. Camera follows the actual full-hull midpoint XY,
never vertical heave or yaw. Global z bounds remain[-3,3]m. Equal physical metre
scales, true wave heights and literal rotations are retained. Global axis ticks,
globalCG telemetry and the complete six-motion plot show translations despite
horizontal camera tracking. Native surge/sway/heave and roll/pitch/yaw
perturbations are displayed without first-sample rebasing. The header explicitly
states the mean global bow heading and the chosen incoming-relative offset.

The65x65transparent incident grid is an analytic sampling of the full global
bank; it is not a scattered/radiated solution. All751displayed fields and hull
poses are preflighted for clipping. These are displayed-grid bounds, not a
continuous-space or continuous-time physical envelope. A fixed color scale
of±3sea standard deviations changes colors only; elevations are unmodified,
and saturated color counts are reported. Hull overlay/alpha0.40 follows the
accepted baseline's legibility convention; visual review must confirm it.

## Separate preview and encoder commands

After root allocation, input binding and normal claim/separate verification:

```powershell
.venv/Scripts/python.exe reference/wecsim/render_generated_reference.py --config artifacts/reference/generated-render/bound-controls.json --output artifacts/reference/generated-render/preview-attempt1 --mode preview
```

Preview writes0/30/60sPNGs, `six-motions.png`, complete preflight/provenance and
output hashes. No encoder runs. Coordinator/independent reviewer must inspect
the literal images, including physical heading, scene legibility, telemetry,
full boat visibility and equal spatial units.

Encoding is a separately allocated fresh output:

```powershell
.venv/Scripts/python.exe reference/wecsim/render_generated_reference.py --config artifacts/reference/generated-render/bound-controls.json --output artifacts/reference/generated-render/encode-attempt1 --mode encode --preview-review artifacts/reference/generated-render/preview-review.json
```

The preview review receipt contains `decision="preview_pass"`, nonempty
`reviewed_by`, `preflight_json` path and its `preflight_sha256`. The reviewed
preflight must be complete preview output with exactly matching source, controls
and input identities. Encoder records that receipt's identity too.

Each display frame takes every fourth native sample, with no motion/time
interpolation:751frames at12.5fps, t0..60 inclusive. This display cadence is
chosen for the incident band up to0.533Hz and expected boat motion around1Hz;
it is not a low-pass filter or a time-convergence claim. All3001native samples
remain in the six-motion plot and source data. Encoded duration60.08s includes
the last t60 endpoint held for one frame. H2641280x900/yuv420p/one FFmpeg thread.
The child encoder PID and exact waited exit are printed for execution evidence.
Input/config/preview-review invariance is checked before marking output complete.

Resource proposal: one allowed logicalCPU, estimated preview10-90s and
encode150-450s, both unmeasured for this full384component path. Declared caps:
preview180s, encode900s; expected memory<1GiB. Baseline251frame rendering took
71.94s, so751frames alone suggests about215s before additional preflight cost.
The analytic field is factored into bounded64frame matrix batches; no component
is filtered. Deadline failures retain incomplete artifacts and require a new
coordinator decision, not an automatic rerun or reduced physical resolution.

After encoding, independently inspect FFmpeg-decoded0/30/60sframes and a frame
at actual roll extremum, ffprobe codec/resolution/frame count/duration, the full
six-motion plot and source/data invariance. This source assignment stops here:
no accepted trajectory, job grant, preview receipt or physical acceptance exists.

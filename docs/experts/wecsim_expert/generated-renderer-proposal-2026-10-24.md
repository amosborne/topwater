# Generated-buoy-sea renderer:10m preview complete, encoder held

Bounded source ownership: `/root/generated_video` owns the renderer, controls
and its notes. Native case/fork/originalrenderer remain untouched. A subsequent
explicitly allocated12m preview1 completed, including its required ledger
bookkeeping; its terminal receipt is preserved in the seeded binding note. No
encoder has run. The10m revision independently passed source review, then
completed separately allocated preview2; receipt is in the seeded binding note.
Independent literal-preview acceptance and a separate encoder grant remain.
The January held-release motion remains rejected.

Source: `reference/wecsim/render_generated_reference.py`, SHA256
`82083d9b5936ce45d93cc01db83c6fd26fcb33f26bbe001254a09cc76d20bb26`.
Controls: `reference/wecsim/render_generated_reference_controls.json`, SHA256
`df16ae8078c54dac66657f9e02953913653230cc34bc3c468e43c131a2fdea37`.
AST parsing and independent10m source review passed. Preview2 completed with
clipping/phase/data-invariance checks; literal-preview acceptance is separate.
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

The close view spans10mx10m; orthographic isometric elevation35.264deg and
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

Resource proposal: one allowed logicalCPU. Actual10m preview2 took1.797s;
encode estimate150-450s remains unmeasured. Declared caps:
preview180s, encode900s; expected memory<1GiB. Baseline251frame rendering took
71.94s, so751frames alone suggests about215s before additional preflight cost.
The analytic field is factored into bounded64frame matrix batches; no component
is filtered. Deadline failures retain incomplete artifacts and require a new
coordinator decision, not an automatic rerun or reduced physical resolution.

After encoding, independently inspect FFmpeg-decoded0/30/60sframes and a frame
at actual roll extremum, ffprobe codec/resolution/frame count/duration, the full
six-motion plot and source/data invariance. The initial source-only assignment stopped before any job. Subsequently root
approved provisional visualization and preview2 ran under its separate grant;
see the binding note. No encoder or physical acceptance follows automatically.

## View-only revision after the first preview

Root inspected the actual12m preview and requested a closer10m footprint. The
initial mean-heading-minus135deg camera was already correct and remains FIXED;
only XY limits/aspect change from12m to10m. Equal units, z[-3,3], motion, sea bank,
heading mapping and tracking behavior are unchanged. The maximum displayed hull
horizontal distance from the tracking center in preview1 was3.7124m, below the
new5m halfwidth. Separately allocated preview2 checked all751new displayed frames and grids.
Previous preview/source archive and ready controls remain untouched. No encoder
is authorized by source review or by a successful preview alone.

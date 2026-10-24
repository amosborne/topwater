# Seeded renderer:10m preview2 complete, encoder held

The initial `/root/generated_video` metadata-only binding ran no MAT-array or
sea-field/plot/render evaluation. Root subsequently accepted PROVISIONAL
visualization and allocated12m preview1, whose receipt is below. Current10m
view revision completed as separately allocated preview2; no encoder grant.
Native motion and
physical data are unchanged; no physical-envelope or numerical qualification
is implied. Independent actual-output implementation review passed; the first
geometry audit is incomplete123/3001 and therefore not a whole-case pass.

New controls: `reference/wecsim/render_generated_reference_seeded_controls.json`.
SHA256 `70f22974dfe47df91fa02b75be70196765fb7d936972475325a1461d4ba744a5`.
Status is deliberately `seeded_inputs_bound_pending_root_visualization_receipt`,
which the renderer rejects before loading/evaluating the motion. The intended
receipt path is present but its hash is not bound. Root must supply the actual
decision; an embedded template has `PENDING_ROOT_VISUALIZATION_DECISION` and
`reviewed_by=null`. Nothing in this binding authorizes a numerical job.

The current10m view revision is source
`82083d9b5936ce45d93cc01db83c6fd26fcb33f26bbe001254a09cc76d20bb26`,
base controls `df16ae80...`; its tiny view-only diff independently passed review.
Original12m source `6f3f5626...` and controls `7212083b...` are preserved in the
preview1 source archive. Kinematics, global sea and physical scaling are
unchanged. The old rejected held-release motion `190deabd...` remains explicitly excluded. New motion must
not be confused with `case-jan90-attempt1` or the calm control.

## Actual candidate identity

Directory:
`artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1`.
Native body record:
[seeded-native-case-proposal-2026-10-24.md](seeded-native-case-proposal-2026-10-24.md).

| Input | SHA256 |
| --- | --- |
| motion.mat | d69652e3ca838afd2df6ddc249762d58321d833324038c51798399a36d055a95 |
| summary.json | 9cb8dcc38b91a77a4c860329061bbe681fe9944487bcc3e6e8d3aa99c501bbd3 |
| native-output.mat | fe30d8eae1e082a44538cb8a40a5300098becd1cd63dbb4c63a07dbc998758b6 |
| native-cic-output.mat | 2d889ef72e1a896fe6d10de9fbc0d826819b8fcf96ff4547bf0655c76c6d7d1e |
| initialization.json | 2140982553f6047dc917931d1a6ccdec50c6e707c592cb2dc1515ed3213ceec7 |
| preflight.json | 1428b93bc99960dcaadc8949276e377ebad1f803a30a19b4f5c6fbedeae2ad5e |
| prepared-case.mat | 61a604243283d5dba712ba6a0257a07ec784ef3caaa1fe38d59847af222ae2e3 |

The actual case reports3001native samples over60s, nonzero initial q, CGvelocity,
worldomega and radiation memory. It preserves actual post0 free motion. The
renderer must retain those native values, never subtract its initial sample.
Observed-record ID, seed42, full384global harmonics and heading-offset90 map
remain the pinned inputs from the unbound source controls. Global mean bow
heading is4.4194798085094735rad; R/t is applied once after native pose mapping.
No motion or sea amplitude, phase, filtering or trajectory change is involved.

## Root completion of the binding

After actual-output/geometry review and a coordinator decision to visualize:

1. Write the intended `render-motion-acceptance.json` in the actual case folder,
   using the controls' embedded template. Set `decision` to
   `accepted_for_visualization`, provide the actual reviewer identity and revise
   physical-limitations text to reflect the completed evidence. Keep
   `trajectory_status="provisional_model_trajectory"` unless separately qualified.
   Exact motion/summary/raw/heading hashes are already filled in the template;
   the CIC hash is additional provenance, not a force used by the renderer.
2. Make a new artifact copy of the controls as
   `case-jan90-seeded-attempt1/render-ready-controls.json`, bind the receipt's
   path/hash into `input_sha256`, and set status to
   `ready_for_coordinator_render_allocation`. Preserve the versioned pending
   binding. Archive this actual ready copy and input/dependency identities in
   the numerical execution proposal before the normal serial claim/readback.
3. Allocate only the preview first. Its future literal images and clipping/
   frame/phase checks need independent review before the encoder receipt/grant.

Exact preview command after those steps and numerical allocation:

```powershell
.venv/Scripts/python.exe reference/wecsim/render_generated_reference.py --config artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-ready-controls.json --output artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-preview-attempt1 --mode preview
```

This creates a fresh directory and0/30/60s previews, the full3001-sample six-motion
plot, preflight and output hashes. Output directories have not been created by
this binding task and must be checked fresh at launch. Preview proposal:
oneCPU, maximum180s, estimated10-90s and<1GiB, unmeasured for this candidate.

Only after preview acceptance and a separate encoder allocation:

```powershell
.venv/Scripts/python.exe reference/wecsim/render_generated_reference.py --config artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-ready-controls.json --output artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-encode-attempt1 --mode encode --preview-review artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-preview-review.json
```

Encoder proposal: oneCPU/oneFFmpeg thread, maximum900s, estimated150-450s and
<1GiB, unmeasured. All751frames take every fourth native sample0:.08:60 at
12.5fps; encoded60.08s includes the final native endpoint hold. The complete
native histories remain in the plot. Viewport remains equal metres10mXY,
tracks physical hull midpoint horizontally, fixed global z[-3,3], and retains
all six native and globalCG telemetry. Actual global sea coordinates follow
the camera without rephasing the immutable bank.

Known limitations remain captioned provisional. No field/scene/clipping result
exists yet. A render is not a geometry, passivity or numerical-accuracy test.
Stopping point: input hashes and receipt template ready for root; NOT READY for
execution. Geometry owns the next separately allocated numerical slot.

## Executed preview1 and pending10m preview2

Separately allocated provisional12m preview1 completed2.344s, terminal0. It wrote
frame-0000.png(t0), frame-0375.png(t30), frame-0750.png(t60), and the complete
3001-sample six-motion plot in `case-jan90-seeded-attempt1/render-preview-attempt1`.
All751displayed hull/sea clipping checks, raw native-grid/pose/knownCG identity,
global/native phase closure2.22e-16m and scalar-vs-factored field2.66e-15m passed.
Displayed hull z[-1.9472,1.6590]m, sea-grid z[-1.6421,1.1717]m and maximum hull
XY distance3.7124m were retained. Source/data and2492dependency implementations
were unchanged; no encoder ran. Author literally inspected all three images and
the plot: labels/telemetry/hull are legible; this is not independent acceptance.
Root's literal inspection requested a closer view, not a physical data correction.

Exactproducer31792/guard9188 terminal0 and authoritative absence were verified.
Idle771639, SEPARATEidle063ddb, FINALrelease261a72 all exit0. Evidence:
`render-preview-attempt1-source/proposal.json` SHA7cf5ad144c01984022af6b11d8379ff4d9c283272b0ff44be7895aa2fb35bfba;
`render-preview-attempt1/render-preflight.json` SHA213166344a7b5457e474f3b03744802653a40803bcd147378083162b08b52d02;
postcheck SHA24069e2f5fb96b631ef8f48577b5bb31b8de334c56bcb0e9642b59dc9c394c19.
These outputs and their archived12m source remain unchanged.

NEW artifact controls `case-jan90-seeded-attempt1/render-ready-controls-v2.json`
SHA861b1282c78c36f8a61d54c5fade3fdd014d64c562ff3f48bb458467161e7596
pin the10m view source and the SAME root provisional-visualization receipt
SHA6cdb258e1aecfa09b78d15a7259dc4c11070e1db1c3f769d87b99c80dea9e4b2.
The old ready controls SHA2811d87b... remain unchanged. The camera was already
fixed at initial global meanheading-135deg (118.21754deg); it is retained and
now explicitly identified in output metadata. Only12m->10m limits/aspect changed.
At that source-only checkpoint no field/render had run for the10m revision;
its later allocated preview2 receipt is below. No encoder has run.

Prospective fresh preview2 command, ONLY after tiny diff review/root grant:

```powershell
.venv/Scripts/python.exe reference/wecsim/render_generated_reference.py --config artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-ready-controls-v2.json --output artifacts/reference/generated-bank-consumer-20261024/case-jan90-seeded-attempt1/render-preview-attempt2 --mode preview
```

Same oneCPU/180s preview cap and independent literal-preview review before any
encoder. Fresh output is required. The original held-release case remains
rejected; provisional visualization still does not accept the whole physical
operating envelope, radiation passivity or numerical accuracy.

## Executed10m preview2 receipt

ONE root-allocated preview2 completed1.797s on oneCPU (guard tool6.098s), terminal0,
using source82083d9b/readycontrolsv2861b1282/unchangedrootmotionreceipt6cdb258e.
Output `case-jan90-seeded-attempt1/render-preview-attempt2` contains literal
frame-0000.png(t0), frame-0375.png(t30), frame-0750.png(t60), and six-motions.png
with all3001native samples. Author inspected all four images: full hull remains
legible/centered with clear telemetry and no header overlap. Independent literal
review/root acceptance is still required before any encoder.

All751displayed frame checks pass: hull z[-1.947155,1.659002]m, sea-grid
z[-1.610300,1.142218]m and maximum hull XY distance3.712425m<5m halfwidth.
Scalar-vs-factored global elevation error2.44e-15m; global/native closure6.66e-16m.
Raw native grid/pose and knownCG identities pass; full384bank, source, controls,
receipt and2492dependency implementations stayed unchanged. Projection remains
orthographic/equalmetres; fixedinitialmeanheading-135deg azimuth118.21754deg,
10mXY midpoint tracking, fixedglobal z[-3,3], no amplified motion or sea height.
These are displayed-grid/pose checks, not physical immersion qualification.

Execution evidence: Claim1dc151, SEPARATEverification6ba147, exactrun727b38exit0.
Producer31120 andguard35800 waited terminal0 and were authoritatively absent
(c25f3d). Idleunderlease1e39d6, SEPARATEidleproofe64dd6, FINALrelease2af20f all0.
No encoder or other numerical child spawned; no next numerical job granted to
this owner. Geometry owns the next separately allocated slot.

- Proposal `render-preview-attempt2-source/proposal.json` SHA89bbaf0d87581dea79e17fd9779d8943d31bf259781d182b78a2d1091d5936a9.
- Preflight `render-preview-attempt2/render-preflight.json` SHAba7d7f9fa80179782a362f1361f71d93d9cb5a1a5be58669d5ff6fcf3353c7d5.
- Postcheck SHAe8e13f500d157af7371490e213e6a1cfd4fbaa8d97cf0c37b12991523bd90bd8.

Stop: preview2 artifacts ready for independent literal review. No video exists
for the seeded case yet, and no physical-envelope/accuracy claim follows.

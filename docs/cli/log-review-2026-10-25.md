# Generated-log review, 2026-10-25

Independent review of existing text only: no numerical imports, tests, solves or
production changes. This review does not qualify the physical results. The
logging owner is concurrently correcting Rich stream capture; this document
does not prescribe another capture implementation.

## Evidence and observed behavior

- [Canonical sea log](../../examples/simple-boat/sea.log): compact and useful.
  It records the observation ID, seed, 384 components, publication and successful
  end. Keep this scale of detail.
- [Interrupted public solve log](../../examples/simple-boat/result.log), beginning
  20:53:42 UTC: the first Rich warning becomes a separate local-time record,
  `WARNING solve WARNING`, then many INFO records such as `RadiationProblem`,
  `(`, `body`, `=`, and `.500`. The radiation warning has the same defect.
  Terminal formatting fragments have become log records; severity no longer
  reliably identifies the warning's content. The duplicated local and UTC
  clocks also distract. This is an observed defect, not hypothetical noise.
- [Archived attempt 1](../../artifacts/workflow/canonical-pipeline-attempt1/output.log)
  ends with a locked `KH.dat` cleanup error and no visible original export
  cause. [Attempt 2](../../artifacts/workflow/canonical-pipeline-attempt2/output.log)
  improves this substantially: failed operation, frame 277, time 22.16 s,
  exception, traceback and retained replay path are all present. Preserve this
  distinction between primary failure and cleanup/recovery details.
- [Attempt 3](../../artifacts/workflow/canonical-pipeline-attempt3/output.log)
  ends after radiation transfers; its wrapper records exit 4294967295 and
  320.329 s. The canonical `result.log` has no END. These are interrupted-run
  evidence, not completion. The process was deliberately stopped externally;
  the log alone cannot establish that reason.
- [Real two-frame render fixture](../../artifacts/workflow/logging-geometry-tests-attempt1/pytest-temp/test_real_render_cli_records_v0/video.log)
  records load/render stages, publication and successful completion. Its single
  `RENDER timings={...}` line includes a repeated output path, Python tuple/dict
  syntax and values such as `0.7179999999934807`. The information is useful but
  harder to scan than necessary. This is fixture evidence, not a completed
  canonical render. [Missing-result](../../artifacts/workflow/logging-geometry-tests-attempt1/pytest-temp/test_render_missing_result_is_0/video.log)
  and [corrupt-result](../../artifacts/workflow/logging-geometry-tests-attempt1/pytest-temp/test_corrupt_result_is_clear_f0/video.log)
  fixtures clearly record ERROR, cause, failed END and no publication.

## Priorities

| Priority | Small improvement | Reason / acceptance |
|---|---|---|
| Must fix | Capture complete logical Rich lines and retain warning severity throughout the warning block. | The current fragmented records obscure both message and severity. Check actual Windows-style segmented flush behavior; a bare `WARNING` plus INFO tokens is insufficient. Do not count individual write fragments as warnings. |
| Must retain | First physical warning plus per-family affected-solve or affected-sample summary; distinct primary error and recovery location; truthful publication/END status. | Existing summaries are already useful. Preserve the attempt-2 frame/time/cause/replay details. Missing END after forced termination must remain incomplete, never reconstructed as success. |
| Next, small | Time the actual long stages explicitly: coefficients, MATLAB preparation/integration/audits, radiation, diffraction, surface geometry/diagnostics, final validation/publication. | The public log currently labels nearly all work `solve motion and total surface`; UTC timestamps help, but readers must infer transitions and durations. Use wall elapsed seconds at existing boundaries, separate from physical simulation time. Avoid timing every helper. |
| Next, small | Give the exact-field/prepared-load discrepancy a WARNING/category and readable force/moment vectors with limited significant figures. | Its current `Surface:` line reports [20.7354, 18.0728, 47.5506] N and [21.6547, 90.0598, 16.2102] Nm but presents as ordinary progress. This is a known provisional consistency limitation; no new rejection threshold or changed loads follows. |
| Next, small | Replace the render Python dict display with one video summary and a short rounded timing line. | Retain frames, fps, dimensions, source endpoint, encoded duration and verification; round timings to useful precision. Preserve both source and encoded durations: two samples at 12.5 fps legitimately encode 0.16 s while their source endpoint is 0.08 s. |
| Next, small | Report progress periodically by elapsed wall time within slow stages, with completed/total work and stage elapsed time. | Current radiation `101/960` counts are useful. Geometry every 100 frames can leave several minutes without a status at the measured old throughput. Do not print every frame. An approximate stage ETA is optional after enough completed work; do not promise a whole-command ETA from one stage. |
| Nice to have | One concise context line with Topwater identity, selected backend/version and input identities. | BEGIN already records directory, duration, backend, or observation/seed. Resolved paths and optional input SHA256 values would connect a log to the exact numeric files without NPZ metadata. Avoid repeated package dumps, environment variables, per-array hashes or a new manifest format. Do not add dependency imports merely to decorate a log. |

For warning counts, keep the documented distinction between warning **records**,
affected BEM solves and affected physical samples. A first warning and its summary
are two records about one family. A family label is more useful than an
unqualified total, but no new nested summary format is necessary.

## Physical meaning and severity

Keep the numerical and physical warnings visible: coarse mesh, irregular-frequency
risk, 42.7 mm maximum / 12.5 mm RMS projection sensitivity, restoring departures,
unavailable calm-reference comparisons, field/load discrepancy and eventual
wet-set departure. Their meanings differ; [warning triage](warning-triage.md)
already provides the necessary distinctions. Neither "expected" nor a successful
END makes these harmless or establishes accuracy. The projection comparator is
not a reference solution or error bound; restoring percentages are not dynamic
error percentages; flat-reference rim intersection is not actual-sea flooding.

Known initialization progress can be INFO after specifically identifying it:
Capytaine tabulation and visualization-only 200 m depth are examples. Do not
demote unknown third-party warnings by a broad text rule. The Windows path-display
defect is fixed separately; the fresh owned scratch-directory warning is lifecycle
noise, not evidence of a physical failure. Refresh the triage document's current
status from completed evidence: its full-sequence-pending/frame-277-blocker wording
is historical now that the cached 751-frame sweep passes, while public workflow
completion and physical qualifications remain pending.

The final successful solve should retain its existing total-minus-incident
metrics with metres, area/time weighting, near-waterline region definition and
complete frame coverage. It should also report initial represented total eta
separately from the full-hull **incoming** eta/velocity gates. No missing metric or
unavailable comparison should read as a pass. These summaries were not reached
in the reviewed public failures.

## Minimal final presentation

Keep the current UTC/level/command prefix. The following is a layout sketch,
not evidence from a completed run; placeholders must come from actual work:

```text
INFO solve BEGIN directory=... duration_s=60 solver=wecsim
INFO solve STAGE radiation transfers started frequencies=960
INFO solve PROGRESS radiation transfers completed=401/960 elapsed_s=...
WARNING solve [mesh-resolution] ... first warning and physical limitation ...
WARNING solve [mesh-resolution] radiation transfers: 5010 solves affected
INFO solve STAGE radiation transfers complete elapsed_s=...
INFO solve STAGE surface geometry started frames=751
INFO solve PROGRESS surface geometry completed=.../751 elapsed_s=...
INFO solve Surface total-minus-incident metrics ... units/weighting/region ...
INFO solve PUBLISHED .../result.npz
INFO solve END status=success elapsed_s=... output_published=true warnings=...
```

On a handled failure, retain an ERROR with stage/frame/time when known, the cause,
whether an old output was preserved, any usable recovery path, and failed END.
On a handled keyboard interrupt, use interrupted status if supported and preserve
the same publication truth. An externally killed process cannot reliably write
its own final line: keep the missing-END interpretation and the external terminal
record; do not retrospectively edit its generated log. No automatic retry, resume
command or recoverability should be promised unless actually implemented.

Stop point: recommendations delivered to the coordinator and logging owner.
Production implementation, tests and completed canonical log acceptance remain
with their assigned owners.

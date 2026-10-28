# Canonical workflow warning triage

## Current accepted 130-second packet run

Implementation `e416ca8` passed 170 focused Python tests in 56.93 s. The new
case uses 10 s rest, 5 s ramp, 100 s dwell, 5 s ramp and 10 s rest, with a
0.01 s native time step. All 13001 native samples and the native phase,
startup, initial state, applied-force, inertia, equation-of-motion and CIC
assertions completed. The ordinary invocation failed afterward, during engine
cache-path restoration, before publishing a new Result.

The saved motion was validated and reused. Recovery2 passed a real restoration
control, matched the rebuilt operator against all 14 numeric native arrays and
all netCDF variables, and computed the complete new total field. Geometry-only
replay passed all 1626 frames and published the accepted Result/PNG pair at
2026-10-28 05:04:00 UTC. The ordinary CLI video passed at 05:21:46 UTC: 1626
decoded frames, 800 by 800 pixels, 12.5 fps, source time 0–130 s, encoded duration
130.08 s and zero video warnings. All eight example files are now accepted.
No body or field solve was repeated for replay. The original ordinary solve
still records its cleanup failure; the successful Result log explicitly records
the private continuation.
See the [packet integration run record](../experts/packet-integration-run-2026-10-27.json).

| Warning or condition | Meaning and disposition |
|---|---|
| Generated packet sidebands | Expected under the approved construction policy. Cutoff is 0.984615 Hz; 59 nodes exceed the original carrier operating limit of 0.533459 Hz. Their amplitude L1 contribution is 0.107014 m, or 0.5749% of total L1. They retain genuine coefficient coverage and are not truncated. This does not qualify their mesh resolution or linear applicability. |
| Mesh resolution / irregular frequencies | The new coefficient producer flags 1962 of 2388 solves for resolution and 822 for irregular-frequency risk. Further field evaluations retain their own warning records. These are real, unqualified numerical limitations under [#32](https://github.com/amosborne/topwater/issues/32), not measured motion or elevation errors. |
| Native direction-coverage heuristic | The native transcript warns about direction coverage/coarseness. Its endpoint heuristic differs from actual coefficient-node membership checks, which passed, as did phase and applied-force audits. This does not establish angular coefficient convergence; retain the limitation under #32. |
| Mean-sea-reference restoring departure | Moment: 8612 samples, first 0 s, last 130 s, peak 195%. Force: 5054 samples, first 13.21 s, last 117.46 s, peak 181%. These compare exact flat mean-sea-reference restoring with linear `-K*q`; they are substantial applicability departures for [#34](https://github.com/amosborne/topwater/issues/34), not nonlinear motion-error percentages. |
| Restoring diagnostic unavailable | 2459 samples, first 18.06 s, last 116.52 s. The mean-sea reference loses supported contact or reaches the opening. This means the diagnostic cannot assess those poses. It does not establish flooding or loss of contact with the actual reconstructed sea. |
| Represented wetted surface / actual opening contact | Wetted-area symmetric difference exceeds 10% from 14 to 116.4 s; peak 25.13% at 43.36 s. All 1626 actual reconstructed-surface opening/contact checks passed. The wetted-area departure remains a fixed-mean applicability concern under #34, not a numerical error percentage. |
| Field/load consistency | The rebuilt operator agrees with the retained operator; this is input identity, not force-field closure. [#33](https://github.com/amosborne/topwater/issues/33) remains separate from mesh accuracy and applicability. |
| Green-function tabulation / fresh scratch export directory | Operational setup messages. The producer creates its own unique directory before the legacy exporter warns that it already exists. Retain freshness checks; neither message alone establishes a physical failure. |
| Root-owned MATLAB cache restoration | The original invocation failed because TMP was beneath its initial cache folder. A recovery1 control then failed a raw string comparison of equivalent path representations. Recovery2 passed the actual production restoration control using sibling cache/code-generation folders and normalized path identity; the owned engine quit before field work. This resolves that control's operational failure, not physical accuracy or a fresh ordinary-solve acceptance. |
| Render verification budget | The first render exceeded a fixed 240-second budget. The CLI budget now scales as `max(240, waterframes)`; 25 CLI tests passed in 3.99 s and the final ordinary render completed in 303.813 s, including 273.828 s of render core work. This was an operational limit, not a physics change or an accuracy gate. |
| Peak-location log coordinate | The reported peak Z, 0.4378107785 m, remains in the computational mean-water frame. Its world Z is 0.4499780082 m after adding the stored DC offset. Published fields and delta metrics already use the correct transformation; this is a nonblocking diagnostic presentation issue. |

The native phase difference is `1.7583157152500917e-14 m`. The preserved static
sea-level term is `0.012167229735326152 m`; replay restores it once to motion and
surface. The complete field cache is bound by SHA-256
`9f3870fb9f3e91deb9b82ea83b25f62ed24e6acb43ad11dfd0d1acf3157660fb`.
New publication logs inherit prior warnings with their source attribution;
private copies preserve the original failed/interrupted logs.

The earlier 60-second case's 1 cm numerical surface qualification does **not**
transfer to this new packet. Passing runtime identities and balance checks does
not establish mesh convergence, surface error against the linear model, or
nonlinear accuracy. Issues #32–34 remain distinct follow-up work.

### Measured new surface disturbance and postprocessing

The following differences are **total-minus-incident surface disturbance**, not
reconstruction errors. Exterior support excludes the posed hull; the one-beam
band uses the existing source-triangle selection and projected-area moments.

| Region | RMS disturbance | Reported incident RMS | RMS ratio |
|---|---:|---:|---:|
| Full exterior patch | 6.83618 mm | 0.384527 m | 1.77782% |
| One-beam band | 8.03576 mm | 0.382852 m | 2.09892% |

The resolved absolute peak disturbance is 51.6636 mm at 24.48 s. Geometry
postprocessing took 934.221 s: wetted-set diagnostics 617.859 s, delta metrics
266.193 s, opening checks 27.780 s and waterline calculation 21.731 s. These
timings explain postprocessing cost; they do not imply another physics solve.

## Historical 60-second canonical result — issue #31 qualification

Before the current rerun, the eight files in `examples/simple-boat/` retained the accepted WEC-Sim motion
and the same boat/sea inputs. The corrected full-complex reconstruction supplies
all 751 total-surface frames; all 3001 motion samples are unchanged. The refreshed
video passes complete decode, dimensions and cadence checks: 751 frames,
800 × 800, 12.5 fps, physical time 0–60 s, encoded duration 60.08 s.

All three [result](../../examples/simple-boat/result.log),
[sea](../../examples/simple-boat/sea.log) and [video](../../examples/simple-boat/video.log)
logs end successfully with no ERROR records. Sea/video have zero WARNING rows;
result has seven: six inherited calm-restoring records and one newly evaluated
wetted-set warning. These are records/summaries, not seven independent events.
The publication log explicitly identifies cached qualified surface publication;
it does not claim another WEC-Sim solve, BEM evaluation or geometry replay.
Its video-pending statement describes publication time; the later video.log
records the successful separate render. Logs are not edited retrospectively.

| Item | Meaning and disposition |
|---|---|
| Surface reconstruction | All 5929 targets × 751 frames pass the 10 mm numerical linear-model requirement: maximum runtime/reference difference 0.341 mm; maximum pointwise difference plus reference uncertainty 5.784 mm. The old 42.7 mm cosine-projection discrepancy belongs to the retired method and is not a current warning. See [qualification and assumptions](../experts/numerical_methods_researcher/issue31-selected-case-acceptance-2026-10-27.md). |
| Calm restoring departure | Six inherited records include moment: 2013 samples, first 9.68 s, peak 183%; force: 1170 samples, first 16.66 s, peak 132%. These compare fixed-mean restoring with a flat calm-water diagnostic, not nonlinear dynamic error. Expected for the unchanged trajectory, but substantial and needing evaluation under [#34](https://github.com/amosborne/topwater/issues/34). |
| Calm diagnostic unavailable | Included in those six records: 338 samples, first 23.08 s. The flat reference intersects the opening rim. This is unavailable diagnostic coverage, not evidence of actual-sea flooding. Every actual reconstructed-surface opening/contact gate passed. |
| Represented wetted surface | One new record: symmetric-difference area exceeds 10% from 29.04 to 53.12 s, worst 18.27% at 44.64 s. Expected fixed-mean applicability departure, retained for further evaluation under #34. This is not a numerical reconstruction failure or a nonlinear motion-error percentage. |
| Mesh resolution / irregular frequencies | Cached publication performs no new BEM evaluation and does not re-emit old producer warnings. Its log links the actual [source evaluation](../../artifacts/workflow/issue31-full-allocation2/child.log). The unchanged coarse mean-hull representation remains unqualified under [#32](https://github.com/amosborne/topwater/issues/32); passing #31 does not resolve that issue. |
| Field/load consistency | The surface and applied-force workflows retain their existing relationship. [#33](https://github.com/amosborne/topwater/issues/33) remains open; unchanged motion and numerical reconstruction accuracy do not establish force-field closure. |
| Native direction/setup information | Prior native transcript retains the periodic-grid direction-coverage heuristic warning and visualization-only 200 m depth information. Actual component/force preservation checks passed; infinite-depth physics is unchanged. No native dynamics rerun occurred in this publication. |
| Encoder messages | Current video.log has no warnings and complete verification passes. Any external FFmpeg duplicate-rate notice in the retained transcript is configuration noise, not a failed cadence check. |

### Measured surface disturbance

These are **total-minus-incident elevation differences**, not numerical or
physical errors. Exterior support excludes the posed hull. Spatial P1 moments
use projected-area weights; time integration includes all saved frames/startup.
The near band selects original source triangles by centroid distance within
one beam of the actual waterline and integrates their exterior portions.

| Region | RMS difference | RMS / incident RMS | Resolved maximum absolute difference |
|---|---:|---:|---:|
| Full exterior patch | 5.08038 mm | 1.90312% | 39.1241 mm |
| One-beam band | 5.99504 mm | 2.22013% | 39.1241 mm |

The resolved peak is at 43.36 s, approximately
(2.800242, 0.207352, -0.171332) m world coordinates. The surface correction
changes these disturbance metrics; it does not change solved boat motion.

The reference uncertainty includes numerical controls and an omitted-frequency
allowance conditional on the explicitly calibrated response envelope. It does
not estimate nonlinear or experimental error. Warnings concerning the underlying
hydrodynamic data and physical applicability remain open. No full-hull total
disturbed-fluid velocity bound has been established.

## Archived attempt-2 assessment

The following is the earlier failed-run assessment, retained with its original
evidence and policy. Its failure is historical; the current successful uncut
workflow status is recorded above.

Reviewed 2026-10-25 against the actual
[attempt 2 log](../../artifacts/workflow/canonical-pipeline-attempt2/output.log)
and its [terminal record](../../artifacts/workflow/canonical-pipeline-attempt2/terminal.json).
Capytaine 3.0.0; maintained WEC-Sim 7.1.0. This was a read-only assessment of
textual source and saved evidence: no calculations, tests, imports or numerical
jobs were performed, and no physical/source correction was authored.

The command exited 2 during required surface export. It did not publish a
successful `result.npz` or render a new video. Existing results are not acceptance
evidence for this attempt. The table distinguishes the current workflow blocker
from limitations that remain visible under the selected provisional policy.

The Windows path-display defect below is now fixed in maintained-fork commit
`9952634f6f2e2746a1626484e27c3003d1e1d5a9`. Its
[three-path regression](../../artifacts/workflow/wave-display-check-attempt3/report.json)
prints Windows, percent-containing and UNC paths verbatim with zero warnings.
That presentation fix does not change the physical/numerical qualifications.

The historical frame-277 sliver has a guarded triangulation correction.
[Saved-frame verification](../../artifacts/workflow/surface-frame-verification-attempt1/frames/report.json)
passes frames 276-278, with unchanged vertices/contact boundaries and only two
triangles replaced at frame 277. The subsequent
[cached sweep](../../artifacts/workflow/cached-geometry-sweep-attempt1/geometry/report.json)
passed all 751 frames in 1711.761 seconds with unchanged fields and motion.
The user has since removed the requirement to cut/remesh canonical water output:
the approved path preserves the grid and computes waterline intersections
separately. Its complete producer replay and public-run acceptance remain pending.

| Message / category | Why it occurs | Evidence and physical implication | Disposition |
|---|---|---|---|
| Capytaine `Mesh resolution` | The unchanged 17-panel hull has panels large relative to the solved wavelength. The first coefficient warning reports radius 3.678 m exceeding wavelength/8 = 3.424 m at 1.5 rad/s; surface warnings begin at 1.45 rad/s. | The log summarizes 1374 coefficient, 5010 radiation-field and 382 diffraction solves affected. This is a resolution warning, not a measured motion/surface error. Sparse panels and short wavelengths remain unqualified. | **Provisional runtime warning; further accuracy evaluation needed.** Preserve first warning and counts. Do not call the result mesh-converged or silently remesh the baseline. |
| Capytaine `Irregular frequencies` / lid recommendation | Capytaine flags frequencies above its estimated first irregular-frequency threshold for this body/representation. | 612 coefficient and 4242 radiation-field solves are flagged; first reported frequencies are 5.3 and 5.25 rad/s. The warning identifies risk of spurious BEM solutions, not proof that every flagged solve is singular or that radiation physics is absent. | **Provisional runtime warning; further numerical qualification needed.** A lid or different discretization is a separate evaluated change, not an automatic cleanup fix. |
| `Surface numerical approximation`: 42.7 mm maximum / 12.5 mm RMS projection difference | The selected causal finite-band radiation reconstruction does not reproduce the complete direct complex near-field transfer. | The independent [method review](../experts/numerical_methods_researcher/surface-reconstruction-review-2026-10-25.md) accepted cutoff 24 rad/s, spacing 0.025 rad/s, full-history lags and zero direct acceleration for this case. The signed-lag comparator itself depends on future velocity; its difference is a material sensitivity, not a true-error bound. | **Evaluated provisional approximation.** Keep the warning. Tiny frequency-step sensitivity does not erase this phase uncertainty or certify physical accuracy. |
| `calm-restoring-moment` and `calm-restoring-force` | At each saved pose, exact restoring in a flat calm-water reference is compared with fixed-mean linear `-K*q`. The 10% threshold is a reporting convention with weight-based denominator floors. | Moment warnings start at 9.68 s: 2013 samples, peak departure 183%. Force warnings start at 16.66 s: 1170 samples, peak departure 132%. These are substantial reference departures, not nonlinear dynamic-error percentages. The flat surface differs from the actual incoming sea. | **Physical-applicability warning; further evaluation needed for accuracy claims.** Continue the explicitly provisional workflow without changing motion/forcing. Do not describe the departures as harmless merely because linearization predicts limitations. |
| `calm-restoring-unavailable`: opening rim reaches water surface | The flat-reference hydrostatics diagnostic cannot evaluate its restoring comparison when that reference intersects the opening or loses supported contact. | 338 samples; first at 23.08 s with reported flat-reference rim clearance -0.00153448 m. This is neither a computed actual-sea flooding event nor proof that actual-sea rim/contact is safe. | **Unavailable diagnostic; retain the limitation.** Actual reconstructed-water contact is a different assessment. Do not fill missing restoring comparisons or infer a physical flooding cause for the export failure. |
| Exact-sideband field / prepared-body load gap | Surface fields are reconstructed from direct sideband solves, while the accepted body is driven by its unchanged interpolated excitation table and low-frequency bridge. | Logged maximum differences are [20.7354, 18.0728, 47.5506] N and [21.6547, 90.0598, 16.2102] Nm. This is a known force/field consistency approximation; internal agreement of the body with its prepared loads does not make the gap zero. | **Evaluated provisional discrepancy.** Report it without modifying accepted loads. Exact dynamically coupled force/field closure and physical accuracy remain unestablished. |
| MATLAB `Escaped character '\U' is not valid` | WEC-Sim `waveClass/listInfo` concatenates the Windows input path into the `fprintf` format string, so `C:\Users...` is interpreted as escapes. | Stack points to `listInfo` line 411; only the displayed wave-elevation filename is truncated to `C:`. The case uses an importer-compatibility elevation file and separately binds/audits the actual six-component force table. The run continued through body and field production. | **Fix display formatting; not the observed solve failure.** Use a literal string format for path display in the maintained implementation. Preserve the actual path in logs; do not infer a failed file read from this display warning. |
| Capytaine `Exporting problem in already existing directory` | The producer creates its own fresh scratch hydrodynamics directory, writes the netCDF and then calls the legacy hydrostatics exporter. That exporter warns whenever its target directory already exists. | The logged path is the unique per-run OS-temp directory. The generic warning does not establish a collision with another case or an unexpected overwrite. | **Noisy lifecycle warning for this path.** Keep freshness/ownership checks; improve local presentation if needed. This warning is distinct from the previous attempt's cleanup/locked-file error. |
| Capytaine `Precomputing tabulation` | The Green-function implementation constructs its numerical lookup table; its source emits this progress message at WARNING level. | Actual work and startup cost are occurring. The message itself reports no failed tabulation or physical invalidity. | **Progress/info presentation.** It can be presented as initialization progress without hiding subsequent physical warnings. |
| WEC-Sim `Set water depth to 200m for visualization` | For infinite BEM depth, `waveClass.setWaterDepth` sets `deepWater=1` and a finite display depth. | The current workflow separately verifies infinite-depth BEM input and uses the external directional force table with nonlinear hydro/drag off. This message does not change that hydrodynamic model to a 200 m seabed. | **Expected informational message.** Label it visualization-only. Reassess its use if future physics paths consume the display parameter. |
| `surface has degenerate triangles`, frame 277 / 22.16 s | Required surface clipping/export produced geometry rejected by its triangle-validity gate. The specific construction defect is under geometry-owner investigation. | The traceback is `export_surface -> validate_surface_topology -> _mesh_arrays`. Attempt 2 preserves completed fields and failure geometry in the printed OS-temp recovery directory. This failure is distinct from the flat-reference diagnostic starting at 23.08 s; no causal link is established. | **Historical fatal in the superseded clipped export.** The guarded repair passed all 751 cached frames. Canonical output now preserves the uncut producer grid; validate its independent contact checks and diagnostics before acceptance. Retain this failed evidence and its original gates. |

The attempt stopped before all 751 surface frames and final surface summaries.
Missing final wet-set/contact or total-minus-incident summaries are incomplete
evidence, not a pass. A surface mesh that validates does not itself establish
nonlinear contact, flooding, breaking, slamming, convergence or accuracy.

The earlier degenerate clipped-export geometry is historical failure evidence.
The immediate task is validating the approved uncut export and its performance.
The physical/numerical warnings remain qualifications on how
the completed result may be interpreted; they are not suppressed and are not
new automatic rejection thresholds. Native dynamics remain held. Logging is
separate from physics: the CLI now reserves `sea.log`, `result.log` and
`video.log`, while this archived attempt log remains the evidence for this triage.

Source checks: [restoring diagnostic](../../src/topwater/wecsim_backend.py),
[surface producer](../../src/topwater/_capytaine_surface.py),
[geometry validator](../../src/topwater/_surface_export.py), and
[native WEC input/audits](../../src/topwater/data/wecsim/run_workflow_case.m).
The installed Capytaine checks/export/tabulation implementations and maintained
WEC-Sim `waveClass.m` were read directly; their executed identities belong to
the attempt's [frozen manifest](../../artifacts/workflow/canonical-pipeline-attempt2/manifest.json).

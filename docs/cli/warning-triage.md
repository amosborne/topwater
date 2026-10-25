# Canonical workflow warning triage

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
passes frames 276â€“278, with unchanged vertices/contact boundaries and only two
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

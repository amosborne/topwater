# Force-kernel and inertia packaging checkpoint

These are qualification results, not changes to production physics, the held
WEC-Sim revision, or canonical motion and surface data.

## Excitation-kernel qualification remains incomplete

The original actual WEC-Sim force comparison failed its initial-calm gate.
An independent replay matched WEC-Sim's convolution to about 1.1e-10 N,
isolating the finite supplied kernel rather than demonstrating a convolution bug.
The full-history fix was included in that run.

The genuine coefficient study uses three available headings, a 12 rad/s upper
band, a 0.02 s lag step, and unchanged reference complex excitation coefficients.
The numerical admission limits are 0.001 for maximum complex response error
divided by each mode's reference scale, and 0.0005 for successive-window and
frequency-refinement changes. These finite-probe limits do not certify the
complete directional sea or physical accuracy.

The single 160 s half-window control passed the reference-error limit in all six
modes. Its change from the refined 80 s control exceeded the stability limit in
pitch: 0.000795018. Other modes passed that comparison. No additional 160 s
frequency-refinement check was performed, and overall qualification did not pass.
Evidence: `artifacts/workflow/native-kernel-t160-attempt1/units/report.json`.

One prospectively fixed symmetric Tukey lag window, alpha 0.25, was then evaluated
on cached kernels. At the refined 80 s half-window, surge reference error was
0.00106943, exceeding 0.001. The 40-to-80 s change exceeded 0.0005 for surge
(0.00141962) and pitch (0.00055174). Frequency refinement passed. The candidate
therefore failed admission. This is a window on the kernel's lag coordinate;
it is not a wave or force startup ramp. No taper search, gain correction or
threshold adjustment followed.
Evidence: `artifacts/workflow/native-kernel-tukey-attempt1/units/report.json`.

Both jobs terminated and retained their frozen-input and process checks. The
next method requires independent review; neither result qualifies new actual-case
forces, motion or a runtime-pin change. Tracking: Topwater issue #35.

The independently reviewed finite-sum replay subsequently passed four small
scalar witnesses and read the exact stored H5 time axes. Its first real-data
attempt rejected the affine-grid shortcut: the conservative axis/interpolation
bound exceeded its prospectively fixed `1e-12` mode-scale guard on the first
40 s case. It stopped before full-array replay or target interpretation. This
does not measure an actual replay discrepancy or establish the surge-error cause.
No threshold was relaxed. The method reviewer is assessing whether a tighter
bound or explicit correction is justified. Evidence:
`artifacts/workflow/native-kernel-commuted-attempt1/replay/report.json`.

The independent reviewer proposed one revised source-only diagnostic: a fixed
centered reference lattice, a Taylor bound retaining interpolation's signed
barycentre cancellation, and an explicit output-phase correction using the
already cached kernel. Both numerical gates remain unchanged. The WEC/BEMIO
expert owns a new V2 source proposal; it requires independent review before any
allocation. This cached correction does not yet establish a cheap evaluator for
a new longer window without a cached kernel. No such control is authorized.

## Inertia library packaging remains incomplete

The supported MATLAB API authored and reopened the intended active rigid-body
products-of-inertia order correction. Direct R2020b export nevertheless changed
53 unrelated archive entries in both the untouched control and corrected
candidate. The unchanged strict auditor rejected both. Neither export was
promoted. Evidence: `artifacts/workflow/inertia-direct-export-attempt1/report.json`.

A separately reviewed private generator attempts to extract only that API-authored
property into the original archive, preserving every other uncompressed payload
and archive metadata. Its first execution stopped before output because the XML
text did not match the assumed lexical representation. The failed attempt and
source remain preserved. Independently reviewed V2 preserves the original XML
entity encoding and appends only the API-verified suffix. Its second private
execution passed: byte-identical control, exactly one intended candidate change,
zero unexplained changes, original R2020b metadata and all other uncompressed
entry bytes and archive metadata unchanged. Candidate SHA256:
`7b0c063c8c4db2dced7e42257002989d69cdbb7a789a21f0a044ca26c793d828`.
Only the target entry's compressed size changed, from 2854 to 2863 bytes.
The single-property audit has not been relaxed.
Evidence: `artifacts/workflow/inertia-parameter-delta-attempt1/output.log`.
Successful extraction and strict audits:
`artifacts/workflow/inertia-parameter-delta-attempt2/report.json`.

The portable, generic extraction helper then reproduced the identical candidate
and passed the same strict audits. A fresh hidden MATLAB process loaded both
absolute library paths, verified old/new expressions, R2020b metadata, lock and
clean state, and unchanged bytes after closing without saving.

The exact committed fork head `67e7c93f7039881a4af9b6f44a4f6232c029429c`
passed eight layout controls, three native tests (twelve tiny simulations), and
the public minimal before/after example (two further simulations). The matrix
tests support their `1e-7` tolerance; no new roundoff-level maxima are claimed
because the full-matrix fixture's raw report was temporary. The minimal example's
raw report is retained. Independent measured review accepted this scope.

Fork PR #9 now contains the correction and is ready for review at that tested
head. No merge, production pin or upstream submission occurred. R2020b execution,
Drag/Flex, actual CIC/dynamic residual forces and whole-boat accuracy remain
outside this component proof. Evidence:
`artifacts/workflow/inertia-portable-delta-attempt1/`,
`artifacts/workflow/inertia-reopen-attempt1/`, and
`artifacts/workflow/inertia-pr-head-attempt1/`.

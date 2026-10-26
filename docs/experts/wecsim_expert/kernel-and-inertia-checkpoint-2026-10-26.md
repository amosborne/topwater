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

The corrected V3 shortcut also failed its unchanged bound guard in the first
40 s case. Its four small scalar witnesses passed, but no full response replay
or attribution followed. The shortcut is stopped; no further bound revisions
or threshold changes are planned. Evidence:
`artifacts/workflow/native-kernel-commuted-attempt3/replay/report.json`.

## Fixed T320 controls passed the declared numerical admission

The next control returned to the official finite quadrature. A dense fine-grid
temporary would occupy 2.289 GiB, with several such temporaries potentially
exceeding available RAM. An explicitly named private implementation processes
256 independent time rows at once, preserving the complete frequency vector,
integrand, interpolation, reduction order and normalization. Its private handoff
copy changes only its function name and selected kernel-call name. Production
and fork implementations remain unchanged.

After independent source review, the 259-row/9601-frequency whole-versus-blocked
fixture passed. Cached T40/.005 and T80/.0025 then reproduced every consumer
kernel and complete complex-response value with exactly zero difference.
Evidence: `artifacts/workflow/native-kernel-blocked-replay-attempt2/units/report.json`.
Attempt1 was a preflight abort before launch on an inherited logging-style check;
its evidence is preserved, and the logging-only correction passed static review.

A separately allocated T320 pair then passed all predeclared gates in all modes:

| Criterion | Worst result divided by mode peak reference scale | Limit |
|---|---:|---:|
| T320/.0025 reference error | 0.0002700813 | 0.001 |
| T320/.00125 reference error | 0.0002837216 | 0.001 |
| T320 coarse minus saved T160 coarse | 0.0003769891 | 0.0005 |
| T320 fine minus coarse | 0.0000136456 | 0.0005 |

The independent method reviewer accepted this finite-probe admission. These are
mode-peak-normalized complex errors, not pointwise relative-error guarantees.
Coverage remains the original query grid, genuine 0/45/90 degree headings plus
seam, cutoff 12 rad/s, lag step .02 s and unchanged quadratic low-frequency
bridge. Full-heading actual-force/startup and motion qualification remain open.
Both jobs terminated with unchanged frozen identities; root verified fresh
MATLAB/launcher/producer absence and released their leases. A working-set sample
during T320 was about 1.1 GiB; no allocator peak-memory guarantee is claimed.
Evidence: `artifacts/workflow/native-kernel-t320-attempt1/units/report.json`.

## Genuine full-heading coefficient candidate is preserved exactly

The separately allocated full-band producer completed all 3606 genuine solves
through the existing API, with excitation .1..12 rad/s at all 25 periodic
headings. Its core producer and physical geometry are unchanged. Actual solve
time was 14.6 s after Green-function tabulation. The new operator manifest is
`f1f788f91cf42ae7d42eb86291f05636f9a1bdfc993fb028ae4d58d5ae28fdd4`.

The independently authored/reviewed preservation comparator passed. Every
recorded numerical comparison is finite and exactly equal, with zero differences:
all old low-band FK/SC/H values at all headings, full radiation A/B/Ainf, original
three-heading full-band excitation, physical/frame metadata, oriented 17-panel
geometry, zero-frequency limit and seam. Complete differences are retained.
The non-author reviewer accepted only the coefficient-candidate handoff.
All vendor `.libs` DLLs were included in this run's external freeze.

Actual warnings remain: 3180 mesh-resolution and 2040 irregular-frequency
warnings. They are not waived or physical-accuracy certification. First mesh
warning: omega 1.5, panel radius 3.678 m versus wavelength/8 3.424 m. First
irregular-frequency warning: omega 5.3, with a lid recommendation. Physical
accuracy remains provisional under the existing issues.

Evidence: `artifacts/workflow/native-fullband-producer-attempt1/` and
`artifacts/workflow/native-fullband-preservation-attempt1/comparison/report.json`.
Both jobs are terminal, frozen-input/process checked and released. Next gates are
actual official H5 sign/normalization/full-radiation readback, all-heading
numerical admission and actual forces/startup/application/motion. No production
settings, dependency pin or canonical physics output changed.

## Directional aggregation and dispatch have a tested draft PR

[Fork PR #10](https://github.com/amosborne/WEC-Sim/pull/10) is a draft at
`eb722540da07c3cba737193f43712bf395ff15fa`, explicitly dependent on unmerged
PR #8. Four group/native-bus tests and eleven history controls pass on this exact
head. The retained identical-input example measures maximum discrepancy 18 on
the PR #8-only baseline and 0 on the paired candidate. Native applied-total
assertions pass at `1e-10`, with measured fresh/regular/restored branch counts
1/2/1 and rewired-clone rejection without deletion. Raw full-bus measurements
were temporary; the claim is the retained test pass at its asserted tolerance.

Attempts 1/2/3 exposed fixture workspace and count-query errors; they remain
preserved. The final fixture uses WEC-Sim's base workspace with owned-variable
restoration, literal block discovery, and regex only on actual names. No force
or branch assertion was weakened. Strict library audits permit exactly one
header change; fresh API reopening preserves its bytes and metadata.

The independent reviewer accepted publishing the draft. Whole-boat/startup
qualification remains pending, and the binary needs reconciliation/retesting
after PR #8/#9 integration. No merge, pin change or upstream action occurred.
Evidence: `artifacts/workflow/imported-group-head-attempt4/`.

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

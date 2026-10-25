# First isolated fork regressions

The sidecar-handle and literal-path candidates were tested independently against
official v7.1.0 (`a23d2993807410fe78cc63546c1ba9d81912024f`). This does not
promote the combined reference branch or change Topwater's held dependency pin.

Coordinator evidence: `artifacts/workflow/fork-regressions-attempt3/` contains
the frozen manifest, stock/fixed logs and JSON test results, imported MAT structs,
comparison report and terminal receipt. Two sequential fresh hidden MATLAB
R2025b processes ran on Windows. The existing shared engine was left alone.
Exact launchers/producers were independently verified terminal, frozen hashes
matched, idle was read back and the numerical lease released.

| Suite | Stock release | Independent fixed candidate |
|---|---|---|
| CAPYTAINE handles | Five expected failures; malformed Hydrostatics control passes | Six pass |
| Literal path display | Char and string path checks fail; other-wave control passes | Three pass |

The handle suite covers repeated sphere and two-body RM3 imports, malformed KH
and Hydrostatics files, a missing later-body sidecar, preservation/usability of
an unrelated open handle, and Windows rename/delete after import. The display
suite covers drive, percent-containing and UNC strings, supplied as char and
string, with warning checks and an unchanged no-wave display control.

All 32 numeric fields in each imported sphere/RM3 struct have exactly equal
values before and after the handle correction. This uses SciPy
`assert_array_equal`, with matching struct keys; it is not a byte-level identity
claim. The archived report's `numeric_fields_bitexact` key is an overly strong
label for that value comparison. The existing sphere warning about removing
three low-frequency NaN entries occurs in stock and fixed imports and is not a
new warning from the lifecycle correction.

Attempt1 reproduced stock leaks but its coordinator wrapper failed while
appending into an untyped empty struct. Attempt2 reached display discovery and
exposed a helper accidentally named `test_paths`; MATLAB treated it as a test.
The wrapper was corrected and the helper renamed `fixture_paths` before the
successful new immutable attempt3. Failed evidence remains under its original
attempt directories. Neither correction changes production parsing or physics.

An independent reviewer approved source and harness acceptance, with these
limits: malformed-input tests accept any thrown import error; only the selected
fixtures and Windows/R2025b were exercised. No dynamics simulation, full WEC-Sim
test suite, older MATLAB release or additional OS was tested. Upstream submission
and integration into the reference branch remain separate actions.

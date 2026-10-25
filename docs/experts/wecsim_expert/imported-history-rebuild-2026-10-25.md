# Isolated imported-excitation history correction

Own-fork issue [#1](https://github.com/amosborne/WEC-Sim/issues/1) is rebuilt
independently from official v7.1.0 (`a23d299`). This correction concerns the
first imported wave's excitation preprocessing. It does not change directional
aggregation, applied-force dispatch, inertia, radiation or the held runtime pin.

The complete supplied elevation record contributes to
`F(t) = dt * sum_tau K(tau) * eta(t-tau)`. The original
`min(kernelTime):dt:max(kernelTime)` sample grid is preserved; source-time
queries remove the centering assumption of `conv(...,'same')` for asymmetric,
one-sided and off-grid lag supports.

Beyond the complete file, elevation is assumed zero and a named
`WECsim:ImportedExcitationRecordSupport` warning reports missing margins once
per preparation call. Short records are permitted; malformed/nonfinite data
inside the record remain errors. This deliberately avoids imposing Topwater's
calm-start/full-support gate on generic WEC-Sim users.

Ramping follows the existing absolute-time law at source time, exactly once.
Zero ramp preserves supplied negative-time data. Nonzero ramp retains the
existing negative-time cosine behavior, with no new pre-zero clamp. This is
compatibility behavior, not a recommended physical startup model. The original
seven-argument method continues to accept already-ramped samples without a
second ramp. These policies differ from the historical combined `2f43f54` fix.

Independent review approved candidate `15a38c68b176cd121317a70d0571b0d4460370b1`
and its scalar interval-interpolation/nested-quadrature oracle before execution.
The helper subclass exposes protected excitation methods and stubs only unrelated
mass/radiation preparation. Key witnesses invoke actual `hydroForcePre` and
`waveClass.setup` with synthetic files through stock-supported APIs; stock
failures establish numerical cropping/lag defects rather than merely rejecting
new optional arguments.

Root's serial `artifacts/workflow/history-regressions-attempt1/` contains the
frozen source manifest, commands, logs, stock/fixed JSON results, report and
terminal receipt. Two fresh hidden Windows/MATLAB R2025b processes ran
sequentially. All eleven candidate tests passed; stock failed ten and passed
the regularCIC control, with no incomplete tests. The selected clean stock
witnesses were `testFullFileHistorySixChannels` and
`testAsymmetricAndOneSidedLags`. Static checks of the new test/helper and root
wrapper also passed.

Coverage includes six distinct channels, impulses before/inside/after output
times, two time steps, asymmetric/one-sided/off-grid support, source-time ramps,
negative history, nonzero simulation start, already-ramped legacy input, partial
records/warnings, invalid/nonfinite inputs and an analytic regular-wave control.
Exact processes were independently verified terminal, frozen hashes matched,
idle was read back and the execution lease released. The shared MATLAB engine
was left alone. Documentation-only amendments preserve the tested production
and test blobs.

This is synthetic excitation-preprocessing evidence on Windows/R2025b. No boat
motion, multi-import behavior, full WEC-Sim suite, other MATLAB/OS combination,
large-record memory/runtime or physical accuracy is qualified. Public example
artifacts remain unchanged. No merge, dependency promotion or upstream submission
is part of this checkpoint.

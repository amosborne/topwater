# Native force resource checkpoint

The imported-wave force transition remains unqualified for the complete boat.
The coefficient controls and one physical component have now been measured;
neither result substitutes for total-force, startup, actual-bus or motion checks.

## Completed checks

`artifacts/workflow/native-kernel-shared-phase-replay-attempt1` reproduced all
retained-heading kernels and complex responses exactly using shared trigonometric
blocks. A separately allocated `native-kernel-allheading-admission-attempt1`
then admitted all25 headings under the original per-mode scales and budgets:
reference error .001, window change .0005, frequency refinement .0005.
The actual `units/report.json` has `all25_numerical_admission_pass: true`, SHA256
`a503f95ca9ce57a49067f2310e9dec4a661499d4f58b84d1dc2112702e666a51`.
These remain finite-probe numerical checks, not physical-accuracy guarantees.

`artifacts/workflow/native-force-component-preflight-attempt1/probe` exercised
the held995 native method on the first active physical positive branch, carrier16,
with unchanged amplitude, phase, heading and origin. It used3001 output samples
over0..60 s and35003 padded elevation samples with the admitted T320 kernel.
All identity, coverage, finite-output and declared resource guards passed.
Common-support elevation error against the old input was4.86e-17 m.

The actual native call took4.7419052 s. Windows reported a lifetime whole-process
peak working set of4,002,521,088 bytes (3.73 GiB); this is not an isolated
allocation measurement. The largest of six component relative RMS differences
against the separate offline spectral oracle was6.7063e-5 (0.00671%).
The component's initial heave force is not subject to the total-force calm gate:
physical directional components can cancel. No total startup conclusion follows.

Receipts are `probe/report.json`, `probe/native/resource.json`, and the job's
`terminal.json`. Root verified terminal processes, unchanged frozen identities
and lease release. Independent method review accepted the limited evidence.
The existing unrelated shared MATLAB engine was left untouched.

## Next bounded proposal

A simple736-component extrapolation is about58 minutes per force pass; that full
runtime is unmeasured. Own-fork [issue11](../history/wecsim-retired-workflows-2026-10-27.md#issue-11)
tracks a private standard linear-convolution candidate, based on PR8/PR10's
`eb722540`. The isolated candidate is
`artifacts/dependencies/wecsim-imported-convolution`; the unchanged detached
control is `artifacts/dependencies/wecsim-convolution-baseline`.

For eligible uniform clocks, sample the same source-history interpolant once,
apply the existing source-clock ramp, and use zero-padded FFT linear convolution.
Unsupported or ambiguous clocks use bounded-row direct evaluation. Preserve
the candidate's original lag origin, negative-time ramp behavior, legacy
already-ramped inputs, warnings, direction interpolation and physical weights.
No Topwater-applied forces or frequency-response shortcut is introduced.

Before any full-case run, require independent source review, unchanged11 history
and4 group/native-bus controls, independent scalar witnesses, and all3001x6 actual
component samples replayed against the retained native result. Both maximum
error/reference peak and time-weighted RMS error/reference RMS must be <=1e-10.
Zero channels require explicit justified handling, without denominator floors or
physical clipping. Measure runtime and fresh-process lifetime memory separately.
Source authors run no numerical jobs; root owns serial execution.

No candidate has yet passed these replay gates. No production settings, runtime
pin, canonical physics artifact or fork merge changed. Complete force/startup/
application/motion verification and user-held dependency promotion remain open.
Mesh-resolution and irregular-frequency warnings also remain unqualified.

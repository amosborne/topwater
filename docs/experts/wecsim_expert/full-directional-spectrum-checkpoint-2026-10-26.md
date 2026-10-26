# FullDir spectral correction checkpoint, 2026-10-26

User requires one own-fork issue/PR. Issue12 covers explicit phase preservation
and consistent directional weights; root owns folding13. Current combined
source topic is `wecsim-spectral-spectrum`, `fix/full-directional-spectrum`.
Combined scientific source was tested at8ee8dcb199f21a08cf81e61bfab141d0b0624364.
Old separate topics and their evidence remain preserved. Runtime995, reference
and public artifacts are unchanged.

First cheap Windows/R2025b allocation passed its declared acceptance:
`artifacts/workflow/native-spectral-small-attempt1/acceptance.json`.
The baseline reproduced both failures: phase discrepancy6.2628106249rad,
physical field discrepancy0.7862334196m and weight-only six-channel numeric
discrepancy0.4107283522. Phase-only150a66df26914308f1a284b187f5f4e2480e647c
preserved phases exactly, with field discrepancy1.11e-16m. Weight-only
f72c5ccd01846e9cb0258258e1419b4be6a03fa9 had weight-only discrepancy1.11e-15
using actual prepared phases. All4 phase and3 first-order weight controls passed.
Per-arm JSON, raw arrays, logs and `terminal.json` are retained. Root verified
frozen identities/exact process termination and released the lease.

These first results establish synthetic first-order API/native-function defects
and separate corrections. Combined validation subsequently passed below.
The packet's exact low-frequency/interpolation/seam handoff, actual forces,
startup and motion remain separate.

The combined delta adds one zero-oscillatory/nonzero-MD control against the
independent sum of physical squared amplitudes times dimensional MD coefficients,
with changed phases to prove invariance. The combined suite keeps all original
assertions in one actual native-bus invocation, now also requiring physical-force
agreement. Nine controls total:4 phase,3 first-order widths,1 MD and1 bus.
Root's independently source-reviewed combined allocation passed all9 controls
and the portable actual native-block example on Windows/R2025b in172.4749907s.
Supplied phases were exact; field maximum difference1.1102230246251565e-16m;
native-function/physical-force maximum difference1.3322676295501878e-15. Actual
native bus/function difference was0, and bus/physical sum difference
1.3322676295501878e-15. MD acceptance is its tested1e-10 assertions, not an
invented raw MD roundoff maximum. Evidence:

- [Combined acceptance](../../../artifacts/workflow/native-spectral-combined-attempt1/check/acceptance.json).
- [Portable example report](../../../artifacts/workflow/native-spectral-combined-attempt1/check/example/report.json)
  and neighboring `arrays.mat`.
- [Transcript](../../../artifacts/workflow/native-spectral-combined-attempt1/check/matlab.log),
  neighboring `test-results.mat`, and outer `terminal.json`/frozen manifest.

Root verified frozen identities/exact process termination and released the lease.
Independent measured reviews are separate from this factual documentation update.

The transcript's coarse-direction coverage warning comes from the generic
coverage heuristic on this deliberately sparse synthetic table. Queries are at
explicit H nodes; no off-grid or physical mesh accuracy is established. Native
infinite-water setup also reports200m for visualization; deep-water remains
enabled, without substituting finite-depth hydrodynamics. Logs retain these
messages. No warning policy or physical gate was weakened.

Qualification stops at this synthetic component/native-bus proof: full packet
low-frequency/seam/interpolation handoff, actual forces/startup/motion, runtime
and broader MATLAB compatibility remain pending. No author scientific imports,
jobs, commits, pushes or runtime/upstream promotion occurred.

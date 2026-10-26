# FullDir spectral correction checkpoint, 2026-10-26

User requires one own-fork issue/PR. Issue12 covers explicit phase preservation
and consistent directional weights; root owns folding13. Current combined
source topic is `wecsim-spectral-spectrum`, `fix/full-directional-spectrum`,
head748032696b3234bfefa58ce4609283da418df2bb before the pending test/doc delta.
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

These establish synthetic first-order API/native-function defects and separate
corrections, not combined-head, native-bus, mean-drift or physical accuracy.
The packet's exact low-frequency/interpolation/seam handoff, actual forces,
startup and motion remain separate.

New source-only delta adds one zero-oscillatory/nonzero-MD control against the
independent sum of physical squared amplitudes times dimensional MD coefficients,
with changed phases to prove invariance. The combined suite keeps all original
assertions in one actual native-bus invocation, now also requiring physical-force
agreement. Nine controls total:4 phase,3 first-order widths,1 MD and1 bus.
Independent review and root's serial allocation precede execution. No author
scientific imports, jobs, commits, pushes or promotion occurred.

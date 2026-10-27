# Issue 31: selected-case numerical acceptance

The accepted simple-boat, 60 s WEC-Sim trajectory passes the user's **10 mm
maximum surface-height error including numerical reference uncertainty**, against
the chosen fixed-mean linearized model. This is selected-case engineering
qualification, not a nonlinear or experimental accuracy claim.

| Quantity | Result |
| --- | --- |
| Largest pointwise difference plus reference uncertainty | 5.78442494698 mm |
| Largest actual runtime/reference difference | 0.34118305099 mm |
| Largest time-weighted RMS runtime/reference difference | 0.08964879218 mm |
| Checked surface | 5,929 targets, 751 times, full 0--60 s |
| Native field evaluation | 5,899 exterior queries; 384 targets use explicit exterior traces |
| Physical geometry/contact/startup checks | All 751 frames passed |

The combined quantity adds each location's measured difference and numerical
uncertainty. It is not the sum of unrelated global component maxima, nor an
estimate of how far the linear model departs from reality. Future-continuation
sensitivity is the largest individual uncertainty component.

The actual-difference and combined-budget maxima both occur at **44.48 s**,
target3283, native mean-hull XY **(2.757895, 0.821091) m**. This continued
location samples its admissible exterior trace at **(2.823489, 0.864523) m**.
At that point the runtime surface is 0.441289854 m and the reference is
0.441631037 m: a signed difference of -0.341183 mm. Its reference uncertainty
is **5.443242 mm**, giving the reported **5.784425 mm** combined quantity.
An independently source-reviewed saved-array reporting pass verified both
argmax locations separately; no surface, motion or uncertainty was recomputed.

## What changed and what the reference means

Runtime surface reconstruction now preserves the complete complex radiation
response, including phase. It applies that response offline to the accepted
saved velocity. The WEC-Sim motion, forcing, mean hull and chosen hydrodynamic
operator were not modified or reintegrated. Incident and diffraction fields
were verified as common; old-total and per-term identity residuals were about
1e-15 m. A moving hull does not feed back into this linear surface calculation.

Validation uses a separately declared C1-Hermite interpolation through saved
velocity/acceleration. It separates the analytically known instantaneous
acceleration term in the reference only; runtime still uses intact full H.
Padding, time step, input interpolation, startup, future-prefix dependence,
frequency spacing and cutoff controls are charged pointwise. The known input
tail uses an independently checked integral enclosure, including frequencies
above Nyquist. Failed causal-cut and loose-tail attempts remain preserved.

One explicit condition remains: at every query and mode, the regular response
above 48 rad/s stays below twice its largest calibrated magnitude from dense
36--48 plus 64/96 rad/s data, with known-32/prior-probe witnesses. This assumes
no unseen larger resonance or faster growth. It is a declared engineering
envelope, not a proof of infinite-frequency behavior or model resolution.
It must accompany this acceptance; no generic 10 mm guarantee for other boats,
seas or settings follows from this case.

The retained wetted-area warning reports 18.27% symmetric difference at 44.64 s
(threshold exceeded from 29.04 to 53.12 s). This flags fixed-mean linear-model
applicability, not failed numerical reconstruction or a quantified nonlinear
motion error. Keep it and inherited solve warnings with their provenance.
Issues #32--34 and broader physical/model accuracy remain separate.

## Evidence and publication checkpoint

The focused reconstruction suite can be run with:

```powershell
uv run --no-sync python -m pytest -q -p no:cacheprovider tests/test_offline_radiation.py tests/test_surface_reconstruction.py
```

The exact final integration cohort is listed in the tracked record's
`tests.cohort` array; pass those node IDs to the same pytest command. It covers
publication/rollback, logs, saved-surface error recovery, archive contracts,
CLI integration, real short video encoding, visual layers and shared plot ticks.
On this Windows host the executed checks used isolated workspace TEMP/TMP and
MPLCONFIGDIR, Agg rendering, one numerical thread and disabled plugin autoload.
The physical qualification itself consumes the retained operator and saved
trajectory; this small regression command does not repeat that campaign.

The [tracked acceptance record](issue31-acceptance-evidence-2026-10-27.json)
contains coverage, numerical budgets, unchanged-data identities, canonical file
hashes, integration-test results and terminal receipts. The larger raw arrays
and execution logs referenced below remain preserved local diagnostic evidence.

- [Final method acceptance](../../../artifacts/workflow/issue31-final-method-review1/acceptance.json),
  SHA-256 `778e72871f57f90b01ba7940809cfb2561e1dd2350aaa5c16f4ec4552933c318`.
- [Independent full saved-array audit](../../../artifacts/workflow/issue31-independent-check-attempt1/results/report.json),
  SHA-256 `55575ea13936296500202d9234e75e88f45c5c6e03fe74a94ab9a618ffe3eb50`.
  It recomputed every target/time budget and verified motion, hull, connectivity
  and every packed water coordinate against the validated field.
- [Completed physical replay and private artifacts](../../../artifacts/workflow/issue31-recovery-attempt1/results/report.json),
  SHA-256 `b1e8f511d5d3833ab816bb3ac054d47c2f4f84b61da13cf42c2c6fd5d3071ae3`.
  The method receipt binds source context, operator, policy, terminal/guard
  proofs and unchanged-motion Result/plot identities.
- [Published canonical video log](../../../examples/simple-boat/video.log),
  SHA-256 `9b26e817b980dd266f9743d99f85737bfda705d74c4dc10b26970d4e2cadd2f3`.
  Decoded verification passed: 751 frames at 12.5 fps, 800 by 800, through source
  time 60 s. Encoded duration is 60.08 s because the final sample is included.

The canonical Result, plot and logs are published, and the video has been
regenerated. Readiness was promoted after execution cleanup; only the Boolean
changed and numerical-function ASTs were verified unchanged. No image inspection
was performed by this reviewer. Final integration checks subsequently passed
50 tests in 24.20 s with no skips, and their independent receipt accepted the
coverage and unchanged numerical code. Repository checkpoint and issue closure
are recorded by the coordinator in the handoff and GitHub issue.

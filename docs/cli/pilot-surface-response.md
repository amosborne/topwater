# Motion-driven surface diagnostic

Developer-only `scripts/pilot_surface_response.py` verifies a selected saved
trajectory's reconstructed surface. It is not another dynamics solver and
does not publish canonical result/video files. Use the central serial claim.

```powershell
uv run --group wecsim python scripts/pilot_surface_response.py --physical-case CASE.json --native-mat NATIVE.mat --motion MOTION.mat --train TRAIN.npz --output NEW_DIRECTORY
```

## Inputs and outputs

Use the frozen physical mesh/mass/hydrostatics, exact accepted native motion
(world/modal velocities and accelerations), and unchanged signed-kappa train.
Output is a fresh diagnostic directory containing retained native/field data,
sampled time histories and a report. Existing inputs are read only. No body
simulation, forcing adjustment or render is performed.

## Assumptions and hard-coded choices

Capytaine3.0.0 on the same mean hull. Approximately5800radiation solves and one
combined diffraction solve per active train frequency; a600s budget is checked
between operations. Compare cutoffs12/24rad/s and spacings0.05/0.025rad/s over
the actual60s trajectory. Include near/far exterior points, explicit trace-offset
sensitivity, and exact incident evaluation at target points. Trace continuation
is an export approximation; no nonlinear contact or uniform error bound is
claimed. An independently evaluated infinity acceleration coefficient is a
diagnostic candidate; baseline convolution does not silently adopt a fitted term.

## Fatal errors and warnings

Wrong source identities/shapes, nonfinite data, invalid zero-prehistory,
unsupported interior field evaluation, failed numerical solves or exceeded
checked budget are fatal. Capytaine resolution/irregular-frequency warnings
remain relevant limitations. Reports measure numerical sensitivity, signed
direct-field versus prepared-load differences, and initial total elevation.
They do not tune amplitudes/phases or require the total-minus-incident delta
to be small.

## Detection limits

Pointwise pilot results are not full-domain error bounds, complete initial
fluid-velocity proofs, nonlinear moving-waterline physics or mesh/time
convergence. The common result publication gate stays closed until the selected
case reconstruction and export policy have been reviewed and implemented.

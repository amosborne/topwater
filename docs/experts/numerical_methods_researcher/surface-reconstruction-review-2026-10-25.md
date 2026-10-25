# Independent review: selected-case total-surface reconstruction

Date: 2026-10-25. Reviewer: `/root/startup_domain_policy`
(`numerical_methods_researcher`). Coordinator accepted the bounded policy below.
This review inspected source and saved reports only. The reviewer authored no
solver or physical correction and ran no tests, imports, BEM solves or body jobs.
The coordinator owns execution, artifact acceptance and integration.

## Decision and exact scope

Approve a **provisional computed linear total-surface reconstruction** for the
selected stronger-buoy 60 s trajectory on the unchanged 17-panel mean wetted
hull. Use radiation cutoff 24 rad/s, frequency spacing 0.025 rad/s, causal kernel
lags through the entire requested history, saved native modal velocities and
zero prehistory. The instantaneous acceleration-to-elevation term is zero.
Use the exact signed-kappa incoming packet and combined-direction diffraction
at its actual sideband frequencies. Keep the accepted body motion and prepared
forcing unchanged. No further refinement campaign is required by this decision.

This is approval of a declared numerical approximation for this case, **not**
physical-accuracy, spatial/time convergence, full-complex equivalence or a
universal surface-error certification. It does not qualify other boats, seas or
durations, authorize native dynamics, or mark historical diagnostic reports
`total_surface_ready=true` retrospectively.

## Method and evidence

Capytaine 3.0.0 uses unit generalized displacement with exp(-i omega t), so
the elevation transfer per modal velocity is `Hv = i * eta_q / omega`.
Combining incident normal-velocity data before each diffraction source solve
is valid by linearity. The finite-band cosine inversion supplies a causal
approximation; it must not be described as reproducing the complete direct
complex transfer.

- [Dense transfer diagnostic](../../../artifacts/workflow/refinement-checkpoint-attempt1/refinement/report.json):
  the earlier far-point worst per-mode recovery discrepancy was about 3.35%.
  Its normalization is a per-mode maximum across sampled points/frequencies,
  not a pointwise bound or a total-surface error percentage.
- [Actual-motion response diagnostic](../../../artifacts/workflow/surface-response-attempt1/field/report.json):
  0.05 to 0.025 rad/s spacing changed sampled radiation elevation by at most
  0.00000326012 m; cutoff 12 to 24 rad/s changed it by at most 0.000753226 m.
  Near-waterline full-complex recovery discrepancies remain about 17–30.4%.
  Shared-Green evaluation agreed with the public evaluator to 1.72e-16.
  The independently evaluated infinity acceleration term was exactly zero.
- [Signed-lag projection diagnostic](../../../artifacts/workflow/surface-projection-attempt1/report.json):
  direct finite-band full-complex inversion versus the causal cosine response
  differed by up to 0.0426650 m, with sampled RMS up to 0.0125110 m. The maximum
  persists before 55 s. The raw inverse also had a future-velocity contribution
  up to 0.0772652 m, so it is not a causally valid replacement reference.

These centimetre-scale phase/causality sensitivities are material and must stay
visible. The much smaller spacing/cutoff changes do not establish micrometre
or submillimetre accuracy. No sign/factor error or omitted wave contribution
was found in the selected construction. No fitted endpoint term is approved.

The endpoint normalization follows `D_a = -lim(eta_q / omega^2)` because
`a = -omega^2 q`. At regular exterior z=0 points, the infinite-frequency Green
function is Rankine minus reflected Rankine and vanishes; this explains the
zero acceleration term. Finite-frequency growth alone cannot identify a
physical derivative feedthrough. See the official [endpoint API](https://capytaine.org/stable/user_manual/post_pro.html#pressure-velocity-free-surface-elevation)
and [Green-function limit](https://capytaine.org/stable/theory_manual/inf_depth_green_function.html#singularities-extraction).

## Geometry, force and startup limitations

The linear field belongs to the exterior of the mean hull. A moving hull can
expose part of the mean-interior footprint. Continue the disturbance from an
explicit admissible exterior trace there, evaluating incident water at the
actual target location. Never silently sample interior Green functions.
This is an export-geometry approximation, not nonlinear contact physics.
The two tested exterior offsets, 0.025 and 0.05 m, changed sampled continued
total elevation by 0.00188311 and 0.00304654 m at the two targets. These are
sensitivities, not bounds. Sharp chines and waterline events preclude a blanket
second-order continuation-error claim without additional evidence.

Exact-sideband field loads differ from unchanged prepared loads by sampled
maxima [20.7354, 18.0728, 47.5506] N and [21.6547, 90.0598, 16.2102] Nm.
Retain this discrepancy; do not quietly alter forcing or claim exact dynamic
force/field closure. The surface is a reconstruction driven by the saved motion.

Initial radiated elevation was at floating-point zero. The seven sampled total
initial elevations were submillimetre, but these samples do not prove full-hull
total-field calm or fluid-velocity bounds. Ordinary complete-artifact gates
remain: initial-state checks, finite synchronized data, correct frames, valid
waterline/export topology, hull exclusion and complete output coverage. No full
artifact or physical startup acceptance follows from this methods review alone.

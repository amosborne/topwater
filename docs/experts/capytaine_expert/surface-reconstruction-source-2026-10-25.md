# CAP-SURFACE-001: retained-field source and bounded pilot

Status: source plus root-executed diagnostic evidence. Root reports the seven
analytical tests PASS in3.30s and the sparse field pilot PASS in66.125s;137solves,
three exterior points. Combined-source/individual-source errors were3.523657e-19m
for elevation and4.54747e-13 for wrench. Root froze and released the numerical
lease; this author launched no jobs. Evidence is under
`artifacts/workflow/surface-pilot-attempt1`. These checks do not qualify a total
transient field or high-frequency mesh accuracy.

The producer in `src/topwater/_capytaine_surface.py` accepts frozen physical
mean-wet triangles, mean CG, full dry mass matrix and exact setup hydrostatics.
It retains Capytaine 3.0.0 indirect-source results while exporting the original
radiation/import coefficient policy (0.1:0.1:12 rad/s and three import headings)
and the separate full-heading excitation policy (0.1:0.1:3.4 rad/s, periodic
25-heading table). No independent mesh change, mean-pose change or force fit.
The three-heading native import and full-heading applied-force table remain
separate; the incident train's load H0 bridge is owned by WEC setup.

The Capytaine API is verified from installed 3.0.0 source and the
[official postprocessing documentation](https://capytaine.org/stable/user_manual/post_pro.html#pressure-velocity-free-surface-elevation).
`compute_free_surface_elevation` needs the retained indirect source distribution;
the old coefficient-only exports used `keep_details=False` and cannot supply it.

## Field convention

Every harmonic field is `Re(A exp(-i omega t))` in the fixed native mean axes,
about the mean CG. For each actual nonzero train sideband, sum the weighted
incident Neumann conditions across all carriers and both signed-kappa branches.
Linearity permits one scattered-field solve per sideband. Negative kappa uses
positive omega with heading shifted by pi; no negative-omega outgoing Green
function is used. Add the exact incident spectral elevation once. Return the
integrated direct FK+diffraction spectral wrench separately for comparison with
the existing applied body loads; never replace those loads silently.

Capytaine radiation problems impose `-i omega` times modal displacement. Thus
`H_eta,v = i H_eta,q / omega`. The candidate causal cosine inversion subtracts an
explicit instantaneous velocity term first. Any further derivative, pole,
low-frequency or high-frequency contribution needs analytical treatment and
independent acceptance before inversion can be relied upon. Full complex
transfer recovery uses `integral h(tau) exp(+i omega tau) dtau`, not a check of
only its cosine/real part. Held-out frequencies are required.

Radiation convolution uses actual six fixed-mean-CG native radiation-input
velocities, zero prehistory and the complete requested lag horizon. It does not
inherit the body's 10 second CIC cutoff. Candidate kernels cannot be used by
`convolve_radiation` without explicit acceptance. `require_surface_readiness`
currently fails before body execution; no total-field-success path exists yet.

## Bounded first execution proposal

1. Run `tests/test_surface_reconstruction.py` first under root allocation. Its
   exponential-plus-instantaneous-term fixture has an analytical complex transfer
   and analytical transient sinusoidal response. It checks causal time ordering,
   spectral normalization and refusal of unaccepted/short kernels. The cosine
   inversion itself has an analytical finite-band tail/trapezoid bound and a
   wider-band/finer-grid refinement check; it does not assume a finite band is
   exact.
2. Run `scripts/pilot_surface_reconstruction.py` with explicit canonical physical
   JSON and accepted exact native MAT. It requires the unchanged 17 wet triangles,
   probes three exterior field points and 13 frequencies from .025 through96
   rad/s. Expected127 coefficient solves plus10 sideband-linearity solves;
   budget180seconds, checked between solves. No MATLAB or body integration.
3. The pilot stores raw complex endpoint data and verifies combined-source
   linearity against eight individual direction/branch solves. It deliberately
   records `total_surface_ready=false` even when these pilot checks pass.
4. Independent review determines endpoint/singular terms, frequency coverage,
   sufficiently dense field transfer and the complex-recovery tolerance. Only
   then may a dense kernel pilot and total-field exporter be implemented/accepted.

The sparse pilot is not quadrature coverage, near-waterline convergence or a
complete total-field result. Incident-plus-scattered startup calm still needs
checking, as does mismatch between new exact-sideband field pressures and the
accepted table-interpolated body loads. Linear mean-body reconstruction does not
enforce nonlinear instantaneous contact; mesh clipping alone does not fix that.

Both field entrypoints reject points inside or on the mean-waterline convex
footprint before evaluating Green functions. For a concave hull this deliberately
over-rejects; it does not admit any interior points. Near-hull export may later
continue an explicit exterior trace to the geometric contact edge, while keeping
the incident field exact at the target. That needs selected-case offset
sensitivity and does not carry an O(a^2) accuracy claim for a sharp-chine hull.

## Endpoint and recovery follow-up source

`scripts/analyze_surface_pilot.py` defaults to stored-data endpoint/scaling and
finite-band complex-recovery diagnostics. It labels sparse interpolation
explicitly and never accepts a production surface. Its `--refine` mode proposes
direct same-mesh transfer evaluation at .05rad/s spacing through24rad/s, smaller
low-frequency probes, four high-frequency probes and six held-out frequencies,
at the same three exterior points. About3000tiny source solves,300s checked
budget; no MATLAB/body job. No follow-up analysis or solves have been executed
by this author. Root must allocate even the stored-data analysis.

The driver reports a diagnostic real-coefficient high-frequency decomposition
`D+i*A/omega+B/omega^2`, full real/imaginary held-out errors, cutoff12/24 and
spacing.1/.05 sensitivity over the full60s lag horizon. These are inputs to a
selected-case finite-band approximation decision, not a universal-convergence
requirement or an automatic endpoint policy. The producer also accepts an
explicit job-temp Green-function cache directory, with original GF numerical
settings preserved, to avoid tabulation rebuild between instances in one job.

## Actual-response checkpoint source (not executed)

Root reports the dense three-point refinement PASS83.188s, with finite-band
complex recovery materially better at .05rad/s than .1rad/s over60s. This is
provisional transfer evidence, not an elevation-error bound. No production
acceptance is inferred from those percentages.

`scripts/pilot_surface_response.py` now proposes one selected frozen-motion
checkpoint: actual native `velocity_world` and `acceleration_world`, cutoffs12/24,
spacing.05/.025, three direct near/far points plus two newly exposed-strip targets
at .025/.05m exterior trace offsets. It compares exact sideband field wrench
against unchanged prepared loads and reports incident+diffraction startup.
The bounded source uses no new high band beyond existing96rad/s probes and has
a600s checked budget. All computations remain coordinator-allocated.

The exact symbolic Capytaine infinity result provides a *candidate*
`D_acceleration=-lim(eta_q/omega²)`. Its mathematical sign follows directly from
`acceleration=-omega²*q`; it is not fitted from observed high-frequency growth.
The pilot reports it and its actual-acceleration contribution independently;
the baseline remains zero direct-acceleration term pending review. High-frequency
growth alone does not establish a physical derivative term. Shared S-matrix
evaluation of six retained source columns has an explicit public-field-evaluator
equivalence check in the pilot before any production use.

Exterior trace mapping continues disturbance fields only; incident elevation is
evaluated exactly at target XY. The resulting selected-case offset sensitivity
is an export-approximation diagnostic, not an O(a²) sharp-chine error guarantee.
The readiness gate remains closed while the source checkpoint is committed.

# CAP-SURFACE-001: retained-field source and bounded pilot

## Current failure-recovery source checkpoint

Root's actual canonical attempt1 reached frame201/751, then exited2 after
1989.285s with no public result. A Windows KH.dat cleanup lock masked the
original geometry/metric exception, so the exact primary failure is unknown.
It is not evidence that flooding, topology or the wet-set gate was the cause.
No contact or applicability gate is relaxed by this source correction.

Independent review approved physical-output snapshot compaction and safe
invalid-time diagnostics. Root executed the current recovery/interface suite:
**209PASS,1optional skip,13.05s**, evidence under
`artifacts/workflow/recovery-interface-tests-attempt1`. This includes injected
primary-error/evidence preservation, geometry-only replay with physical field
producers forbidden, physical output free of unused cap vertices, and
missing/empty/nonnumeric cached-time causes preserved. Successful scratch removal
and failure retention are also covered by the WEC owner's interface tests.
No full new artifact is accepted; actual pipeline attempt2 remains pending.
Source/tests are frozen for root's checkpoint/retry. These documents alone are
updated at this docs-only checkpoint.

Before any frame export, the producer now writes a private
`surface-recovery/completed-field.npz` beneath its OS-temp workflow directory.
It contains immutable completed total/incident fields, grid, actual full motion,
sampled poses, body/rim/cap identities, reference plane data and input policy
scalars. It is an internal disposable intermediate, never an example output.
The ordinary WEC scratch context deletes it on success. On a
`SurfaceExportError.recovery_directory`, that context retains its own workspace
and preserves the same primary exception; secondary cleanup failures cannot
replace the physical/geometry error.

Per-frame errors immediately print their original traceback plus stage, frame
index and physical time. They save a small `failure-frame.npz` containing the
untrimmed graph, posed closed hull, physical faces, rim/caps, incident field,
reference plane and clipped mesh when available. `failure.json` records the
original error and exact cache location. Evidence-writing errors are secondary.
If completed-cache validation fails before a frame/time is established, the
diagnostic time is explicitly absent; error handling never re-reads the bad
time axis and masks the original validation exception.

`replay_surface_cache` runs only clipping, contact/startup checks and metrics on
the completed data. It invokes no Capytaine, MATLAB, body integration or wave
synthesis. A root-frozen replay may optionally write the standard Result NPZ
without changing its motion arrays. The public physical hull snapshot is
compacted/remapped to remove vertices used only by artificial caps; the cached
full closed hull, rim/cap indices and CG frame remain unchanged. Example
root-allocated replay command:

```powershell
uv run python -c "from topwater._capytaine_surface import replay_surface_cache; replay_surface_cache(r'C:\Users\alex\AppData\Local\Temp\topwater-wecsim-CASE\surface-recovery\completed-field.npz', diagnostic_directory=r'artifacts\workflow\surface-replay-attempt1', output_path=r'artifacts\workflow\surface-replay-result.npz')"
```

Both diagnostic/output destinations must be fresh. The original private cache
is never overwritten or implicitly deleted by replay. Freeze and claim normal
serial execution before replay; it is a numerical geometry job, not read-only
activation. This author has executed no cache replay or test; root's suite above
executed the focused analytical/geometry replay controls. An actual cached-field
replay from a new pipeline failure has not yet been needed or accepted.

## Current production source checkpoint

Source implementation is complete and coordinator-selected provisional policy
is enabled. Focused current-source tests passed as recorded above. Actual
generic WEC-Sim execution and the complete total-surface artifact still require
the next root execution and acceptance after the masked first attempt.
No job was launched by this author. Earlier source-only/closed-gate paragraphs
below are historical checkpoints superseded by this current section.

`reconstruct_total_surface(hydro,setup,motion,surface_time_s)` returns packed
vertices/faces and their per-frame offset arrays. It uses the actual native
six radiation-input velocities, with the body trajectory and applied loads
unchanged. Selected numerical policy is cutoff24rad/s,
`delta_omega=min(.025,.98*pi/T)`, full requested lag horizon, zero instantaneous
velocity/acceleration terms, and explicit causal cosine projection. The exact
symbolic infinity candidate evaluated to zero in root's actual-response pilot,
as expected from Rankine/reflected-Rankine cancellation at regular exterior
points on z=0.

Root's selected frozen-motion near/far checkpoint PASS226.234s: frequency-step
.05-to-.025 response difference at most3.26012um, cutoff12-to24 difference at
most.753227mm, trace-offset .025-to-.05 differences1.88311/3.04654mm at two
exposed-strip targets. Shared/public field evaluator and direct-source
linearity checks passed. These are selected-case sensitivities, not physical
accuracy bounds or universal mesh qualification.

Full complex transfer recovery has **not passed** near the hull. The independent
stored-data signed-lag-versus-causal-cosine diagnostic differed by at most42.7mm
and12.5mm RMS at tested points; its noncausal alternative is not a reference
solution or error bound. Independent researcher/coordinator accepted the causal
projection as an explicit provisional linear export approximation. Production
prints this limitation once; it must not claim that fine frequency spacing
removed this phase approximation. No coefficients/tails were fitted to motion.

Disturbance evaluation batches six retained source columns through one standard
Capytaine S matrix per frequency. The public-equivalence check is in root's
pilot evidence. Kernel inversion/convolution is point-batched to bound memory.
Repeated Capytaine mesh/irregular-frequency warnings retain the first warning
and print a count/limitation summary instead of flooding each sideband.

The total field is sampled on an upward .1875m P1 grid. Its native patch includes
the hull trajectory plus one-beam margin and the exact default fixed WORLD-view
corners, transformed to native axes, at the renderer's Slerp surface poses; one
grid cell pads the union. Native-axis extent alone is not used as a substitute
for world-camera coverage. Custom larger framing can fail explicit coverage.
Disturbance targets inside/near the mean footprint use a nearest exterior trace
at `min(.025m,.01*beam)`. Incident waves remain exact at the original targets.
At every saved frame the total graph is trimmed against the actual posed full
convex hull by the independently tested exporter. Physical opening-rim/cap
contact fails as unsupported flooding. Counts/connectivity may vary by frame.

The initial total-P1 gate bounds local represented elevation on all exported
triangles touching the initial hull XY box, including diffraction and radiation.
It uses the existing epsilon*.sum(amplitudes) elevation criterion. Shared
full-hull **incident** eta/velocity bounds still run separately. No total
disturbed-fluid velocity bound has been established, and output says so.

Reported total-minus-incident metrics use the same represented P1 incident
surface interpolated at new cut vertices. They include projected-area/time mean,
RMS, resolved peak/world-location/time, incident RMS and relative RMS where the
denominator is nondegenerate, for the full patch and a one-beam band. The band
is classified by triangle-centroid horizontal distance to the actual clipped
total-surface/hull **waterline**, not the full hull's horizontal outline.
Incidence-one outer-patch edges are excluded using the known native rectangular
domain bounds. Remaining contact segments must lie on support planes of the
actual posed hull. Missing/invalid contact fails rather than substituting an
outline. This centroid spatial classification is stated explicitly.

Wet-set warning uses exact geometric areas of the represented P1 graph, not
sampled point proxies: intersect each physical hull triangle with water-cell
vertical prisms, clip by current water and transported flat-mean wet planes,
then sum `A_current+A_mean-2*A_common`. The denominator is the initial **flat
mean-equilibrium hydrostatic** wetted area, not the tiny-tail total-P1 t=0 wet
area. Shared prism-edge ownership handles vertical facets without double count.
The partition must reproduce the known initial area; missing coverage fails.
Above10% symmetric-difference area, report first/worst/last time. This is a
fixed-mean-operator applicability warning, not nonlinear motion accuracy.

Root evidence: `artifacts/workflow/surface-response-attempt1/field/report.json`,
`artifacts/workflow/surface-projection-attempt1/report.json`, prior sparse/dense
pilot records. New analytical/geometry integration test source covers actual
area symmetric difference, transported material wet set, vertical partition
ownership, missing coverage, P1 moments, and packed calm geometry at an oblique
global heading. Those additions have not yet run at this checkpoint.

## Earlier source and diagnostic checkpoints

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

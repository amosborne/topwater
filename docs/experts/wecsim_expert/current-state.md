# Current state: WEC-Sim / BEMIO expert

Updated 2026-10-28 UTC. Authority: [handoff](../../checkpoint-handoff.md),
[physical contract](../shared-contract.md), [coordination](../coordination.md).
Root owns acceptance, integration and serial execution. Consult the
[execution ledger](../execution-state.json) for live/terminal status; this page
confers no execution permission. All eight example files are accepted, including
the final ordinary CLI video. The final producer is terminal with matching frozen
inputs and an empty guard Job. Root independently confirmed idle state, no active
job and no numerical lease after release.

## Current 130-second packet integration

Implementation checkpoint `e416ca8` passed 170 focused Python tests in 56.93 s.
The directory solve contract takes positional dwell, with rest/ramp defaults of
10/5 s: dwell 100 means total time 130 s, native time step 0.01 s and surface
cadence 0.08 s. Boat/sea carrier arrays stay unchanged. Phases refer directly to
dwell midpoint 65 s in the existing negative-time physical convention.
`_wecsim_setup` calls the utility with negative phase and conjugates positive
phasors; each sideband keeps original heading and its own dispersion. Aggregate
1 cm envelope-exceedance budget is allocated proportional to carrier amplitude.
Original carrier operating eligibility is retained; generated sidebands above
it are warned about and retain genuine H coverage, without mesh-accuracy claims.

The approved DC reference shift keeps all computational BEM/CG/K geometry fixed.
Native initial heave `-d` supplies the full stiffness heave-column coupling;
physical motion and surface receive `+d` once. Surface geometry runs in the
computational frame; startup tests physical dynamic elevation plus `d`.
Native net-load/startup, applied-force, full inertia, CIC and equation-of-motion
audits remain. No Topwater force history or fork change was introduced.

The ordinary solve completed all native MATLAB assertions before saving
`motion.mat`: phase, startup net load, initial pose/rest/memory, complete clocks,
applied forces, full inertia, Newton–Euler balance, added-mass delay and full CIC
replay. The invocation then failed only during restoration of the root-owned
engine's initial cache path. It took about 13.3 minutes; no public Result was
published by that failed invocation. Native motion contains 13001 samples at
0.01 s, with phase difference `1.7583157152500917e-14 m` and preserved static
sea level `d = 0.012167229735326152 m`.

Recovery2 passed a real production restoration control with initial cache and
code-generation directories beside TMP and normalized path identity checks.
Its owned engine quit before field work. The saved native arrays passed their
original audit tolerances. Rebuilt hydrodynamics matched all 14 numeric
`native.mat` arrays and every retained netCDF variable at relative tolerance
`1e-12` and absolute tolerance `1e-9`.

The complete new incident, diffraction and radiation field is cached before
geometry processing. Cache SHA-256 is
`9f3870fb9f3e91deb9b82ea83b25f62ed24e6acb43ad11dfd0d1acf3157660fb`.
Root stopped recovery2 to avoid its wall-time ceiling. Geometry-only replay
passed all 1626 frames and published the accepted Result/PNG pair at
2026-10-28 05:04:00 UTC. No body or field recomputation occurred. Ordinary CLI
rendering passed at 05:21:46 UTC: all 1626 decoded frames, 800 by 800 pixels,
12.5 fps, source time 0–130 s and encoded duration 130.08 s, with zero warnings.
CLI elapsed time was 303.813 s; render core time was 273.828 s. The earlier
60-second surface accuracy bound does not qualify this new packet.

The first render attempt exceeded the old fixed 240-second verification budget.
The CLI now uses `max(240, waterframes)`; 25 CLI tests passed in 3.99 s and the
final ordinary render passed. This is an operational budget correction, with no
new science. The original ordinary solve remains a failed cleanup invocation;
the successful Result log truthfully identifies the private continuation.

Actual represented wetted-area departure exceeds 10% from 14 to 116.4 s, with
a peak of 25.13% at 43.36 s. All reconstructed-surface opening/contact gates
passed. The total-minus-incident disturbance is 6.83618 mm RMS over the full
exterior patch and 8.03576 mm in the one-beam band; the reported incident-relative
ratios are 1.77782% and 2.09892%. Its resolved absolute peak is 51.6636 mm at
24.48 s. These are disturbance measures, not numerical or physical errors.
The peak-location log's Z is computational: 0.4378107785 m corresponds to world
Z 0.4499780082 m after adding the stored DC offset. Published geometry and delta
metrics are correct; this is a nonblocking diagnostic-coordinate presentation issue.

Geometry postprocessing took 934.221 s: wetted-set diagnostics 617.859 s,
delta metrics 266.193 s, opening checks 27.780 s and waterline 21.731 s.
This is measured postprocessing cost, not another body/field solve.

The packet has 384 carriers (368 active), 128 positive frequency nodes through
0.984615 Hz and 272 unique headings. Its summed certified origin envelope
exceedance is 0.00942218 m; full-hull incident startup also passed. Generated
sidebands above the original 0.533459 Hz operating limit are retained with a
warning: their amplitude L1 contribution is 0.107014 m, or 0.5749% of total L1.
These checks do not establish mesh accuracy or nonlinear applicability.

Evidence: [tracked run record](../packet-integration-run-2026-10-27.json),
[focused tests](../../../artifacts/workflow/packet-production-python-tests-attempt1/results/report.json),
[failed ordinary invocation](../../../artifacts/workflow/simple-boat-packet-run1/results/report.json),
[native/operator recovery](../../../artifacts/workflow/packet-native-recovery2/results/report.json),
and [geometry replay](../../../artifacts/workflow/packet-field-replay1/plan.md).

## Adopted reference and workflow

**The default pin is6b774228e1f6b31dee473e626204b0ef65e5ec38.** The user merged
PR6/7/9/14 into the release-based `topwater/reference`. The default checkout is
`artifacts/dependencies/topwater-wec-sim`, clean on that branch; no path override
was used for actual default-path preflight.

Adapter `9539192bb733663584483f33ef2fb3de273d6fcf` is the historical adoption
checkpoint underlying the reference evidence below. The active implementation
also includes the new packet/DC integration described above. The old source
worktree and detached reference checkout remain for provenance, not runtime use.

WEC-Sim computes forces and motion through its EXISTING public
`bodyClass(hydroStruct)` spectral input. Official BEMIO import/normalization and
`rebuildHydroStruct` prepare ordinary excitation fields; consumed Ainf/CIC and
restoring inputs remain unchanged. Unused raw A/B/T are omitted rather than
fabricated. No custom coefficient-view/accessor or BEMIO writer extension remains.
Topwater supplies waves/coefficient data and reconstructs the total surface.

**The accepted eight-file result now contains the new 130-second packet case.**
Boat/sea inputs are unchanged. The earlier reference adoption checkpoint
`7e1b854` and its cleanup evidence below are historical. **Cleanup actions are
complete**. Independent current-state audit confirms source/output identity,
remote/local retirement and recovery preservation.
Only fork branches `main`0753b2e and `topwater/reference`6b remain; exact995 is
preserved by the verified tag `archive/topwater-legacy-9952634`, not an active
legacy branch. [Completion/recovery record](fork-cleanup-plan-2026-10-27.md).

## Historical accepted reference evidence

- Tiny public-input discrimination passed three low/off-grid/seam and missing/
  empty-unused-A/B cases; consumed CIC/Ainf/mass/restoring were exact.
  [Acceptance](../../../artifacts/workflow/standard-input-review-v1/acceptance.json).
  Narrowed PR14 also passed ten native regressions, including the real force block.
- Exact6b/953 passed38 Python and21 MATLAB cases (13 validators +8 handoff).
  [Acceptance0cffbbc3](../../../artifacts/workflow/reference-handoff-review-v1/acceptance.json).
- Full canonical3001x6 forces/startup/phase, once-only application, rest/negative-
  history certificate, full CIC replay, clocks and motion passed.
  [Acceptance286f3a82](../../../artifacts/workflow/reference-adoption-motion-review-v1/acceptance.json)
  binds [case3 reportec4d8d8d](../../../artifacts/workflow/reference-adoption-case3/report.json).
  Maximum position difference59.758 micrometres and SO3 angle0.00156985deg meet
  1mm/.01deg. Controlled runtime1212.407s includes simulation745.503s and native
  force replay305.104s; it is not all integration time. Failed prior attempts remain.
- Fresh TOTAL incident/diffraction/radiation field, standard Result/PNG, full
  ordinary CLI render and exact-eight audit passed from the saved NEW motion.
  No body rerun, coefficient campaign or old-field substitution.
  [Private acceptance4622b301](../../../artifacts/workflow/reference-adoption-private-review-v1/acceptance.json)
  binds final report2104ca1f,3001 motion/751 surface samples and751 decoded frames.
  Same operatorf1f788f9,17 panels and centroid quadrature were retained.
- The active Topwater checkout passed38 focused tests in2.69s and real default preflight: clean6b
  checkout, no override, hidden fresh R2025b, all three required licenses, zero
  loaded models and verified owned-server shutdown.
  [Default-path acceptance3712ab33](../../../artifacts/workflow/reference-adoption-main-review-v1/acceptance.json).
  This readiness check performed no further body integration or public solve.
- [Promotion report](../../../artifacts/workflow/reference-adoption-promotion1/report.json)
  verifies all eight canonical files byte-for-byte against the accepted candidate:
  Result2b615b56, PNG275ac4c5, videoe5fce0e4. Physical inputs are unchanged and the
  verified rollback copy of all eight original files is retained.

- [Cleanup execution](../../../artifacts/workflow/wecsim-cleanup-execution1/) verifies
  native archival8/10, deletion of obsolete issues1/2/11/13 and nine remote heads,
  fork-only Actions disable, plus removal of11 old worktrees/13 trial refs.
  Merged PR6/7/9/14 and closed issues3/4/5/12 remain. The shared-store bundle,
  shallow boundary, six dirty and24 ignored snapshots passed recovery checks.
  Production6b, stocka23, the clean common parent/linked validation and source953
  remain; original failed studies are preserved.

## Limits and preserved history

This is implementation/numerical compatibility acceptance, not physical accuracy.
The prior selected-case surface qualification remains historical evidence;
the new 130-second packet is not yet numerically qualified to that bound.
Issues #32–34 separately track mesh/coefficient accuracy, force-field consistency
and fixed-mean applicability. The native coarse-direction warning is a coverage
heuristic for the repacked nodes; matching forces does not establish genuine
angular coefficient resolution. Mesh/irregular-frequency, field-phase and fixed-
mean applicability warnings are not error bounds or blanket waivers.

Earlier e24 optional-view qualification is historical evidence, not the current
architecture: [motion](spectral-motion-qualification-2026-10-27.md),
[end-to-end](spectral-end-to-end-qualification-2026-10-27.md).
[Inertia proof](inertia-pr9-qualification-2026-10-26.md) covers active Rigid Hydro,
not Drag/Flex or general accuracy. [Engine lifecycle](workflow-engine-lifecycle-2026-10-25.md)
and [retired workflows](../history/wecsim-retired-workflows-2026-10-27.md) retain
original evidence/policies. Do not resume retired studies from historical plans.
Generic fixes remain upstream candidates; no upstream issue/PR/contact is
authorized without explicit user permission.

[Independent completion audit](../../../artifacts/workflow/reference-adoption-completion-review-v1/audit.json)
records the current source, API and recovery checks.

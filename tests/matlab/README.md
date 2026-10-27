# Workflow coefficient validation

Run from the repository root in MATLAB R2025b:

```matlab
results = runtests('tests/matlab/workflowHydrodynamicsTest.m');
assertSuccess(results);
```

No BEM solve, model, shared engine, external fork or binary fixture is needed.
The suite calls the same coefficient guard used before the phase comparison and
before H5 export. It preserves valid infinite water depth and checks finite
coefficients, compatible shapes and their frequency/heading/time reference grids.

The minimal pre/post witness is
`testNaNRealExcitationRejectedWithLegacyWitness`. It injects one NaN into the
same high-frequency excitation coefficient in an imported/raw pair. The exact
old residual-max predicate accepts the pair because `max` omits NaNs; the new
production guard must throw `Topwater:NonfiniteHydrodynamics` and name
`imported.ex_re`. A finite control with `h.h = Inf` must still pass. The
imaginary-channel witness repeats the check with the opposite phase sign.

Independent source review and coordinator-owned Windows/R2025b checks passed:
all eleven hydrodynamics tests and both engine-setting tests pass. The same
setting tests fail against the old guard. Evidence is
`artifacts/workflow/current-path-checks-attempt2/`. The witnesses evaluate the
old finite-data predicate directly; they do
not claim a pre-fix full solve or that the canonical result contains nonfinite
coefficients. Physics, transforms and normalization are unchanged.

The later native-excitation handoff adds two passing controls for a zero
excitation-kernel frequency endpoint. Raw coefficient and radiation grids remain
strictly positive; this is not a blanket exception for nonfinite coefficients.

## Native excitation coefficient handoff

`workflowHydrodynamicHandoffTest.m` uses the official BEMIO functions from the
coordinator-frozen `TOPWATER_WECSIM_PATH` checkout, falling back to the repository
dependency checkout. It generates small synthetic six-mode netCDF/MAT/sidecar
inputs in owned temporary folders; it performs no BEM or Simscape/body solve.

```matlab
results = runtests('tests/matlab/workflowHydrodynamicHandoffTest.m');
assertSuccess(results);
```

The synthetic radiation table contains nonzero damping only at a frequency
above the common excitation raw grid. Its preserved kernel must still match a
full-band BEMIO transform exactly after writing and reading H5. Other controls
check the single native-to-BEMIO phase conversion against an independent kernel
formula, the explicit zero-frequency bridge, genuine component sums and periodic
seam, native physical identity, finite records and refusal to invent raw
frequencies that were not recorded. These are exchange/transform tests, not
physical hydrodynamic accuracy or native applied-force equivalence evidence.

Production `prepare_workflow_hydrodynamics(prepared)` requires existing
`prepared.n`, `netcdf_path`, `native_mat_path` and explicit excitation IRF
half-width/time-step/frequency-step controls. It returns a normalized BEMIO
structure and writes the private `canonical.h5`. The full radiation kernel is
computed first from the original wider band. Standard H5 raw tables use only
genuine common frequencies and full periodic headings; the excitation kernel
has its own bridged quadrature grid. No applied force history is synthesized.

Independent source review and coordinator-owned Windows/MATLAB R2025b verification
passed all 19 tests: 13 coefficient guards and six official BEMIO/H5 handoff tests.
The frozen fork source was `25b58bda96ab0272cbfb4b0755539d5bb8ca227e`;
the job's manifest records the exact helper, test and dependency hashes. Evidence:
[report](../../artifacts/workflow/native-handoff-tests-attempt1/report.json),
[MATLAB log](../../artifacts/workflow/native-handoff-tests-attempt1/matlab.log)
and [terminal record](../../artifacts/workflow/native-handoff-tests-attempt1/terminal.json).
The process exited successfully, with frozen inputs verified and its execution
lease released. These synthetic controls demonstrate the exchange contract,
including exact retention of the full-band radiation kernel and direct added
mass through H5 readback; they do not qualify the production caller or WEC force
and motion equivalence. Its three integration review blockers remain separate.

Native WEC force-only comparison and downstream applied-sum checks remain
mandatory before replacing the baseline.
Finite excitation lag, quadrature and record padding can change numerical
forces; source exchange tests do not prove that those differences are small.

## Standard public spectral input

The earlier measured results above bind their original helper/test identities.
The selected path uses unextended WEC-Sim `bodyClass(hydroStruct)` and ordinary
coefficient fields. The maintained default checkout is
`artifacts/dependencies/topwater-wec-sim`, pinned to6b774228; an explicit
`TOPWATER_WECSIM_PATH` must select that same reviewed revision. The old995 runtime
lacks the FullDir corrections and is rejected before physical preparation.

Eight handoff cases retain the six legacy imported-history/H5 controls and add
standard spectral-struct preservation plus invalid-input checks. The spectral
branch performs official `readCAPYTAINE`, normalization, full-band `radiationIRF`
and `rebuildHydroStruct`; then supplies derived excitation samples through existing
`simulation_parameters.w/direction`, `hydro_coeffs.excitation.re/im` and
`hydro_coeffs.mean_drift`. The input table includes the exact native component
nodes and periodic boundaries, retaining the old bilinear/quadratic policy.

Ainf, CIC kernel/time, restoring and physical properties stay unchanged. Unused
raw frequency A/B and raw-period metadata are omitted from this mode-specific
struct instead of relabeled or invented. The genuine producer/raw BEMIO arrays
remain unchanged outside that runtime transport. No excitation IRF or H5 is
needed on the spectral branch. Legacy imported-elevation preparation still
computes a genuine excitation IRF and uses the unchanged official H5 writer.

Run the command above for `workflowHydrodynamicHandoffTest.m`, together with
`workflowHydrodynamicsTest.m` for the existing thirteen coefficient guards.
Python coverage: `uv run pytest tests/test_wecsim_spectral.py tests/test_wecsim_backend.py`.
The21 native validator/handoff cases and38 focused Python checks passed at6b+953.
The canonical full force/application/rest/CIC/motion and fresh eight-file delivery
also passed independent review. These are selected-case compatibility checks, not
mesh/time convergence or general physical-accuracy claims.

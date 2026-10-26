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

## Spectral component view (source-only additions)

The preceding 19-test result binds its original frozen helper/test identities.
It does not qualify these new spectral additions. Three additional handoff
tests are authored in `workflowHydrodynamicHandoffTest.m`, bringing that file to
nine cases; the coefficient-guard file retains its existing thirteen controls.
Authors have not executed the new tests.

Use the coordinator-frozen fork containing the reviewed optional-excitation-IRF
writer support via `TOPWATER_WECSIM_PATH`; held995's writer does not support the
new no-IRF case. Do not silently change the runtime pin. The existing legacy
tests cover all-three-IRF-fields-present. New tests check spectral omission with
every other official raw/radiation/physical readback field exactly preserved,
partial IRF rejection before an existing output is truncated, and invalid
derived-view finite/zero-drift/H0 failures before H5 creation.

The selected spectral handoff accepts `prepared.component_view` and returns
`[h,nativeHydro]`. It builds no excitation IRF or dummy kernel; the official
writer/reader exchange raw coefficients and full radiation, and the view is
attached only to the returned public body struct. This is a coefficient/IO gate,
not native harmonic force, startup, applied-bus or motion acceptance.

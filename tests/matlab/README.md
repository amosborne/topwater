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

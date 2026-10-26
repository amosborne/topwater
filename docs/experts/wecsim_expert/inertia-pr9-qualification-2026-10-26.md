# Active Rigid Hydro inertia: tested PR9 checkpoint

Updated2026-10-26 UTC. Coordinator-authorized source, packaging and serial checks;
no merge, runtime repin or upstream submission. Our fork
[PR9](https://github.com/amosborne/WEC-Sim/pull/9), “Correct rigid hydrodynamic body
inertia-product order,” is ready for review at
[`67e7c93f7039881a4af9b6f44a4f6232c029429c`](https://github.com/amosborne/WEC-Sim/commit/67e7c93f7039881a4af9b6f44a4f6232c029429c),
branch `fix/rigid-hydro-inertia-order`, based on combined `25b58bda`.
The production pin remains `9952634f`. Older failed studies remain preserved.

The one physics correction appends `([3 2 1])` to the active Rigid Hydro File
Solid's `adjustedInertiaProducts` expression. Public products `[Ixy Ixz Iyz]`
enter File Solid's `[Iyz Izx Ixy]` order without a sign change. Public properties,
production mass adjustment and residual-added-mass code remain unchanged.
Drag/Flex are excluded.

## Packaging and API evidence

Supported API authoring/export/reopen produced the authoritative old/new property
values. Its untouched control and candidate each contained53 unrelated serialized
changes; both strict audits failed and were not waived. The narrow extractor
transplanted only the verified property text onto the original R2020b archive,
preserving lexical entities, all other uncompressed payloads and non-derived ZIP
metadata. Original library SHA256 is `0ef996e4ef79f1eef3bcf451719a39e397820937f3e4586a3e99bcf8c61beaeb`;
corrected candidate is `7b0c063c8c4db2dced7e42257002989d69cdbb7a789a21f0a044ca26c793d828`.

Local evidence below lives under ignored `artifacts/`; these are workspace paths,
not GitHub-hosted files:

| Check | Outcome | Evidence |
|---|---|---|
| Narrow extraction | Control byte-identical; candidate exactly1 intended expression,0 unexplained/added/removed entries; R2020b metadata | `artifacts/workflow/inertia-parameter-delta-attempt2/{control-audit.json,candidate-audit.json}` |
| Portable generic helper | Reproduced exact candidate SHA above; unchanged strict control0/candidate1 audits pass | `artifacts/workflow/inertia-portable-delta-attempt1/report.json` |
| Fresh API reopen | Exact absolute FileName, old/new expression, original lock on, Dirty off, R2020b; close without save preserves bytes/SHA | `artifacts/workflow/inertia-reopen-attempt1/reopen/report.json` |
| Exact PR head |8layout controls and3native tests pass;12tiny simulations in native setup | `artifacts/workflow/inertia-pr-head-attempt1/{acceptance.json,fork-head.json,test-results.mat}` |
| Public minimal example | Two additional simulations complete; retained raw before/after report | `artifacts/workflow/inertia-pr-head-attempt1/minimal-example/report.json` |

Each job also retains manifest, log and terminal evidence in its named directory.
Root independently verified frozen identities and exact process termination;
leases were released. Cap accepted the measured qualification. Static checks
preceded execution. The tested fork README's “pending” wording predates the run;
the PR body records the current result. No prose-only amendment invalidated the
tested head.

## What the native checks establish

The actual active File Solid is copied into an identity/rest/zero-CG/zero-gravity
fixture. Six independent unit forces/torques identify its full6x6 matrix from
direct t=0 acceleration. Both dry and production-adjusted bodies run in each
before/after fixture; the actual nonsymmetric residual is added algebraically.
The exact-head tests passed their `AbsTol=1e-7` matrix/acceleration assertions and
unchanged-products/residual assertions. Detailed full-matrix setup output was
temporary, so this checkpoint does not claim newly measured roundoff errors.
Earlier standalone roundoff values remain in
[the earlier proof](inertia-adjustment-proof-2026-10-25.md), with their own scope.

The retained public example applies the same1Nm roll torque to the same7kg body
with rotational tensor `[10 1.1 -1.3;1.1 12 1.7;-1.3 1.7 14] kg*m^2`:

| Initial angular acceleration, rad/s² | Release expression | Corrected expression |
|---|---:|---:|
| Roll |0.104094244641134 |0.102615384041695 |
| Pitch |-0.015746134614160 |-0.010944563702830 |
| Yaw |0.010903090436361 |0.010857554110644 |

These are retained exact-head example measurements; corrected values agree with
the independent `M*a=w` result within the tested tolerance. No dynamic CIC,
simulated residual-added-mass force, full boat, Drag/Flex or older-MATLAB runtime
qualification follows. Fixture construction uses R2025a+ physical-port APIs;
actual execution coverage is Windows/R2025b. An R2020b SLX label alone establishes
neither older runtime nor operating-system compatibility.

## Portable reproduction and checks

The fork ships `tests/helpers/prepareRigidHydroInertiaLibrary.m`,
`extractRigidHydroInertiaDelta.py` and the unchanged strict auditor. Use a separate
stock release library: API-author separate control/candidate exports; run the
generic extractor with that stock, exports, MAT receipt and author source;
strict-audit its control/candidate; freshly API-reopen the accepted result without
save. The generic helper computes runtime hashes and validates receipt paths and
expressions. The MAT receipt is execution provenance, not authenticated historical
byte evidence. Python with SciPy is needed for receipt parsing; archive/XML work
uses the standard library. Exact portable commands are in the tested fork
[README](https://github.com/amosborne/WEC-Sim/blob/67e7c93f7039881a4af9b6f44a4f6232c029429c/tests/README-rigid-hydro-inertia.md).

In a fresh owned MATLAB process, with `forkRoot` set to the exact proposed head:

```matlab
shape = runtests(fullfile(forkRoot,'tests','probeInitialSampleTest.m'));
assert(all([shape.Passed]) && ~any([shape.Incomplete]));
native = runtests(fullfile(forkRoot,'tests','rigidHydroInertiaTest.m'));
assert(all([native.Passed]) && ~any([native.Incomplete]));
addpath(fullfile(forkRoot,'tests','examples'));
rigidHydroInertiaExample(newPrivateOutputDirectory);
```

Requires licensed MATLAB, Simulink, Simscape and Simscape Multibody. Source/library
identities and outputs must be frozen and process trees verified terminal under
the coordinator's serial lease. Root owns further integration and acceptance;
merge or runtime-promotion work requires its separate authorization.

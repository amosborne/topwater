# WEC-Sim excitation handoff: component evidence

The replacement source passes its focused Python and MATLAB startup controls.
It is not yet qualified by an actual same-input force or motion comparison.
The held fork pin and canonical example outputs remain unchanged.

Topwater now prepares physical directional elevation records and coefficient
data. WEC-Sim computes excitation and motion. Both signed packet branches retain
their directions, amplitudes, phases, spatial origin and absolute time. The
production path supplies no precomputed wave-force history. An offline old-force
calculation remains available only for independent comparison.

## Reviewed source and checks

The private source manifest `native-excitation-comparison-v1/source-manifest.json`
has SHA256 `d22f5ad4d94a10ed86ad22c61899ce9da4d42eadfce03cad4a0d8bd447f8f184`.
It records the source review checkpoint; the later test-only teardown correction
is identified separately below. The Capytaine and separate BEMIO experts reviewed
the consequential boundaries independently of their authors.

| Check | Measured result | Evidence under `artifacts/workflow/` |
|---|---|---|
| Backend and private comparison gates | 36 Python tests pass, 5.12 s | `native-python-tests-attempt1/` |
| Original isolated startup controls | 2 MATLAB tests pass, with Windows fixture cleanup warnings | `native-startup-tests-attempt1/` |
| Startup controls after teardown correction | 2 MATLAB tests pass, no warnings | `native-startup-tests-attempt2/` |

The Python controls cover directional origin-phase reconstruction, padded times,
packet frequency support, backend boundaries, independent current-setup metadata,
oriented geometry/file identities and missing native output. A carrier inside the
coefficient band cannot conceal an unsupported packet sideband. Missing output
fails even when MATLAB returns normally, retaining any original MATLAB error.

The startup controls check the diagnostic carrier scale against an independent
formula and reject a noncalm initial force using the original limits. They do not
compute or apply an alternative force history. The cleanup fix only moves cwd
restoration ahead of temporary-directory deletion; test math is unchanged.
Corrected test SHA256:
`01f60c0c654e9d70c46e1aceadc8595954bdd1968c311bdc1303a4ece05ee9d2`.

Each job used a serial root allocation. Frozen identities, terminal process and
MATLAB launcher absence were independently checked before releasing its lease.
MATLAB ran in fresh hidden batch processes, without using the shared engine.

## Remaining gates

The new fork clone cleanup compares internal signal edges and ports before any
deletion. Its rewired-clone rejection and actual applied-bus tests are source
reviewed but not executed. Private supported-API library preparation succeeded,
but its strict whole-archive audit failed: 188 changed entries and two additions
beyond the intended one mask-header change. Evidence:
`rigid-dispatch-preparation-attempt1/`. That library is not accepted or promoted.
The failed artifact and original inputs remain preserved for diagnosis.

Retained coefficient reuse also needs independent historical lineage and actual
operator identity checks. The comparison manifest verifies consistency; a new
manifest alone does not prove provenance. A genuine dense FK/SC bank was found,
but has not been accepted as a drop-in replacement for current producer output.

Next compare actual WEC-Sim excitation against the unchanged offline reference
on identical coefficients and wave records, quantify all six channels and kernel
window/quadrature sensitivity, then compare all six motion histories and the
total-surface export. The first force-only comparison may use the explicitly held
fork; that does not qualify the new combined fork or authorize changing its pin.

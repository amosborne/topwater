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

## Fresh coefficient run

The reviewed historical-bank identity check passed physical metadata, provenance
and oriented-panel checks. Its exact common-grid gate rejected only a floating
representation difference of 4.440892098500626e-16 rad/s. Rather than introduce
historical repacking or relabeling, root ran the current producer once with the
unchanged canonical inputs and production defaults.

`native-coefficient-producer-attempt1/` completed all 1,800 genuine BEM solves
(8.5 s after Green-function initialization). Its private `coefficients/` outputs
include the full native FK/SC contract and source/input/geometry/output manifests.
The independent reviewer verified all nine manifest file hashes, producer input
and output hashes, and all 69 recorded Capytaine implementation/native file hashes.
Operator manifest SHA256:
`c57a86c780730f2fdc235cdc7503843583a5fce3aecb322bffd07ea26425826f`.

Five `capytaine.libs` DLLs were omitted from the coordinator's frozen-input list;
the producer's recorded pre-run hashes match their independently checked current
hashes. This limitation is preserved, and subsequent manifests include those DLLs.
No rerun or coefficient correction is justified by this omission alone.

The producer, source and recorded dependency checks passed; exact process absence
was independently verified and the lease released. The existing 1,374 mesh and
612 irregular-frequency warnings remain unresolved physical qualifications.
The existing-directory export warning concerns a fresh owned scratch directory,
as covered by the CLI warning triage. No motion, surface or force history was
computed in this job, and no canonical file was changed.

## First actual force comparison

`native-force-comparison-attempt2/` computed WEC-Sim's excitation using the
accepted current coefficients and the held995 full-history correction. The
origin elevation check differed by only 2.33e-15 m. The +/-20 s finite excitation
kernel used dt=0.02 s and frequency step=0.005 rad/s; radiation memory stayed10 s.

| Channel | RMS difference / reference RMS |
|---|---:|
| Surge | 1.052% |
| Sway | 2.703% |
| Heave | 0.252% |
| Roll | 2.446% |
| Pitch | 1.301% |
| Yaw | 4.353% |

The initial heave force was -37.257 N versus -5.962 N in the offline reference;
the native force exceeded the unchanged24.525 N initial force-norm limit. The
startup gate correctly rejected the run, with no body motion produced. Cleanup
was verified; exact process-tree absence and frozen identities were checked
before releasing the lease. This is a measured discrepancy, not accepted
equivalence. Its cause remains unestablished; a stored-data convolution/kernel
discriminator is next. No ramp, retiming, amplitude change or limit relaxation
was introduced.

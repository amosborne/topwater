# Corrected spectral path: private end-to-end qualification

The corrected WEC-Sim spectral path now has a complete private simple-boat result,
motion plot and video. The original boat and sea are unchanged. WEC-Sim computes
wave excitation and integrates motion; Topwater reconstructs the total free surface
afterward. Neither a reference force history nor a reference trajectory drives the
new body model. Dependency promotion and merges remain held.

## Requirements and evidence

| Requirement | Current evidence |
|---|---|
| Preserve phases, directional amplitudes and raw radiation data | Exact e24c1a7: 14 native MATLAB, 38 Python and 9 official BEMIO/H5 controls passed |
| Full canonical force agreement and calm startup | Independent acceptance f460dc37; all 3001 samples in all six channels pass unchanged peak/RMS gates |
| Actual native bus and once-only outer application | Bus acceptance 547becf2 and motion acceptance d934f305 |
| Body rest and zero radiation memory | Fresh-process/native-initializer certificate, observed zero initial force buffers and full retained CIC replay |
| Six-DOF motion agreement | Maximum position difference 0.059758 mm; relative orientation 0.001569854 degrees, below 1 mm/0.01 degree limits |
| Fresh total incident/diffraction/radiation surface from the new motion | Assembly a6c096b4 uses unchanged production source 240198ee and operator f1f788f9; completed cache 09d39c73 binds the full saved new motion |
| Standard result and plot | Final Result validation, finite geometry, exact accepted poses, 3001 motion samples and 751 surface frames, including both 0/60 s endpoints |
| Requested visualization | Renderer 141046c5 has 50 passing focused controls; full public CLI renders the new Result |
| Common plot scales and seven symmetric nice ticks | Plot source f39a7572 has 23 passing focused controls; normal production plot exporter is used |
| Complete video | 751 decoded frames, 800x800, 12.5 fps, 60.08 s encoded duration |
| Exactly eight workflow files and consistent logs | Final audit binds all eight hashes; every command log ends successfully |
| Warning triage | WEC final report 147c92ad classifies all 34 result warning rows; sea/video have zero application warnings |
| Preserve public reference and held dependencies | Full frozen pre/post identity checks and protected-public hashes pass |
| Terminal process/serial execution proof | Every allocation exited zero; exact owned processes absent through host CIM before idle readback and final lease release |

Independent final-delivery acceptance: artifacts/workflow/private-delivery-review-v1/
acceptance.json, SHA76f5616aac1a515ba26d8ee49b7f9277ba938b07e175bb96d96b44df94326215.
The non-author adversarial reviewer accepted the complete technical artifact scope,
including all eight actual file hashes, final report624ed7bb and final warning triage.
The earlier 286-test workflow run predates the final renderer/plot changes; it is
not evidence against those later source hashes. The focused checks above cover
the final versions. No images were viewed during full rendering or final audit.

The final style retains water, full gray hull mesh, black hull opening outline and
solid blue moving waterline in that layer order. A fixed dashed mean-waterline
outline has no label. The README explains the mean-hull linearization and explicit
exterior-trace approximation. Separation of the two waterlines is qualitative,
not a numerical accuracy certificate.

## Artifacts and physical limits

Private folder: artifacts/workflow/native-spectral-private-candidate1/. Exactly:
boat.npz, sea.npz, result.npz, result.png, video.mp4, sea.log, result.log, video.log.
Reports and recovery cache remain outside the folder. Result NPZ hash starts
6d298d9e, plot 275ac4c5 and video e5fce0e4. Complete hashes are in
native-spectral-private-evidence1/final-report.json.

New total-minus-incident projected-area/time RMS is 5.98768 mm over the full patch
(2.24251% of incident RMS), and 7.99367 mm in the defined one-beam band (2.96028%).
Maximum resolved difference is 65.0892 mm. These quantify disturbance contributions,
not numerical error. Direct old/new pointwise subtraction is marked incomparable
because saved XY grids differ; no fluid regridding or fitting was performed.

Known warnings remain substantial: calm restoring departures of 183% in moment
and 132% in force, 338 flat-reference unavailable samples, and 29.15% represented
wet-set departure. Coarse-panel/irregular-frequency and provisional causal-field
projection qualifications remain open. Issues 31-34 are not resolved by this
case-compatibility result. A missing force-field diagnostic is not evidence of
force-field consistency. The total local startup elevation bound is 0.645412 mm
on the explicitly stated conservative exterior support; no total disturbed-fluid
velocity bound has been established.

Higher radiation-warning counts reflect rebuilding an empty field source bank
rather than reusing old coefficient-stage sources. The operator did not change,
and the count explanation does not waive its physical limitations. The raw video
transcript's duplicate frame-rate option notice has no observed cadence failure;
actual full decoding and frame-rate checks passed.

## Runtime and holds

Measured coupled-motion validation: 23.0 min; surface/export: 16.1 min; ordinary
render CLI: 153.8 s. Within surface geometry, wetted-area diagnostics take 335.4 s
and delta metrics 143.6 s; waterline intersections take 11.5 s. Native compilation,
integration and audits lack separate retained timers. See the complete
[timing breakdown](../../workflow-timing-2026-10-27.md).

Motion used exact WEC topic e24c1a7 plus ONLY the independently qualified PR9
inertia library in a detached private worktree. This does not merge PR9 into PR14.
Own-fork issue12/draft PR14 carry the accepted force/bus/motion status. The held
production pin, topwater/reference branch, public physics and upstream submission
holds are unchanged. Topwater-native remains a nonworking prototype.

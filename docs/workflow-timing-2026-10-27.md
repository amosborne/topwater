# Simple-boat validation timing

Measured on the current Windows host during qualification of the corrected
WEC-Sim spectral path. These are validation and delivery allocations; they include
checks that are not a pure body-integration benchmark. No performance change was
made during this validation.

| Final delivery allocation | Controlled wall time |
|---|---:|
| Coupled-motion validation | 1382.281 s (23.0 min) |
| Total-surface reconstruction, diagnostics and NPZ/PNG export | 968.172 s (16.1 min) |
| Ordinary full render CLI, including process overhead | 155.531 s (2.6 min) |
| Final eight-file audit and complete video decode | 29.219 s |

The separate force-only replay took 497.594 s and actual signal-bus qualification
386.656 s. Those are additional validation runs. Review time, preparation outside
these allocations, earlier component tests and failed studies are not included.

## Surface/export breakdown

| Stage | Wall time |
|---|---:|
| Radiation transfers | 96.3 s |
| Radiation convolution | 17.1 s |
| Combined diffraction transfers | 77.8 s |
| Total-field synthesis | 137.5 s |
| Surface geometry and diagnostics | 505.6 s |

Within geometry/diagnostics, wetted-area departure checks consumed 335.4 s,
waterline-based field-delta metrics 143.6 s, opening-clearance checks 14.7 s and actual
waterline intersections 11.5 s. These diagnostics dominate this stage; the waterline
intersection itself is a small fraction. About 134 s of the complete export
allocation remains outside the above stage timers, including setup, identity
checks, archive/plot publication and comparison. Those operations do not yet
have separate complete timers.

## Motion and rendering detail

The native transcript records 38.763 s of WEC preprocessing and 0.279 s of mass/bus/
template setup. There is no retained timer separating final compilation, body
integration and postchecks. File timestamps bracket roughly 6 min of force replay
and 14 min of simulation plus audits; they are approximate checkpoints, not stage
stopwatches. Do not label the full 23 min allocation as pure solver time.

The render's own timer reports 139.766 s; the normal command log reports 153.797 s,
and the enclosing controlled process 155.531 s. Internal rendering includes 94.507 s
of geometry/artist preparation, 39.026 s of drawing, 0.904 s of encoder sending,
2.000 s of encoder closing and 2.641 s of complete video verification. Sending can
include backpressure and is not a pure encoder CPU measurement.

The next performance investigation should first separate native model compilation,
integration and audits. The measured surface costs identify wetted-area and delta
diagnostics as candidates for later profiling. No speed change is implied by this
accounting, and physical warnings must remain visible.

Evidence: command logs in the private candidate folder and root allocation
receipts under artifacts/workflow/native-spectral-private-{assembly,render,final}1.
Machine-readable accounting is in native-spectral-private-evidence1/timing-report.json.

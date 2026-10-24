# Frozen-source Duffy–Gauss correction: prospective diagnostic

Scope: correct the under-resolved composite pressure integration, retaining the
same864 Capytaine indirect source densities, mesh, physical model and Green settings.
The earlier centroid rule and native body-refinement gates remain failed; no new
reference, source or runtime qualification follows this conditional experiment.

Use Gauss–Legendre nodes/weights mapped to[0,1], barycentric coordinates
((1-u)(1-v),u,(1-u)v), relative triangle weights2(1-u)wu*wv. Declare native centroid
control plus orders8/12/16:864/55296/124416/221184 targets. Test constants, linear
and polynomial moments, positive weights, per-face area and strict ownership.
Verify paired physical geometry. Reflect TARGET positions/owners only, evaluate
the entire frozen S*sigma field independently at both halves, and use actual
coarse normals/areas. No mirrored fields, matrix averaging or forced zero entries.

First reproduce the native potential and all36 loads at1e-10 scaled tolerance.
Keep the existing0.001 relative +1e-6 normalized stabilization gates on the two
finest Gauss increments/trend, plus every failure mask/ratio. Do not assume a
spectral rate or Richardson extrapolation for edge-log traces. Retain all targets,
weights, six complex potentials and all36 raw complex/signed/absolute terms in
A_f-A_c=(A_Q-A_c)+(A_f-A_Q). The remainder includes changed source/collocation
solution and finest-output quadrature error; it is not uniquely source error.

Archive actual producers, controls, pinned complete paired-study bundle and
physical input before native work. Serial job, one native thread, <=100target
blocks by864sources,600wall seconds including preparation. Check time before/after
native blocks and final completion. This is offline reference work, not Topwater
simulation speed. Missing required orders, invalid targets, centroid mismatch,
nonfinite values or failed imaginary/load checks cannot pass.

Deliverables: complete diagnostic or concrete rejection, focused moment/ownership/
mirror tests, independent review, durable producer/input/result evidence, tracker
and handoff. Stop after interpreting quadrature stabilization and attribution.
No BEM re-solves, new formulation, remeshing, optimization campaign or motion run.
Full reference and source-discretization flags remain false even if this subgate
passes. The broader canonical time-domain validation goal remains active.

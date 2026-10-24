# Canonical Capytaine / BEMIO / WEC-Sim reference

The minimum 17-panel reference now has an independently reviewed ten-second,
six-DOF WEC-Sim trajectory from actual rest, with waves arriving from upstream.
It contains 501 native samples at dt=0.02 s. Final centered isometric video and
six-motion plots passed independent literal encoded-frame audit. Temporal and
mesh motion accuracy are unestablished; this is the first minimum-case result,
with no convergence claim.

Read the [workflow](../../docs/capytaine-wecsim-workflow.md),
[current authorization](../../docs/checkpoint-handoff.md) and
[physical contract](../../docs/experts/shared-contract.md) before choosing work.
Topwater implementation and refinement are outside this producer milestone.

[canonical-boat.json](canonical-boat.json) preserves the shared 24 ft hull,
minimum wetted triangles, equilibrium pose, full dry inertia and three-component
sea. Consumer startup uses the frozen, independently checked dispersive packet
in local artifacts/reference/minimum-20261024/packet-attempt3. Its noncompact
tails passed the prospectively frozen local-calm/force/arrival gates, including
the full physical-hull planform. Signals and body startup were not clipped;
there is no temporal ramp, warmup, PTO, mooring or empirical damping.

The reproducible producer pipeline has two stages, each requiring a bounded
assignment and the serial numerical lease:

1. [export_case.py](export_case.py) performs native minimum-mesh hydrodynamics
   and retains exact geometry/M/CG/Ainf identities. Original completed results
   are in local artifacts/reference/minimum-20261024/native-attempt2.
2. [correct_hydrostatics.py](correct_hydrostatics.py) reuses every frequency
   coefficient unchanged and computes exact degree-two hydrostatic moments with
   supported native GL2 on those same 17 triangles. It adds the explicit
   body_name compatibility alias needed by official WEC-Sim v7.1 BEMIO.
   [check_hydrostatics.py](check_hydrostatics.py) independently verifies analytic
   triangle moments and actual full hull CB/K. Accepted output is native-gl2;
   the original centroid hydrostatics remain preserved as diagnosed evidence.

Official BEMIO import/H5 writing preserves exact native CG/CB/volume/full K and
Ainf with reviewed units, DOFs and one phase conversion. WEC-Sim runtime reads
that H5 and uses the full dry-plus-added inertia with independently reviewed
consumer boundary corrections. Actual accepted motion is in local
artifacts/reference/minimum-20261024/case-attempt1. Its native Transport Delay
uses zero initial residual-acceleration feedback then forward extrapolation;
this identified dt=0.02 discretization is documented, not a frame/inertia error.
Small Newton-Euler residuals establish internal consistency, not motion accuracy.

Accepted video/plots are in case-attempt1/render-attempt3. Video preserves native
0..10 s states at 25 fps (251 frames, 10.04 s encoding including one held final
frame), with equal physical units and no exaggerated motion or wave height.
The displayed surface is explicitly prescribed incident elevation. Original
draft and rejected render iterations remain preserved.

Evidence and limitations:

- [Producer/export and exact hydrostatic review](../../docs/experts/capytaine_expert/minimum-producer-2026-10-24.md)
- [Joint interface, inertia and frame review](../../docs/experts/capytaine_expert/consumer-interface-review-2026-10-24.md)
- [Independent native motion review](../../docs/experts/capytaine_expert/native-motion-review-2026-10-24.md)
- [Current maintained WEC-Sim fork/pin](../../docs/experts/wecsim-fork.json)
- [Capytaine current state](../../docs/experts/capytaine_expert/current-state.md)

[history/](history/) preserves the abandoned Python time-domain prototype under
its original policies. It does not supply the accepted WEC-Sim dynamics. Future
mesh variation is judged using all six whole motion histories; intermediate
coefficients remain diagnostics. No automatic study follows this checkpoint.

# Fixed-world boat boundary reference

This experimental geometry check clips the actual 24 ft, 2500 kg reference boat
against a fixed piecewise affine free surface. Its 15 parent nodes sample the
observed SoCal baseline realization with seed 42. The boat retains its independently
computed mean-water equilibrium pose. This snapshot does not initialize a nonlinear
incident/scattered field or qualify a fluid solver.

Run from the repository root:

```powershell
uv run --locked python -m tests.experiments.free_surface.boat_world_boundary --output artifacts/boat-world-boundary-seed42.json
uv run --locked pytest tests/test_reference_boat.py tests/test_boat_parent_clip_proof.py
```

Use a fresh output stem: the driver refuses to overwrite existing evidence. It
saves the prescribed elevations, both actual-height and mean-water incident
potential traces, complete explicit wave realization, body pose, boundary mesh,
face ownership and exact-parent interpolation weights. Failed runs retain their
partial polygons and error metadata.

The frozen successful report is
`artifacts/boat-world-boundary-seed42-checked.{json,npz}`. It contains 185 vertices
and 366 triangles: 110 hull, 130 free-surface, 78 wall and 48 bed triangles. Closed
volume, first moment and affine rigid-motion normal work agree with independent
column/solid calculations to approximately 3e-14 m³, 1.3e-13 m⁴ and 4e-15,
respectively. Parent height and potential interpolation inherit exactly. The
reported acceptance checks are geometry roundoff checks, not force accuracy gates.

The first attempt is preserved as
`artifacts/boat-world-boundary-seed42-failed-symmetry-roundoff.{json,npz}`. A
floating transverse mass-moment residual created two distinct exact intersections
separated by less than one float64 coordinate ulp. The fixture now enforces its
authored transverse symmetry only after checking integrated residuals against
summation roundoff. No mesh vertices were welded and no closure guard was relaxed.

Fraction arithmetic supplies an ownership reference, not a production mesher.
The next physical test requires the current elliptic field and nonlinear pressure
on this boundary; this geometry result alone provides neither.

For directional pressure audits, `construct(output, *, body=None, pose=None,
parent_xy=None, parents=None, parent_eta=None, parent_psi=None)` returns the mesh,
saved arrays and report. Explicit coefficients remain attached to the fixed
world parents. Every construction retains its own fresh checkpoint stem.

`body_material_samples(arrays, barycentric)` returns original material XYZ and
facet-group IDs with shapes `(body faces, quadrature points, 3)` and `(body faces,
quadrature points)`. `interpolate_body_trace(offset_arrays, offset_phi, points,
groups)` evaluates the offset P1 trace at these same material points. It searches
only matching-group wet triangles and rejects absent coverage; it does not align
gauges or assume offset node correspondence. Interior quadrature samples should
be used for material derivatives: current waterline vertices generally become
dry in one directional offset.

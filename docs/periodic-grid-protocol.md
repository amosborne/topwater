# Periodic grid and field candidate: prospective checkpoint

Scope: construct a body-fitted mean-plane reference sampling grid, evaluate the
archived periodic sources on all four meshes, and test the proposed radial
pullback geometry before adopting any nonlinear restart. This advances the same
canonical one-second case; it does not replace full motion acceptance.

Candidate geometry: polar annulus about the frozen world CG in x/y, outer radius
24 m, with 32/64/128 uniform angular divisions plus every mean-waterline corner.
Radial rows are the exact mean waterline, then offsets 0.15/0.4/1/2.5/6/12 m and
the outer circle. A proposed finite TD cylinder has depth 24 m; radius/depth
32 and 48 m remain required domain controls before runtime qualification. No
truncation boundary condition is qualified by this sampling-grid choice.

Never evaluate the reference on the waterline itself, inside the mean hull or
above z=0. For the boundary limit use explicitly retained exterior radial
stencils at 20/10/5/2.5 mm. Linear extrapolations from each adjacent pair to zero
are diagnostics; retain both finest differences and do not silently treat them
as exact boundary values. Free-surface interior sampling points are at z=0 and
strictly exterior. Body traces already exist in the source-recovery artifact.

Compose total potential from incident plus the archived combined scattered
source, with all original complex phases. eta=i*omega*phi/g. Store per-component
and summed initial fields, every raw mesh-level value, the non-canceling sum of
max(abs(level4-level3),abs(level3-level2)), and the two finest boundary-limit
increments. Do not use Richardson extrapolation or fit wave/body amplitudes.

Prospective diagnostics (not motion qualification): each field's mesh and
boundary-limit envelopes must be finite and non-increasing within 1e-9 SI, and
their component sum must be <=0.1 mm in eta and <=0.001 m^2/s in phi. Those limits
screen field stability only; they do not convert field error to motion error.
Interpolation checks use the nested angular grids and actual retained fields.
Unresolved mapping/projection/linear-limit or motion uncertainty cannot pass.

The proposed pullback maps each mean waterline ray onto the actual posed hull
at the reconstructed boundary elevation, then blends radial displacement to zero
at the outer circle. It must insert hull-facet switches, including chine crossings,
and verify a one-to-one surface, positive projected triangle areas, exact contact
with the physical hull, no triangle intersection with the hull interior, and
an exterior inverse-map target at every reference evaluation point. A collection
of valid contact vertices alone does not establish valid contact edges/faces.
Reject invalid candidates before adopting a restart. A pinned material contact
line or an interior/above-plane potential continuation is not an acceptable fix.

Exclusions: a qualified TD truncation condition, new fluid formulation, tuned
loads, injected future Capytaine trajectories, an arbitrary-startup claim, major
optimization, and a complete restart claim from field diagnostics alone.

Archive producers, controls, physical hull/pose and source-bundle identity before
the serial job. Single native thread, <=1800 s including Green preparation and
evaluation; evaluate points in bounded chunks rather than constructing N-point
by N-panel matrices for the whole grid. All source/geometry identities must match.
Missing arrays/dependencies or invalid points fail. No numerical producer edits
until the exact process handle is terminal; use fresh output stems.

Deliverables: grid/field/map evidence (or a concrete rejected candidate), focused
geometry/phase/admissibility tests, independent review, complete durable bundle,
issue #30 update and fresh-context handoff. Stop after the field/map result and
review; keep full qualification false until compatible initialization and motion
uncertainty are demonstrated. An actual coupled preparation/RHS cost screen must
precede the one-second run and reject a route projected above 200 wall s/s.

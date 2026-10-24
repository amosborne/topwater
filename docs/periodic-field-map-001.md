# Periodic field and surface-map candidate 001

The candidate completed all 12 source-mesh/frequency evaluations at `05cc94e`.
Exact session 86292 is terminal exit 1: complete=true, passed=false is intentional.
Offline evaluation and mapping took 124.014 s, not Topwater simulation time.
The [prospective protocol](periodic-grid-protocol.md) and
[controls](../validation/periodic-grid-controls.json) were committed before launch.

All 1,518 Capytaine evaluation points are strictly outside the frozen mean hull
at z=0. The 138-angle reference grid includes every mean-waterline corner.
Waterline values use explicitly retained 20/10/5/2.5 mm exterior stencils, never
an interior or above-plane continuation. Sources were reused from recovery 003;
no new radiation/diffraction solves were needed. Green evaluation used blocks
of at most 100 by 7,700 entries and one native thread.

All three surface maps (full/half/quarter amplitude) passed geometric checks:
145 contact vertices, seven inserted facet-switch angles and 2,030 triangles.
Every contact edge shares a physical hull plane, projected triangle areas are
positive, and clipping whole triangles finds no hull-interior overlap. Minimum
conservative opening clearance is 0.584 m. This establishes the represented
surface geometry only, not a complete fluid domain or compatible fluid state.

The original field diagnostics failed. Mesh and exterior-limit trends are
non-monotonic under the declared 1e-9 SI floor. Their combined elevation envelope
is at most 11.076 micrometres and potential envelope at most 0.000051836 m^2/s,
but small values do not waive the failed trend requirement. Total-field angular
representation is also inadequate: its retained envelope reaches 9.667 mm in
elevation and 0.034599 m^2/s in potential. No tolerance was relaxed.

An archived array-only diagnostic separates incident and scattered contributions.
At the outer row the incident elevation envelope is 9.667406 mm versus the
total's 9.667362 mm. Across all rows the scattered envelope is at most 0.065735 mm
and 0.000301874 m^2/s. This points to preserving the already-known incident waves
analytically when representing initial fields, rather than blindly refining the
whole far grid. It does not qualify a revised representation or remove geometric,
source-trend, projection, domain or field-to-motion uncertainty. Angular
interpolation is grid/initial-state representation error, separate from Capytaine
source-mesh uncertainty.

Initial elevation ranges from -0.002671 to 0.033148 m; potential ranges from
-0.319319 to 0.304157 m^2/s. Stored arrays retain complex phases at every level,
incident/scattered samples, boundary-limit estimates, empirical envelopes and
mapped scalar potential. Future reference loads/trajectories remain external.

[Complete bundle](../validation/evidence/periodic-fields-001.zip),
[summary](../validation/evidence/periodic-fields-001-summary.json).
SHA256: `cc35d61aa1febeb6b9389d9f756304fa15a411fdbd14c01b6ae932eeabc45d65`.
Includes exact producers/controls/physical inputs, complete nested source/raw/freeze
evidence and the diagnostic's actual code/input identities. All 246 local tests
passed (session 88691 terminal exit 0). No numerical job remains live.

Independent adversarial review verified all identities, finite samples, exterior
points, incident phase/elevation signs, every envelope and failed trend, mapped
scalar potential and all complete contact-edge/triangle checks. No blocking
finding remains for preserving this explicitly unqualified candidate.

Reproduce with the locked correlation environment, restored nested inputs and a
fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.periodic_fields --raw artifacts/validation/reference-study-003.json --frozen artifacts/validation/frozen-fixed-001.json --sources artifacts/validation/periodic-sources-003.json --sources-bundle validation/evidence/periodic-sources-003.zip --output artifacts/validation/periodic-fields-next.json
```

Full reference qualification, complete nonlinear restart and actual TD simulation
remain unavailable. Body-boundary compatibility, a closed fluid domain,
field-to-motion uncertainty and nonlinear-limit checks are explicitly false.
Original coefficient-study failures remain failed. Radius/depth24 m are candidate
controls, not a qualified truncation or runtime mesh selection.
At 0.1 Hz the deep-water wavenumber is 0.040243 m^-1 (wavelength156.13 m).
The incident vertical factor exp(-kH) is 0.3807 at H24 m and 0.1449 at H48 m.
A zero-normal-velocity seabed at24 m would not preserve the infinite-depth case.
Declare outer/bottom conditions and initial compatibility before a coupled RHS;
neither candidate depth has earned a domain-error claim.

Next: correct known-incident representation and construct the minimum mixed-
boundary pressure/body RHS needed to rank initialization errors and measure real
coupled cost. One t=0 acceleration difference cannot bound one second of motion;
an unvalidated RHS's sensitivity cannot be the sole qualified reference-error
estimate. Total-surface geometry and discretization still require controls even
when the incident field is analytic. Stop before a long run if coupled preparation
and step costs project above 200 wall s/s; major optimization remains held.

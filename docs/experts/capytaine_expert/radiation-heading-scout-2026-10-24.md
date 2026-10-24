# Combined radiation and heading diagnostic proposal

Status: independently source-reviewed, allocated and executed once, terminal0.
Original proposal below is retained; executed evidence and limits follow.
Coordinator settled this implementation scope after the stored consumer comparison
revealed a material native/CIC discrepancy. Earlier exports and baseline remain
immutable. This note authorizes neither a correction nor a trusted canonical case.

## Trigger and question

Stored evidence:
`artifacts/reference/generated-bank-consumer-20261024/compare-attempt1/summary.json`.
April particular histories differ native/CIC40 by0.39567deg roll and0.05393m heave.
Native17-panel Ainf_heave=10541.39kg; A(12)_heave=2311.19kg, difference-8230.20kg.
B(12)_heave=3864.21kg/s remains substantial. The consumer separately bounded the
40s numerical quadrature/window effects; their smallness does not attribute the
larger finite-band/native discrepancy to either missing tail or mesh error.
These particular histories contain no held-release transient/free drift trajectory.

The17-panel exporter uses centroid hydrodynamics, native default GF/solver and no
lid. Logs explicitly warn about mesh resolution and possible irregular frequencies,
including omega5.25. Warnings are heuristics, not proof of an observed spurious mode.
Largest panel radius3.678m exceeds wavelength/8=.280m at5.25rad/s. At12rad/s,
wavelength/8=.0535m; even uniform68 panels remain far underresolved at high frequency.
At64rad/s wavelength is about.015m; the17-panel64 point is an asymptotic stress
probe, emphatically not resolved radiation-tail data.

## Exact bounded scope

| Output | New solves | Stored reuse |
| --- | ---: | --- |
| radiation-17.mat |16,24,32,48,64rad/s x6=30 |1,2,12rad/s A/B and direct Ainf |
| radiation-68.mat |1,2,12,24,48,inf x6=36 |none |
| jan05-offset90-pairs.mat |384 exact diffraction pairs |unchanged global realization/heading frame |
| Total |450 |no auto cascade |

All36 raw radiation couplings survive with original asymmetry. Infinite-frequency
added mass is separate from finite A/B; no meaningless infinite B is fabricated.
Uniform four-way subdivision preserves each physical planar parent triangle;
it does not change the supported boat. No lid, alternate quadrature, stiffness,
mass, damping, phase, amplitude, direction fitting or new physics is introduced.
No272-panel branch is implemented. That would require a separate decision/grant.

The second output uses ONLY existing
`ndbc-46086-20250105T1140Z/heading-offset-90deg.json`. This is an internal rejection
hypothesis about strong global roll with smaller relative-wave geometry, not a
canonical heading freeze. All384 components, including zero amplitudes, retain
global realization, local phases and amplitudes. Pair schema matches the reviewed
first producer, native exp(-iomega*t), same meanCG and frame conventions, no extra
CG phase or BEMIO conjugation. No full heading matrix or additional RIRF is exported.

## Source, command and resources

Source `reference/capytaine/export_radiation_heading_scout.py`; controls
`reference/capytaine/radiation-heading-scout-controls.json`. It imports only
the existing `triangles_at` subdivision helper from `export_case.py`, whose
identity is frozen. Archived draft manifest:
`artifacts/reference/radiation-heading-scout-20261024/proposal-draft.json`.
Syntax inspection passed; no source execution or BEM launch performed for readiness.

```powershell
.venv/Scripts/python.exe -u reference/capytaine/export_radiation_heading_scout.py --controls reference/capytaine/radiation-heading-scout-controls.json --output artifacts/reference/radiation-heading-scout-20261024/attempt1
```

One CPU thread, maximum300s checked between solves. Conservative estimated20-90s
is unmeasured for this68-panel/high-frequency mix; not a performance guarantee.
Estimated memory<1GiB for these small matrices, also unmeasured. One native Green
initialization/solver is shared across both meshes and the paired bank. Preserve
partial output/checkpoints/logs on error. No automatic retry or refinement.

Before launch: independent source review and coordinator READY allocation, archive
actual source/controls/helper/inputs/dependencies/Python identities, then CreateNew
lease and successful claim, SEPARATE successful matching readback, then run. Reuse
the unchanged reviewed guard. Inspect every actual command exit. Exact process and
child terminal, idle write under lease, SEPARATE idle readback, final release.
Stored consumer signed-gap work may finish before this claim; do not overlap jobs.

## Expected checks and interpretation

Geometry tests are per physical parent, not merely total cancellation: oriented
area, area and exact first/second moments sum over children; child vertices remain
coplanar and barycentrically inside parent; normals preserve winding.68 triangles
must equal4x17, with identical CG/M/K. Stored17 A/B at1/2/12 and Ainf are copied
by pinned file hashes and exact omega indices; output readback must be exact.
All finite new A/B/Ainf arrays and all36 entries are preserved without repair.

Pairs: exact bank arrays/indices/CG/rho/g/yaw/R/t/DOF/source hashes, finite384x6
FK/diff/H, H=FK+diff, analytic imaginary quarter-cycle sign, independent actual
centroid incident-pressure FK normal/moment integration, all-component nonzero-phase
scalar-versus-complex summation, and global/local eta closure at three points/times.
No consumer applied-force or dynamics validation is claimed by these checks.

Compare low-band1/2 changes across17/68 to identify immediate spatial sensitivity;
compare direct Ainf across meshes and A(12/24/48) at both meshes. At fixed mesh,
high-frequency A approaching its own Ainf suggests a cutoff contribution but
cannot establish physical asymptotic accuracy on this unresolved mesh. Changes in
Ainf show discretization sensitivity; neither pattern alone uniquely isolates
irregular-frequency effects. Few points cannot integrate the omitted KK tail or
qualify the40s radiation model. Preserve ambiguity and all raw outputs.

WEC owner later reports native/CIC particular/global/relative proxies for January90
without treating proxies as actual wet set or held-release trajectory. The geometry
study's5/10deg relative flat-water departures remain material. No fitted Ainf,
offdiagonal zeroing or added damping is an acceptable inference from this scout.

Stopping point: independent reviewed source/manifest/resource proposal ready for
coordinator allocation; after any later allocated run, terminal/raw outputs/checks,
not a production correction or full sweep.

## Executed result

Source remains exactly ff75f6bc210ce6aaf4b47825b2e2ba6d06e7d23826b1cd85c509328ebd86c7c9.
Researcher source PASS preceded allocation. Only status changed in READY controls,
SHAdef7dd6487a5bbe564136613415969d51f1e70fa29a662e630525b56316977f8;
frozen proposal SHA51c5a3b66f50aeae4f7e487571144f5566cbbd659a32a34431d85e2452839f7c.
Output `artifacts/reference/radiation-heading-scout-20261024/attempt1`.
Manifest94c3c29dd0de40cdc6302726a585519e3f90c9964caa0cf61ff4287f5127e9de;
postcheck1907a437f7a7d9f256935a0ecd2cc2861c063e1a33c918f4f07c86cfb1048296.

Separate Claimfbf8a3 exit0 and Verify1f59c9 exit0 preceded Run session29495;
producer39156/guard5000, UTC launch18:26:08.189316.450 calls completed56.9944896s,
one thread, producer/guard/tool terminal0 (chunk8052b0); exact guard wait and both
process identities absent after tool terminal (82101d exit0). Source spawns no child.
Idle written under lease393eed exit0, SEPARATE idle readback7e0a65 exit0, final lease
removal e7cecc exit0. No ledger writes after release, no further numerical allocation.
Failed early Get-Process sampling found producer already absent before tool output
drained; it was not used as standalone completion proof. Actual wait/tool proof retained.

All geometry/per-parent moment/plane/winding/barycentric checks pass. Stored17
A/B(1,2,12) and Ainf exact; full raw finite shapes/M/K/CG and frozen/output hashes
pass. Pair schema/indices/frame exact; independent FK max1.4553e-11, all-component
scalar/complex force sum5.8208e-11, spatial eta closure1.6654e-16m.

| Quantity |17 panels |68 panels |
| --- | ---: | ---: |
| Ainf heave,kg |10541.392 |9169.628 |
| Ainf pitch,kgm2 |6031.705 |16660.115 |
| A heave at1rad/s,kg |19767.677 |19242.534 |
| B heave at1rad/s,kg/s |5311.986 |5518.611 |
| A heave at2rad/s,kg |12966.470 |12927.137 |
| B heave at2rad/s,kg/s |15108.767 |16523.086 |
| A heave at12rad/s,kg |2311.191 |-506.527 |
| B heave at12rad/s,kg/s |3864.212 |-48117.677 |
| A heave at24rad/s,kg |2589.950 |8894.439 |
| A heave at48rad/s,kg |2546.980 |7498.892 |

17-panel Aheave at64=2619.575kg remains far from its direct Ainf even though
Bheave=3.011kg/s. Thus this finite/direct mismatch is not resolved by merely
raising the same coarse-mesh finite probe frequency.68-panel A(24) is closer
to its Ainf, but A(48) moves away again. Severe underresolution prevents using
either pattern as a physical tail model. Negative68-panel diagonal heave B at12
violates passive pure-heave radiation behavior under the established damping
convention; retain it as a numerical/pathology diagnostic, not a physical fluid
coefficient or repaired value. Strong mesh sensitivity also appears in direct Ainf.
No unique attribution to irregular frequencies, finite cutoff or discretization
is established. No272 run, lid change, coefficient fit or production replacement.

Independent executed review requested directly with exact artifacts and terminal
evidence. January90 table is available to the consumer for separately allocated
native/CIC particular/global/relative proxy comparison, not canonical adoption.

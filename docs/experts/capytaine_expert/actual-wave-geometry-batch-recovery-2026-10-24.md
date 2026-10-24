# One bounded geometry recovery: source only

Coordinator authorized source/proposal implementation only. No fixture, benchmark
or geometry execution; renderer preview owns the next numerical slot. Old scalar
source, attempt1 outputs/policies and all physical inputs remain unchanged.

New files:

- `reference/capytaine/actual_wave_geometry_batch.py`, SHA61c3598f3478b8d1cf98fcb915c138a87919758f7f722200c4926476e348bf26.
- `reference/capytaine/actual_wave_geometry_batch_core.py`, SHA3b04e8f28e8c9787b78c1d622fbfd7a232c15375e79e21555ed05cb6a0582d1a.
- `reference/capytaine/actual_wave_geometry_batch_config.json`, SHAa411d7aa4310bed0d89546e1c96ec9630d0caa8afd069085ace6226792f3e420.

Draft archive `artifacts/reference/actual-wave-geometry-seeded-20261024/proposal-batch-source-only.json`,
SHA96788898801033208ded227d9d387766d2a966caa3b3dbc6daa2840b31da3f97.
AST inspection passed. Independent source/math/oracle review PASS at the exact
runner/core/config hashes above; no READY
status or allocation. Python and existing NumPy/SciPy/threadpoolctl identities
pinned; no new dependency or Topwater solver work.

## Preserved mathematics and changed execution

Analytic384-component eta/gradient, spectral Hessian bound, two-sided roundoff
allowance, physical material faces, mean-wet intersection/XOR inequalities,
rim clearance and actual RzRyRx/CG/frame semantics remain unchanged. All3001 native
samples and strict grid identity remain required; no decimation/frame selection.

New core implements the same convex polygon clipping order and<= comparison rule
over bounded padded arrays (capacity6 vertices after two half-plane cuts). Area is
the norm of the planar polygon area vector; mathematically it equals the scalar
convex fan areas. Lower/upper current wet sets and mean-wet intersections retain
the same XOR lower/upper set construction. Crossing division occurs only on actual
opposite classifications; equal/parallel edges avoid0/0. Degenerate polygons yield
zero area and cannot manufacture a positive wetted area.

Triangle diameter/roundoff calculations and all clipping passes are batched;
spectral constants are cached once. Refinement chooses up to256 parent triangles
at a time, with parent-major child ordering and repeated material-face identifiers.
Hessian/interval math and depth8 stop rule are unchanged. Original scalar rim code
is retained; cached analytic field/roundoff evaluations also reduce its overhead.

Prepared oracle fixtures compare new batch bounds with untouched scalar
`area_bounds` for flat/translated material triangles, zero-area triangles,
all-wet/dry, parallel and equal-valued edges, exact contact, tiny signed values
and ordinary cuts at several epsilon values. Original independent mean-area,
flat-plane, reflection-parity and global/native field fixtures run first in any
future grant. All fixtures remain UNEXECUTED at this source-only checkpoint.

## Prospective reporting precision and caps

Area/XOR interval width target1% of mean wet area (0.1397955m2), rim-clearance
width10mm. They are prospective numerical reporting widths, not user physical
requirements or a change to any physical gate. Earlier0.1%/1mm prefix is preserved.
Footprint curvature retains explicit bounds and no hidden adequacy claim.

ONE prospective job:1CPU,600s,10million triangle evaluations, depth8,65536 active
triangles per sample and2048 wave-point batches. Numeric leaf arrays at that cap
occupy tens of MiB; fixed batch polygon arrays are small. Together with unchanged
libraries/inputs this is designed below1GiB, not a measured peak-memory result.
Time/triangle/depth/memory exhaustion preserves coarse bounds/unknowns and exact
unprocessed indices. No precision relaxation or automatic retry.

```powershell
.venv/Scripts/python.exe -u reference/capytaine/actual_wave_geometry_batch.py --config reference/capytaine/actual_wave_geometry_batch_config.json --output artifacts/reference/actual-wave-geometry-seeded-20261024/batch-attempt1
```

Before launch: independent source review, coordinator allocation, READY status
and fresh exact source/config/input/dependency freeze. Use unchanged claim guard,
separate successful matching readback and actual returned command exit checks.
Exact terminal/child proof, idle under lease, separate idle proof, release final.

## Cost rationale and acceptance

Attempt1 measured123 poses/300s,740489 triangle evaluations (~6000/pose). Full3001
coverage in600s needs about12.2x saved-pose throughput; use>=13x as a planning
margin. Coarser numerical widths should reduce refinement demand; batching removes
per-leaf Python polygon loops and repeated spectral reductions. Neither gain has
been measured. Source includes separate field, batched clipping and adaptation
timing counters; their scope excludes scalar fixture/rim overhead and is reported
explicitly, so unexplained wall time is not attributed to a fabricated hotspot.

Acceptance for a later run: all prepared fixtures/oracle checks pass; unchanged
actual-input/frame/source identities; all3001 samples processed or explicit
incomplete status; valid finite ordered area/XOR/rim bounds and exact coverage;
observed stage/overall timings, resource/limit events, terminal evidence. Processing
all samples does not itself establish physical linearity, kernel adequacy or
continuous-time behavior. No nonlinear force, trajectory model or BEM is added.
Independent reviewer additionally requests an existing-output postcheck after any
future allocated job: new area/XOR enclosures at matching0..2.44s prefix times
must overlap preserved higher-precision attempt1 intervals. Containment is not
required because refinement order differs. This compares stored outputs only,
not another geometry execution. No performance or executed-fixture PASS exists yet.

Stopping point: independently reviewed source-ready proposal pending parent grant.
No hours-long run, measurement campaign or additional followup is implicit.

## Reviewed reporting-cost correction

Coordinator identified that the initially reviewed batch runner inherited full-prefix
JSON rewrites every20 samples. That was quadratic reporting work and potentially
significant after the numerical backend sped up. A reporting-only source correction
now appends each row once to geometry-results.jsonl, flushes compact manifest/range
checkpoints every20, and writes geometry-results.json plus explicit unprocessed
indices once at terminal. Final elapsed includes full geometry serialization;
progress/checkpoint cost is not hidden. Core/mathematics/precision remain unchanged.

The source now performs the stored-prefix postcheck itself. Old geometry file
SHA6bc8bea487932835538a1c5dfbfe1203c1f8325ec53ab0f1357296a67c9b0196 is frozen in
config/proposal; each matching wet-area/XOR/available-rim closed interval must
overlap. Disjoint bounds raise an error and mark complete false. Missing intervals
or coverage remain unknown, not proof. Independent reporting/postcheck source PASS
received at the final runner/core/config identities above. No numerical test or
fixture was run for this correction. READY status/new HEAD freeze remain pending.

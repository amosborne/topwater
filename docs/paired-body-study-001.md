# Paired body-only study: uniform family rejected

The complete study at source `223050d` ran all four Topwater levels and all 24
Capytaine infinite-frequency calls; the numerical stage took 9.930 s. Exact session 17509 is terminal
exit 1 intentionally: complete=true, passed=false. This is offline reference
work, not Topwater time-domain throughput. The [prospective protocol](paired-body-protocol.md)
and [controls](../validation/paired-body-controls.json) preserve the original gates.

The physical hull, calm pose/world CG, full dry inertia, density 1025, gravity 9.81,
infinite depth and zero speed are unchanged. A declared y>=0 half is reflected
before solves; matrices are never averaged or forced to zero. Four exact midpoint
levels have 54/216/864/3456 panels. Every level preserves all seven facets, moments
and shared edges. Actual Capytaine radii match the farthest-vertex convention.
Both resource preflights passed; the finest Topwater S/K allocation is 191,102,976
bytes. All geometry/source/input identities are retained.

Reflection greatly reduced the previous seed's odd-parity numerical couplings:
maximum normalized values are 3.34e-7 in Capytaine and 1.08e-14 in Topwater across
these levels. That corrects a mesh-symmetry defect; it does not establish physical
accuracy or convergence. All 18 even-parity entries fail both reference-increment
gates and both candidate gates. Reference trend also fails for Surge/Heave and
Surge/Pitch entries. Maximum reference previous/latest tolerance ratios are
2091.4/1150.0; candidate mesh/comparison ratios are 63.69/309.45. These ratios are
not percentages of true error. Every failure mask, ratio and physical uncertainty
entry is stored, with all 36 raw entries for each solver/level.

| Capytaine diagonal |54 panels|216|864|3456|
|---|---:|---:|---:|---:|
| Surge, kg |108.82|98.06|93.90|91.72|
| Sway, kg |692.02|626.97|594.88|577.79|
| Heave, kg |10360.37|9139.44|8437.75|8126.08|
| Roll, kg m^2 |1073.83|1088.91|1042.18|1011.03|
| Pitch, kg m^2 |5631.90|16273.03|16947.63|16518.09|
| Yaw, kg m^2 |545.97|1270.67|1365.02|1335.52|

All six auxiliary potentials are now retained at every level for both solvers,
alongside Capytaine normalized indirect sources and boundary data. Symbolic
infinite-frequency traces explicitly cancel their -i*omega factor. Original
S/K operators verify these traces before the cache becomes LU; all eight full
36-entry load matrices reconstruct independently from the saved potentials.
Maximum recorded internal operator/trace ratio is 0.0004997 of its declared
tolerance. Internal consistency does not rescue the failed physical gates.

Uniform refinement preserves a poor seed shape: 48.15% of panels have aspect
indicator >1000 while covering 0.2116% of wetted area. An archived array-only
diagnostic exactly redistributes each measured matrix increment into child-cell
pressure contributions. In the latest Capytaine increment, 68.7% of the failed-
entry-weighted absolute contribution is below 0.1 m depth; the slender-panel subset
accounts for 0.0172%. Topwater gives 67.2% and 0.1403%. These are measured increment
allocations, **not** continuum error bounds, causal localization or proof that
slivers are harmless. Potential errors propagate nonlocally.

The initial study attempt stopped after its 54-panel Topwater solve, before any
Capytaine call, because of a wrong mesh-property name. It remains incomplete and
preserved. The corrected faces_radiuses API is tested during geometry preflight
on every level. Both geometry freezes have bitwise-identical numerical arrays.
No result from the incomplete attempt is treated as a reference.

[Complete evidence](../validation/evidence/paired-body-001.zip),
[summary](../validation/evidence/paired-body-001-summary.json).
SHA256: `052df7c50865340e7ef4ebaddfc48adfef191b24057646bd395a889515dc98f2`.
After extracting, run `python restore_shared_inputs.py` to restore the five
identical parent payload paths with verified SHA. Fresh restoration and complete
owned producer/input graphs were verified. All 261 tests passed; JUnit and the
prelaunch memory observation are included. Independent review reconstructed all
pressure matrices, failure masks and six increment accountings. No stage blocker
remains, no numerical jobs are live, and no reference or simulation is promoted.

Reproduce with fresh stems (restore the paired freeze before the study):

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.paired_body freeze --parent-bundle validation/evidence/mixed-pressure-001.zip --output artifacts/validation/paired-body-freeze-next.json
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.paired_body study --frozen artifacts/validation/paired-body-freeze-next.json --output artifacts/validation/paired-body-study-next.json
```

Next bounded diagnosis: hold the 864-panel Capytaine sources fixed and refine only
pressure integration at strictly interior descendant centroids. Retain all 36
signed matrices and absolute per-entry diagnostics, quadrature stabilization and
the exact telescoping A_f-A_c=(A_Q-A_c)+(A_f-A_Q). The remainder is unresolved:
it includes changed source/collocation solution and finest-output quadrature error.
Use Capytaine normalized indirect sources, not Topwater's direct potential trace.
No new BEM solves, runtime formulation change, remesh or optimization campaign.
Declare controls/limits before jobs and preserve every existing failed gate.

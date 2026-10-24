# Minimum producer interface and export attempts, 2026-10-24

Physical input unchanged: canonical SHA256
ea18844904e03bbfd72abc3396240c97b1994b3584e36c8149f8ae5f60a99e4f.
Original untracked exporter preserved under artifacts/preserved-producers with
original SHA256 a1feb6e8581695c1a7562bf1e129bbd609f96b54ab7dcac25da9b5ecdaa92942.

Consumer reviewed and approved hydrodynamic producer interface in
reference/capytaine/interface-contract.json. Incident forcing generation was
removed from active producer ownership; the old FFT-packet prototype is archived,
not accepted startup physics. Consumer separately owns and reviews actual-arrival
and excitation-IRF padding/timeline.

Producer controls: minimum17 wet triangles, omega0.1:0.1:12rad/s at propagation
0/45/90deg; six direct Ainf solves plus120*(six radiation+three diffraction)=1086.
One thread is set before NumPy/Fortran imports with OMP/OPENBLAS/MKL environment
and threadpoolctl context. Capytaine3.0.0 official export uses scipy netCDF backend.
Outputs: canonical.nc, Hydrostatics.dat, KH.dat in same root, native.mat/.npz
with dimensional full M/K/Ainf and exact worldCG/CB/volume/rho/g; boat.stl depicts
full physical hull at mean pose. Native legacy hydrostatics rounds CG/CB to4
places and KH to6 significant figures: exact auxiliary metadata is authoritative
for consumer mapping checks. No manual excitation conjugation in producer.

Attempt1: artifacts/reference/minimum-20261024/native/manifest.json and sibling
native-run.log. Terminal exit1 before solves, 0.4895053s internal elapsed. Invalid
prototype Delhommeau import was corrected after inspecting installed package
exports; no numerical result. Terminal evidence in ledger before final lease
release. Attempt1 exact input/dependency proposal remains preserved.

Attempt2 proposal: artifacts/reference/minimum-20261024/producer-run-proposal-attempt2.json,
SHA2569139bc1ffa6909fcd0d162868d227aeda8d72448915e46f9385d8d02af071e5e.
Exporter SHA2565ca23fb1d5ec2bb61fa1818a072fd0ddbe4db713c0991019626395cc1d28c981;
contract SHA25627bb72a21c46c227a63441ab274f6e203ac2f89f4882ce786cc6019a7e727744.
Exact package source/binary hashes, environment versions and Python identity are
in proposal. Coordinator allocated only this corrected attempt; sources/input/
dependencies frozen until exact process terminal. Runtime native manifest snapshots
all three executable inputs. No automatic rerun or refinement.

Command from project root:
.venv/Scripts/python.exe -u reference/capytaine/export_case.py --case reference/capytaine/canonical-boat.json --contract reference/capytaine/interface-contract.json --output artifacts/reference/minimum-20261024/native-attempt2

Soft1800s budget is checked between solves, not a terminal-process guarantee.
Coordinator execution ledger owns current liveness; this note is evidence only.
No motion/video or mesh-convergence claim is established by producer export.

Attempt2 terminal: tool session84329 returned exit0; elapsed69.9800751999719s,
all1086 solves complete. Only child launch in producer is synchronous git
rev-parse, already returned before solves; no asynchronous children are created.
Post-export check command also terminal exit0. Fourteen finite/input-snapshot/
output-hash/DOF-order/omega-heading-coverage/full-M/CG/auxiliary-matrix checks
passed in native-attempt2/producer-checks.json. Full roll/yaw dry mass entry
-70.03519269408712 retained exactly; volume2.4390243902439m3. Ledger idle was
recorded while claim held, followed by lease removal as final operation.

## Original native hydrostatic quadrature limitation (resolved separately)

Computed native CB=[-0.013037926443176352,0.0011121301769291333,
-0.10790348796654825]m. The nonzero lateral value and Kroll,pitch=19607.59 are
unphysical for the symmetric physical surface, though the triangulation partition
is asymmetric. Source-level cause: installed capytaine/meshes/meshes.py:481 uses
one face centroid with its area when quadrature_method is None.
meshes/surface_integrals.py:79 integrates n_i*x_i^2 to compute CB; this is
degree2 and centroid quadrature is not exact. bodies/hydrostatics.py:210 onward
likewise uses quadratic waterplane moments, affecting stiffness and cross terms.
Native Kroll,yaw=319.7551 and Kpitch,yaw=-27.27499 reflect CBxy errors.
These are numerical coarse-integration artifacts, not a changed canonical CG.
Linear WEC restoring about mean pose does not necessarily apply a constant CB
torque, but erroneous K injects spurious cross-coupling into motion.

No outputs were changed, symmetrized or used to redefine equilibrium. Coordinator
authorized exact degree2 triangle integration of physical hydrostatics only,
preserving the same17-panel hydrodynamic solve and all original native results.
Installed meshes.py:454 pads triangles to [a,b,c,c] in the public faces array.
Thus native Gauss-Legendre2 in meshes/geometry.py:88 has a valid collapsed
bilinear triangle map; earlier caution based on its4-vertex signature alone was
resolved by tracing that property. For r=c+u[v(a-c)+(1-v)(b-c)], Jacobian is
proportional to u. Degree2 physical integrands times Jacobian have degree at
most3 in u and2 in v, within native2x2GL exactness. No custom quadrature injection
or forced-zero entries is needed.

Separate correct_hydrostatics.py uses supported native mesh/body API for
hydrostatics only, preserves every original frequency-domain field and Ainf/M,
and writes unique native-gl2 outputs. It adds body_name as an explicit alias of
existing body coordinate: Capytaine3 exports body, while official WEC-Sim v7.1
readCAPYTAINE requires body_name. Original centroid outputs stay unchanged.
check_hydrostatics.py independently integrates triangle area/first/full-second
moments on unit and tilted translated triangles through all6 vertex permutations,
rejects the old centroid rule on x-squared, then derives actual17face V/CB/fullK
from analytic moments. Independent startup_review approved derivation and both
source files without author edits; see
../numerical_methods_researcher/hydrostatics-review-2026-10-24.md.

Correction/check proposal: artifacts/reference/minimum-20261024/hydrostatics-run-proposal.json,
SHA256e275e5ca57bace763daeb7de2aa98fb455a69d37c86b3c4aae17a4a8d3508775.
Zero new BEM solves, one thread, expected10s/max120s. Coordinator allocated the
correction/check after consumer's serial job was terminal.

## Corrected native hydrostatics and executed signoff

Both commands returned exit0 under synchronous exec chunk54917b; correction
elapsed0.8899721000343561s, total shell7.0368692s. No child processes are launched
by either script. Eight analytic checks passed, including12 triangle/permutation
regressions and the unchanged17face exact CB/fullK comparison. NativeK differs
from independent analytic moments by at most4.3655745685100555e-11.
Corrected CB=[4.0456028284341404e-16,1.8025240112851e-17,-0.13007045936669082]m.
Heave/pitch=-10123.376648462312 is retained; roll/pitch=-7.448357295059162e-12
and yaw couplings are roundoff residuals, with no projected/forced zeros.
Independent startup_review inspected executed evidence/output hashes and returned
PASS. This accepts exact hydrostatic moments and coefficient reuse, not motion
or mesh accuracy. Ledger idle recorded while lease held, release final.

Executed source and prelaunch identities:

- correct_hydrostatics.py SHA25641b85ca58ba7ad63d2276399c68c751d66d1c373ff4ff4dd9cbb298e3a898ec0.
- check_hydrostatics.py SHA25653c88c500a98e9a80b0c123197c447f749d35f0541b72742c02b5fbec734b44c.
- Exact correction/check snapshots archived in
  artifacts/reference/minimum-20261024/hydrostatics-executed-source; match proposal.
- Output manifest: artifacts/reference/minimum-20261024/native-gl2/manifest.json.
- Analytic result: native-gl2/hydrostatics-checks.json; log: ../hydrostatics-run.log.
- Corrected canonical.nc SHA256a287ce55ffcc50363a83e812e1860b5748f9e011d9e2735c74190da3a4e7fb6e.
- Corrected native.mat SHA2563262befdfbc3effd9bac3c5a9688de613481f8846325d93f1fa38b0f60983705.

Reproduce with .venv/Scripts/python.exe -u followed by each command below, under
a fresh allocated lease and unique output directory:

    reference/capytaine/correct_hydrostatics.py --case reference/capytaine/canonical-boat.json --source artifacts/reference/minimum-20261024/native-attempt2 --output artifacts/reference/minimum-20261024/native-gl2
    reference/capytaine/check_hydrostatics.py --corrected artifacts/reference/minimum-20261024/native-gl2

All original frequency-domain and other native variables are retained exactly
except intentional hydrostaticK/CB/disp_mass replacements; all reuse flags true.
Ainf and dryM untouched. body_name alias is explicit official-importer version
compatibility. Consumer must use native-gl2, not original centroidK. Next is
read-only executed consumer phase/normalization/full-inertia mapping signoff.

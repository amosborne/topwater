# Mild buoy train: required excitation coverage

Source-only proposal,2026-10-24. Latest coordinator assignment authorizes these
new producer files and this note; no numerical launch or radiation change.

The new immutable collected-buoy parent is observation
ndbc-46086-20230319T0840Z, seed42,384 carriers, Hm0=0.5351635264103863m,
bow yaw3.51348990965954rad (201.308143deg). Inputs:
`artifacts/reference/calm-buoy-20261024/export-attempt1/{case.json,components.npz}`.
Their SHA256 are b63d82d280b5fd3a70ef5dab93e707e300d1ca7cdbf156efcf66931fc15b26ee
and814644c993a957528432a6b5aed92d16d1ee94d52119e65e2e7dc26916b7eaa3.

WEC owns the physical incoming train:30s flat-top carrier windows, signed-kappa
Fourier propagation with omega=sqrt(g*abs(kappa)), full60s/rest/zero memory.
Negative kappa uses heading+pi. The384 original MEM directions and sidebands
require coverage beyond the old0/45/90deg coefficients. Direct solving each of
roughly500000 derived nodes is unnecessary for this required prototype coverage.

## Producer contract

`reference/capytaine/export_calm_excitation.py` performs exactly816 diffraction
solves:34 omega nodes0.1:0.1:3.4rad/s times24 headings0:15:345deg. FK uses native
analytic incident-pressure integration; total H=FK+diffraction. No radiation,
Ainf, hydrostatic, BEMIO, H5, CIC, mesh, lid, density, mass or geometry changes.
Same original17 physical triangles, centroid hydrodynamic quadrature, CG and six
rigid work-conjugate modes, rho1025/g9.81/depthinfinity/forward0, Capytaine3.0.0
native indirect/default low_freq Delhommeau settings as existing exporters.

Outputs `excitation.mat` and `excitation.npz`:omega×heading×influencedDOF.
The stored heading vector has25 entries:0:15:360deg. Entry360 is an exact copy
of0 and consumes no solve. Consumer MUST NOT append another seam. NPZ requested
keys are omega_rad_s,heading_rad,H_complex,FK_complex,diffraction_complex;
canonical MAT names omega_rad_s,heading_native_rad,excitation_complex,
fk_complex,diffraction_complex are also present in NPZ. Vectors are columns.
DOF order Surge/Sway/Heave/Roll/Pitch/Yaw. Units N/m and Nm/m about meanCG,
native exp(-i*omega*t), no conjugation and no extra CG phase. Full M/K/CG,
triangles, yaw, global_from_native_R/t and both parent hashes are included.

Consumer real/imaginary bilinear interpolation is explicit and UNQUALIFIED for
accuracy; heading modulo2pi, no omega extrapolation. Consumer prospectively
excludes derived omega<0.1 as numerical coverage policy and omega above
2pi*0.53345874095 by project's wavelength>=2beam rule, keeping both signed
branches and reporting removed sideband norm without renormalization. These
policies belong to the new train, not a silent change to original carrier data.
Coverage does not establish excitation interpolation/mesh accuracy or repair the
unchanged radiation operator. No band qualification is claimed.

## Prepared checks and run proposal

Prepared but UNEXECUTED: frozen inputs and original Capytaine source/binary
identity; exact original triangles/CG/M/DOF order; all816 finite full6 FK/diff/H;
independent centroid pressure/normal/moment FK check at EVERY node; exact seam;
H sum; quarter-cycle sign; MAT and NPZ array readback; input and original H5
hash invariance after export. Original H5 is a read-only frozen input.

Exact prospective command, from repository root:

```powershell
.venv/Scripts/python.exe -u reference/capytaine/export_calm_excitation.py --controls reference/capytaine/calm-excitation-controls.json --output artifacts/reference/calm-buoy-20261024/excitation-attempt1
```

1CPU, max180s checked between solves, expected60–100s including native Green
initialization (estimate, not measured), prospective<1GiB. Previous450-call
producer took56.99s and2838-call producer67.12s, dominated by initialization;
this estimate is intentionally bounded without a benchmark job. Warnings and
console output must remain in the eventual outer execution log. Existing serial
claim/separate verify/run/exact terminal/idle-under-lease/separate verify/final
release procedure remains mandatory after separate coordinator allocation.

Source SHAe7991250866fa10aa7e89aa07a123a8b87c3e8cbb3384f086f8810a2b1aa11e8;
controls SHAee4c986887ef4a7d2f62cc3359a61b1afe131069846b7a072396bf45c135ed43.
Archived draft proposal:
`artifacts/reference/calm-buoy-20261024/excitation-proposal-source-only.json`,
SHAcebd6e0e86bb3448423126ded401384a2e45d1ee9295f4ad34461273f704ed76.
Actual package-file identities and Python executable identity are frozen there;
syntax parsed with stdlib AST only. Controls remain draft pending joint consumer
and independent source review, then coordinator grant/READY-only refreeze.

Source basis: existing exercised `export_generated_banks.py` (SHA3c15c37e...),
installed Cap3 `bem/solver.py`, `bem/airy_waves.py`, `bodies/dofs.py` pinned through
the original-native-manifest package hashes. No new hydrodynamic API is assumed.

Consumer joint source/interface approval and independent source/freeze PASS
received at the original draft hashes above. All9 inputs and listed dependency
file hashes independently match. Included seam, signed-branch coverage, native
phase and unchanged original radiation interface accepted. The quarter-cycle
check is an algebra fixture, not independent BEM accuracy evidence.

Coordinator authorized READY-only control update, with all physical inputs and
producer algorithm unchanged. READY controls SHA256:
327d940ab0514eb1d069321cd1c45a78e6e1d209fbd402204bec7088acc67c6b.
Archived `excitation-proposal-ready.json` SHA256:
74aeba4dcf11e3712f0cd6ad96c10b8544bb6c2960d610ff82a8f3c8226d5209.
Actual HEAD4cf9405137cf1e838817cb1b45fa0a173e725df8. Source SHA unchanged.
No numerical allocation or startup/motion approval is implied by READY.

## Executed producer evidence

Coordinator subsequently allocated exactly this single816-call export. Completed
816/816 diffraction solves, zero radiation solves,1CPU,69.3083133000182s. Exact
producer29652/guard30896 launched2026-10-24T22:24:48.424287UTC; guard waited exact
child exit0 and tool session27899 completedexit0 (chunk0012cc). No rerun.

All prepared producer checks PASS: finite full6 tables, H=FK+diffraction, exact
copied seam, original triangles/CG/M, MAT/NPZ readback, all frozen inputs and
dependency identities, unchanged original radiation H5. Independent analytic
incident-pressure FK reduction at every node maxerror5.922200433305882e-11.
Quarter-cycle check remains an algebra fixture, not BEM accuracy evidence.

Stable output `artifacts/reference/calm-buoy-20261024/excitation-attempt1`:

- manifest.json SHA21bbe8489c166d69ce1896aa5d3dbddb4a7a6c324f4c7bb0f88e970ba39f2123.
- excitation.npz SHA31c0b956f81f8819c01dad7dcf6b6d42f3526881ff6569d41671b3be70894b8f.
- excitation.mat SHAf82aefd4d8f249d5f1cf05b418fe97a634bd2be91a1b3d6f267652d5e9fa002a.

Parent folder contains `excitation-attempt1-{claim-readback,idle-readback,
terminal-record}.json` and `excitation-attempt1-run.log`. Log uses PowerShell
UTF16 encoding; initial read-only terminal verification used UTF8 and failed
before any state write. Correct UTF16 reread succeeded; no producer retry or
source change occurred. Warnings preserved, including largest panel3.678m vs
wavelength/8=0.667m at3.4rad/s. No interpolation or mesh accuracy claim.

Claimc400bd and separate verifyf69b22 both exit0; exactterminal then frozen-input
verification and idle-under-lease67fdfb, separate idle-readback5664bc, final
lease-release1a9392 all exit0. No numerical owner remains. Consumer notified of
stable outputs; independent executed integrity review requested but not yet
claimed complete. No incident preflight, body run or video job by this owner.

Stopping point: this allocated export complete. Next action belongs to WEC
physical incoming-train preflight under separate coordinator scheduling.

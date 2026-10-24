# ONE calm-water free-decay control: prospective, not launched

Later status: ONE coordinator-allocated control executed and independently
reviewed PASS. Actual receipt below supersedes the earlier source-only status.

Source preparation only,2026-10-24. Failed Jan90case and whole384component bank
remain immutable. This is a diagnostic control, not a replacement canonical sea.

Same published fork f962103a, original17panel H5/fullM/K/Ainf/residualAr,
native acceleration delay, radiation40s/dt.02, ode4dt.02, six free DOFs.
Initial q=[0,0,0,pi/360,0,0] (0.5degroll), all velocity/radiation memory0.
Body idealized held at this small heel in calm water before0 then released0.
All six applied excitation columns are EXACTZERO for every sample, with native
scalar input identicallyzero. No added damping, restraints, coefficient repair,
input scaling, delay removal or model fitting. Diagnostic changes only sea/IC.

Source identities: run_free_decay_case.m
`e56a8e48e1e59a78567bb5192c05b4ac42842a35ebfdd6e47fee746ec787e85f`;
prepare_free_decay_inputs.m
`d6cbe1aee40007a1e3663bfcd4fd0c4e42fbe59c22f45b89431f65b5e8c615f2`;
wecSimInputFile_free_decay.m
`fbf0f4a834ff2d656be5c7d9ae32fa46d6e9dc6518fd2924ce39ed2483d0833c`.
Builder and proper API case force hook are reused unchanged. Known old engine
output variables are cleared at script entry so failure cannot save stale Jan90
output. Stored files are never cleared or replaced.

Actual q0/v0/memory/fullparameterM/H5/grid/all3001x6 zeroexcitation checks apply.
Native raw60columns are preserved before assertions. Added fullworld Solid
Newton-Euler residual (rotated I and worldomega cross Iomega) and exact native
residual acceleration-delay recurrence checks run on all3001samples; raw error
arrays are saved before their1e-6N/Nm runtime assertions. Physical outcome is
reported from complete histories, not the assertion result alone. Sustained/severe
growth rejects bounded small-perturbation use; transient coupled/non-normal
amplification alone does not prove asymptotic instability. Decay does not establish
passivity or validate the forced case. No automatic next control or remedy.

Prospective archive command (hash/copy only):

```powershell
.venv/Scripts/python.exe reference/wecsim/prepare_free_decay_manifest.py
```

After root allocation, authoritative engine readiness, serial claim and SEPARATE
verification, exact numerical command:

```powershell
.venv/Scripts/python.exe reference/capytaine/geometry_execution_guard.py run --proposal artifacts/reference/generated-bank-consumer-20261024/free-decay-attempt1-source/proposal.json --owner /root/wecsim_reference
```

Child is run_matlab.py run_free_decay_case.m, fresh free-decay-attempt1.log;
output is free-decay-attempt1 under the consumer root. One engine/oneCPU,
expected30–90s, declared600s cap,60s/3001samples. Exact terminal+child proof,
idle-under-lease, separate readback and final release before reporting/next work.
No renderer or automatic retry. Source/readiness review precedes any allocation.

Future forced-case IC topic is NOT implemented here. A freely floating stationary
initialization would need consistent q0/v0, finite-CIC velocity prehistory, native
acceleration-delay state and incident/scattered-field interpretation. Seeding only
particular q/v while resetting memory is inconsistent. It cannot be used to hide
the proven nonpassive operator or suppress an unstable homogeneous mode; any such
proposal follows the present control's evidence and separate authorization.

## Executed receipt and limits

Output `artifacts/reference/generated-bank-consumer-20261024/free-decay-attempt1`.
Exact final rune56a8e48/inputsd6cbe1ae/configfbf0f4a8 and all archived inputs/
dependencies remained unchanged through terminal. Total30.555s (preprocess17.967,
compile2.814,dynamics6.214). All3001native samples retained, all six applied
excitation columns EXACTZERO, known0.5degroll/v0/memory0 correct. Fullworld
Newton-Euler residual<=6.83e-13N/Nm and native delayed-Ar residual<=1.14e-12.

Roll window maxima decay0.5->.03322->.005096->.0001774->.0001730->.0000902deg;
final roll.00004258deg. Maximum translations[surge,sway,heave]=
[.00026933,.01576186,.00046838]m; maximum rotations[roll,pitch,yaw]=
[.00872665,.00123209,.0000655915]rad. Final displacement is
[-.00005075,.01576186,.0000002006]m. Small unrestrained drift is retained.
This one small unforced perturbation shows no severe growing response over60s;
it does not prove passivity, stability for all states, or credibility of the
forced model. Guaranteed full-matrix nonpassive directions remain, and the
failed Jan90forced case is unchanged/rejected.

Independent numerical-methods reviewer verified frozen hashes/rawIC/zero forces,
all balance identities and window decay statistics without a simulation rerun.
Original engine9168, exact client36696/guard30692 terminal0 (session88564,
chunk4f11dd); authoritative both absent, only child mwdocsearch33628 nonnumerical.
Idle recorded under lease, separate readback961a97 passed, release FINAL. No
encoder, retry or next numerical job. Current next assignment is a readable
SOURCE/PHYSICAL IC PROPOSAL ONLY for consistent q/v/radiation/delay prehistory;
no implementation or numerical grant.

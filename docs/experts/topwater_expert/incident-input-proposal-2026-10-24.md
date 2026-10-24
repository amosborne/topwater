# REF-002 incident-input diagnostic proposal, 2026-10-24

Prepared and then executed ONCE under a separate coordinator allocation after
Capytaine's geometry allocation was terminal and its lease released. This note
retains the proposed scope and appends actual execution evidence below.

Producer: `reference/sea_states/incident_input_diagnostics.py`, config
`reference/sea_states/incident-input-diagnostics.json`. Proposed command:

```powershell
.venv/Scripts/python.exe -u reference/sea_states/incident_input_diagnostics.py --config reference/sea_states/incident-input-diagnostics.json --output artifacts/reference/incident-input-diagnostics-20261024
```

One future coordinator-allocated job, one CPU/thread, expected under5s (estimate,
not measured),60s maximum. It imports NumPy only, launches no subprocesses and
calls no generator, BEM, WEC-Sim or Topwater solver. It reads the three immutable
seed42/384-component banks exported earlier and existing0/45/90deg body mappings.
No rephasing, rescaling, time shift, spectrum truncation or new realization.

For each of nine INTERNAL candidate combinations, sample0..60s inclusive at
dt0.02s. Evaluate the exact harmonic sum and analytic eta_x,eta_y,eta_xx,eta_xy,
eta_yy and eta_t. The Hessian entries are graph second derivatives; they are
not principal geometric surface curvatures. Algebraically grouping the eight
equal-frequency directional components avoids repeated temporal evaluation while
retaining their exact amplitudes, directions and phases. Independent direct
ungrouped derivatives at t=0 provide a numerical identity check.

Points are CG and unique projections of every full equilibrium-hull vertex,
triangle centroid and edge midpoint. The baseline quaternion maps material
geometry into initial boat-aligned equilibrium coordinates. Save all fields,
coordinates, times and departures from the instantaneous CG tangent plane in NPZ;
JSON reports every t=0 value, sampled extrema, Hm0/Tm02, global boat yaw and body
mean propagation/arrival plus directional resultant. CG input slope-angle proxies
are exactly rolllike atan(eta_y), pitchlike -atan(eta_x), in degrees. They DO NOT
predict boat response or define the coupled Euler pose of the boat.

Sampled footprint departures are not claimed continuous extrema. Report a
separate conservative all-time continuous-footprint bound
0.5*max_vertex|r-CG|^2*sum(a*k^2), from the spatial Taylor remainder/Hessian norm.
It covers the entire projected triangulated hull since the squared radius is
bounded by its vertices. eta_t is included to expose nonstationary local startup;
no boat velocities or radiation-memory initialization are supplied.

Acceptance: accepted bank manifest/output hashes match; unchanged canonical
physical input; finite fields; direct t=0 sum closure; sampled departure within
continuous Taylor bound; source/input/NumPy binary hashes unchanged; all output
hashes recorded; explicit input-only/initialization limitations. Syntax/config
parse and diff whitespace check passed before launch. Actual execution and
postchecks all pass; no retry was needed.

Stopping point after a separately authorized execution: exact command/process
terminal proof, immutable input and output hashes, energy/direction/t0/plane
diagnostics returned to coordinator, ledger idle while lease held then released.
No canonical selection from slope proxies. Actual WEC producer/response and the
geometry envelope remain necessary for selecting ONE final case.

## Executed evidence

ProposalSHA25643ca2183fb7ab873c20cd88ee309b1a92cf78656902895f0a186c99324294135;
scriptSHA256628d78bf75514d078fe4fe6939e4b8d608f13e6eee1312e6f4515f71d1c3ff70;
configSHA25673757c2096cfb0dd0490a72c2e7fab6f1d2018c8374d23fd2579ddebbe2be4b6.
Sources/proposal/inputs archived under
artifacts/reference/incident-input-diagnostics-executed-source.

Claim creation/identity/archive/ledger command29fe7f exited0; separate claimed
readback67003c exited0. Exact producerPID7120, launcher11396, Popen.wait returned0;
shell chunkfd528b exited0. All9cases complete in2.2994226s,57points each,
3001time samples. No numerical subprocesses. Empty stderr and post-terminal
Get-Process inventoryd9f937 show producer/launcher/Python/ffmpeg absent.
Source/dependency/outputhash checks6984ab pass; idle ledger written under claim,
separate idle readback0d9ce1 passed, lease removed LAST by984f78 exit0.

Output: artifacts/reference/incident-input-diagnostics-20261024/manifest.json,
nine full-field NPZ and nine per-case JSON summaries. Additional run.log, run.err,
terminal.json, postcheck.json, claim-readback.json and idle-readback.json have the
same incident-input-diagnostics prefix under artifacts/reference.

All t=0 direct-ungrouped identity errors <=7.78e-16 in each field's stated units.
All fields finite. Existing measured/realized energy and phase banks unchanged;
Hm0=2.2784205055/2.2761370785/2.2703303724m. Body mean propagation is180/135/90deg
for relative offsets0/45/90deg respectively; body arrival0/315/270deg.
Resultants remain0.759067/0.690983/0.790152. Default offset0 faces mean arrival.

These substantial prescribed surface variations do not establish roll response.
The initial surface is nonflat, has nonzero vertical velocity and differs from a
CG tangent plane across the hull; merely setting a boat to that plane does not
initialize the radiation/fluid state or eliminate curvature. The conservative
continuous all-time Taylor bounds7.243/10.692/3.778m are loose and must not be
presented as observed departures. Observed sampled departures are in the table.

## Numeric input summary

Dates identify buoy records; every row is the same original seed42 bank under
the indicated heading. Angles below are CG WATER-PLANE proxies, not body angles.

|Record date|Offset deg|etaCG0 m|Roll/pitch proxy0 deg|Peak abs roll/pitch proxy deg|Departure0/peak m|
|---|---:|---:|---|---|---|
|20220424|0|-0.0446807|2.35201 / -0.563043|7.03821 / 8.42104|0.0688379 / 0.372246|
|20220424|45|-0.0446807|1.26571 / -2.06131|8.00429 / 9.84977|0.14939 / 0.459679|
|20220424|90|-0.0446807|-0.563043 / -2.35201|8.42104 / 7.03821|0.234485 / 0.495379|
|20250307|0|0.253747|5.42526 / 1.87256|9.38714 / 11.8541|0.136381 / 0.779795|
|20250307|45|0.253747|5.15833 / -2.52155|11.4181 / 10.5548|0.263884 / 0.523397|
|20250307|90|0.253747|1.87256 / -5.42526|11.8541 / 9.38714|0.235968 / 0.420643|
|20250105|0|-0.250573|1.84453 / 0.894968|4.23683 / 5.71852|0.0570408 / 0.280289|
|20250105|45|-0.250573|1.93689 / -0.671813|4.33571 / 4.68687|0.129571 / 0.264743|
|20250105|90|-0.250573|0.894968 / -1.84453|5.71852 / 4.23683|0.183036 / 0.306707|

Initial CG analytic slopes/Hessian (Hessian units1/m); all values refer to
the incident graph at t=0, not the hull or dynamic restoring curvature.

|Record date|Offset deg|eta_x0|eta_y0|eta_xx0|eta_xy0|eta_yy0|eta_t0 m/s|
|---|---:|---:|---:|---:|---:|---:|---:|
|20220424|0|0.00982726|0.0410734|0.00456878|0.00367507|-0.0434162|0.0441355|
|20220424|45|0.0359922|0.0220944|-0.0157486|-0.0239925|-0.0230988|0.0441355|
|20220424|90|0.0410734|-0.00982726|-0.0434162|-0.00367507|0.00456878|0.0441355|
|20250307|0|-0.0326939|0.0949726|-0.0039444|0.00531535|-0.0572906|-0.0840947|
|20250307|45|0.0440377|0.0902739|-0.0253022|-0.0266731|-0.0359329|-0.0840947|
|20250307|90|0.0949726|0.0326939|-0.0572906|-0.00531535|-0.0039444|-0.0840947|
|20250105|0|-0.0156214|0.0322043|-0.0046114|-0.00398065|-0.0248422|-0.0898148|
|20250105|45|0.0117259|0.0338179|-0.0187074|-0.0101154|-0.0107461|-0.0898148|
|20250105|90|0.0322043|0.0156214|-0.0248422|0.00398065|-0.0046114|-0.0898148|

Sampled extrema over all57points and3001times; curvature columns are maximum
absolute Hessian entries, not continuous maxima or principal curvatures.

|Record date|Offset deg|eta min/max m|max abs eta_x/eta_y|max abs eta_xx/eta_xy/eta_yy 1/m|
|---|---:|---|---|---|
|20220424|0|-1.70129 / 1.02067|0.183299 / 0.131398|0.0551707 / 0.0345449 / 0.0664558|
|20220424|45|-1.71459 / 1.05572|0.180407 / 0.145808|0.0713523 / 0.0387911 / 0.0562436|
|20220424|90|-1.68757 / 1.06108|0.131693 / 0.175963|0.0677276 / 0.0330714 / 0.0547892|
|20250307|0|-1.86976 / 1.10174|0.221785 / 0.195169|0.109206 / 0.0536695 / 0.0750042|
|20250307|45|-2.00544 / 1.09607|0.208701 / 0.225995|0.0866879 / 0.050956 / 0.110554|
|20250307|90|-1.92422 / 1.13072|0.181918 / 0.225328|0.0620383 / 0.0663661 / 0.101242|
|20250105|0|-1.5735 / 1.03414|0.135723 / 0.0821336|0.0330953 / 0.0194248 / 0.0391141|
|20250105|45|-1.58717 / 1.04511|0.0965286 / 0.099898|0.0335742 / 0.0246812 / 0.0292045|
|20250105|90|-1.58717 / 1.05415|0.0937634 / 0.121783|0.0403897 / 0.0223005 / 0.0374964|

# Spectral handoff implementation: component qualification

Root's independently reviewed `spectral-handoff-check-attempt1` passed at exact
fork head e24c1a73372baf37aebe37a8b6bf0cd0d2368b0c on Windows/MATLAB R2025b:
38 Python controls (3.53s),14 native MATLAB controls and9 official H5 handoff
controls. `matlab/acceptance.json` reports passed=true; `output.log` and
`terminal.json` retain the execution evidence under
`artifacts/workflow/spectral-handoff-check-attempt1/`. Terminal exit0, frozen
identities and process absence were verified; root released the lease.
The reviewed integration source manifest remains
`native-spectral-integration-source-v1/manifest-v2.json`, SHAaa154432.

This accepts the tested component interface, including optional-view low/off-grid
node handling, supplied-phase RNG isolation, legacy fallback and official H5
raw/radiation preservation. It does not accept full simple-boat forces, startup,
actual applied bus, motion or physical accuracy. Those checks remain pending;
the generic native-bus evidence is recorded separately in the FullDir checkpoint.
Held runtime995 and reference/pins remain unchanged. No author imports, tests or
MATLAB jobs were run. Root owns scheduling, integration and publication.

`_wecsim_spectral.prepare_spectral_inputs(setup,hydro)` returns physical FullDir
wave data and a separate normalized `component_view`. It coherently adds only
identical physical headings, preserving both signed branches. The existing Hz
file input remains authoritative; derived view frequencies reproduce native
`(frequencies*2)*pi`, with original packet frequencies/deltas retained as
diagnostics. No new normative wave-frequency field was added. C's packet weights
are already present; native dOmega/dTheta storage factors are compensated once,
with phase=-arg(C). Zero-energy rows stay zero without random phases or floors.

The derived view samples the old bilinear complex H plus the declared quadratic
0..0.1 H0 bridge at all actual native component frequencies and physical heading
nodes. It includes H0, genuine upper coverage boundary and exact periodic
endpoints. RE=real(H)/(rho*g), IM=-imag(H)/(rho*g), required md=0. Raw H5 grids,
RE/IM/FK/SC, A/B, full radiation kernel/time/frequency, Ainf and physical metadata
are unchanged. Python never multiplies C by H or creates an applied force history.

MATLAB preparation accepts `prepared.component_view`, returns `[h,nativeHydro]`,
omits excitation IRF construction/requirements on that branch, verifies official
H5 readback and attaches the view only to the returned public body struct.
Legacy IRF preparation remains available for preserved studies/tests. The
separately owned generic writer change omits excitation-IRF datasets only when
all three fields are absent; a partial set fails before output truncation.

Source ownership is split: BEMIO owner owns Python data preparation, MATLAB
preparation/validation and focused tests; the other WEC owner owns native
optional-view access/coverage, writer, backend/input/model/startup and actual
native force/application path. Contract: view fields omega_rad_s/headings_deg
and re/im/md [6,heading,frequency]; wave fields frequencies/directions/spectrum/
spread/phase. Backend supplies prepared.spectral_file and prepared.component_view.

Keep commit purpose clear: PR14's phase import/direction-width fixes are generic
upstream candidates with their existing native witnesses. Optional view and
writer support are generic interoperability extensions, requiring their own
focused review within the coordinated issue12/PR14 process. Sampling the old
operator, C-to-FullDir storage and canonical gates are Topwater integration policy;
they are not additional upstream bug claims. No upstream submission, new issue,
pin change or merge is authorized by this note.

Next: the coordinator's bounded resource probe, then independently executed
full3001x6 canonical force/startup/application gates, followed by motion/output
validation. Component acceptance does not authorize those jobs or promote the
candidate. Preserve all failed studies and original accepted public outputs.

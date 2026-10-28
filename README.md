# Topwater

Topwater is developing a time-domain potential-flow solver for a freely drifting
boat with six rigid-body motions and a changing wetted surface.

**Topwater's own fluid-body solver is not validated yet.** The selected
implementation is the existing **P0 Rankine-panel core**, extracted from
`4e8f511`: physical inputs, calm-water hydrostatics, Laplace boundary solves and
instantaneous added mass. It has **no coupled time-domain wave-response path**.
The archived nonlinear diagnostic runtime is a different, unqualified experiment.

The canonical60s workflow is complete in [examples/simple-boat](examples/simple-boat/):
boat and sea inputs, six-DOF motion and total free-surface geometry in `result.npz`,
`result.png`, `video.mp4`, and one log for each implemented CLI stage. Capytaine
provides hydrodynamic coefficients; WEC-Sim computes wave excitation and motion
through its existing public inputs. Topwater reconstructs the total incident,
diffracted and radiated water surface from the saved motion without feedback.

The maintained WEC-Sim reference is based on official v7.1.0. The standard-input
replacement matches the original canonical trajectory within0.0598mm position
and0.00157degree orientation. These are compatibility results, not physical-error
bounds. Surface reconstruction now meets the selected-case 1 cm numerical target
under its stated reference assumptions. Mesh resolution, field/load consistency
and linear applicability remain open in [the warning triage](docs/cli/warning-triage.md).
An earlier ten-second minimum-mesh reference remains available for
[historical reproduction](reference/wecsim/README.md).

Native solver implementation remains held. See [the current handoff](docs/checkpoint-handoff.md)
for authorization and [AGENTS.md](AGENTS.md) for startup and expert routing.
The broader [Topwater comparison/runtime goal](docs/time-domain-baseline-goal.md)
remains a later objective. [Issue #30](https://github.com/amosborne/topwater/issues/30)
tracks checkpoints.

## Canonical boat and sea

![Canonical 24 ft boat at calm-water equilibrium](docs/images/reference-boat.png)

The single application fixture is a simple **24 ft by 9 ft** open V-hull with
**2,500 kg** total mass: 1,600 kg distributed hull/structure, a 600 kg amidships
point load and a 300 kg stern motor point load. Integration determines its
centre of mass and full inertia tensor. See [the boat definition](docs/reference-boat.md).

The current canonical inputs are [boat.npz](examples/simple-boat/boat.npz) and
[sea.npz](examples/simple-boat/sea.npz): the same physical boat and the actual
46086/2023-12-23T00:40Z spectrum with seed42. One example is one boat/sea pair.
Both are plain numeric archives; [the NPZ contract](docs/npz-format.md) defines
their units, frames and arrays. They do not initialize or qualify native physics.

The intended range is convex 18-30 ft boats, 7-11 ft beam, in ordinary ocean
conditions. Incident wavelengths must be at least twice beam and frequencies
at most 0.6 Hz. Breaking, slamming, flooding and loss of contact are outside the
supported simulation scope. Changing supported boats must require input data
only; the reference mesh is not a solver requirement.

## WEC-Sim linear approximation

The incoming wave field is specified independently of the boat. Capytaine
calculates hydrodynamic responses about a fixed mean hull: its equilibrium
position, orientation, draft and waterline. The boat itself is free to move in
all six degrees of freedom during the WEC-Sim simulation.

Incoming-wave excitation uses those fixed responses; moving the boat does not
recompute its hydrodynamic geometry or diffraction solution. Motion still feeds
back through displacement-dependent restoring forces, added water inertia and
velocity-history-dependent radiation forces. These hydrodynamic responses retain
the mean-hull linearization.

The reconstructed diffraction and radiation fields also remain referenced to
the mean hull. Diffraction follows the incoming waves, while radiation is driven
by the computed motion. These disturbance fields are not rigidly carried along
with the animated boat. The moving waterline is obtained afterward by intersecting
the reconstructed total water surface with the current posed hull; it does not
feed back into the forces.

Motion can expose water locations inside the original mean-hull waterline
footprint, where the exterior hydrodynamic solution is not defined. Today we
approximate the disturbance at each such location using the nearest mean-waterline
edge: find its closest boundary point, step outward by
`min(0.025 m, 0.01 * beam)`, and copy the diffraction and radiation values from
that exterior sample. Nearby exterior locations within this offset use the same
trace policy. The incident wave is evaluated at the actual requested location.
Thus only the disturbance is continued, approximately constant inward along a
nearest-edge normal; the total water elevation is not held constant. The offset
specifies where to sample outside the hull, not a maximum continuation distance.
No blending or continuity guarantee is imposed where the nearest edge changes.

This continuation does not solve a moving fluid boundary. Separation between
the fixed mean waterline and moving waterline is a useful visual indication of
departures from the assumed wetted geometry, rather than an accuracy estimate:
small separation alone does not establish linear-model accuracy. The numerical
reconstruction check is recorded in [issue #31](https://github.com/amosborne/topwater/issues/31);
fixed-mean applicability remains open in [issue #34](https://github.com/amosborne/topwater/issues/34).

The reconstruction target is **1 cm maximum numerical elevation error against
this chosen linearized model**, including reference uncertainty. We hold the
hydrodynamic data and solved motion fixed, calculate a separate reference that
retains the full complex response, and compare elevations at the same published
grid points and times. Numerical resolution checks and an allowance for omitted
frequencies contribute to the reference uncertainty. The difference plus that
uncertainty must fit within 1 cm. This check does not establish that the linear
model predicts nonlinear calculations or real water within 1 cm.

For the earlier selected 60-second case, all 5929 surface targets and 751 saved frames
pass: maximum difference is 0.341 mm, and maximum difference plus reference
uncertainty is 5.784 mm. This relies on a stated high-frequency response envelope;
see [the qualification and limitations](docs/experts/numerical_methods_researcher/issue31-selected-case-acceptance-2026-10-27.md).
That bound does not automatically qualify the new 130-second packet case.

## Current foundation

`HullMesh`, `RigidBody`, `WaveComponents`, `InitialState` and `SimulationInputs`
own validated, immutable NumPy data in SI units. Body axes are bow/port/up;
orientation is an active body-to-world quaternion `(w, x, y, z)`. Position is
world CG, linear velocity is in world axes and angular velocity is in body axes.
Hull faces wind outward from the body. Artificial opening caps carry neither
shell mass nor water pressure.

`topwater.hydrostatics` computes flat-water buoyancy and restoring derivatives.
`topwater.bem` and `topwater.boundary` assemble and solve Laplace problems.
`topwater.inertia.instantaneous_added_mass` computes six instantaneous responses
for a flat zero-potential surface. This is not a finite-frequency wave solver.

## Canonical workflows

All commands take an optional example directory, defaulting to the current
directory. Filenames are fixed:

```powershell
uv run topwater generate-sea examples/simple-boat --record ndbc-46086-20231223T0040Z --seed 42
uv run topwater solve examples/simple-boat 100 --rest 10 --ramp 5 --solver=wecsim
uv run topwater render examples/simple-boat
```

The requested 100 s is dwell: 10 s rest, 5 s ramp, 100 s dwell, 5 s ramp,
10 s rest give 130 s total. Packet construction and WEC-Sim use a fixed 10 ms
timestep. Sea phases refer directly to the dwell midpoint, 65 s here. Rest
defaults to 10 s and ramp to 5 s; the end is fixed, without automatic settling
extension. See [packet construction](docs/wave-packet-search.md).

Generation writes sea.npz. Solve reads boat.npz/sea.npz and writes the complete
result.npz and six-DOF motion plot result.png; an existing video is left untouched.
Render reads only result.npz,
including its exported hull snapshot, and writes video.mp4 without recomputing
waves. Results permit independently sampled motion and changing surface topology.
Validated data/video outputs are replaced atomically; failures preserve the
previous data. Commands also write consistently formatted `sea.log`, `result.log`
and `video.log`, each describing its latest attempt. Logs are replaced at command
start, including failed attempts. See [logging](docs/cli/logging.md) and
[warning triage](docs/cli/warning-triage.md).
Final result/video binaries use Git LFS. No intermediate solver files belong in
the example directory.

Read each command's inputs, hard-coded choices, errors and warnings:
[generate-sea](docs/cli/generate-sea.md), [solve](docs/cli/solve.md),
[render](docs/cli/render.md). The renderer retains the accepted800square
Matplotlib top-down style. Historical incident-only movies and measurements
remain documented under [the reference workflow](reference/wecsim/README.md).

Native remains the default and reports unavailable. To prepare this machine for
the first WEC-Sim backend, use `uv sync --locked --group wecsim`; MATLABR2025b,
licensed toolboxes and the maintained fork are additional requirements. Current
surface qualification limits are stated on the solve page rather than silently
substituting incident waves.

## Development

```powershell
uv sync --locked
uv run --locked pytest
```

Capytaine supplies coefficients and field reconstruction for the optional
WEC-Sim workflow. It also supplies independent comparison data for native-core
validation:

```powershell
uv sync --locked --group correlation
```

Generic mathematical tests do not establish Capytaine agreement. The frozen
[validation contract](validation/canonical-boat.json), [reset instructions](docs/validation-reset.md)
describe historical validation requirements; current scope is in the handoff.

The later requirement is 0.1 degree in each rotation and 1 cm in each
translation throughout 60 seconds, at no more than 20 wall seconds per simulated
second on the development machine. Real time is preferred. Neither this longer
target nor the first one-second/200x milestone has been demonstrated.

Optimization history is preserved at Git tag `archive/optimization-2026-10-23`.
The retired canoe and historical experimental paths do not qualify this boat.

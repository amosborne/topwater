# Topwater

Topwater is developing a time-domain potential-flow solver for a freely drifting
boat with six rigid-body motions and a changing wetted surface.

**Topwater's own fluid-body solver is not validated yet.** The selected
implementation is the existing **P0 Rankine-panel core**, extracted from
`4e8f511`: physical inputs, calm-water hydrostatics, Laplace boundary solves and
instantaneous added mass. It has **no coupled time-domain wave-response path**.
The archived nonlinear diagnostic runtime is a different, unqualified experiment.

The first minimum-mesh Capytaine -> official BEMIO -> WEC-Sim reference is
complete: ten seconds of six-DOF motion from rest, a shared upstream incident
packet and a centered isometric video. The17-panel,0.02s fixed-step result has
independent input/interface/physics/video reviews; temporal and mesh accuracy
remain unestablished. [Reproduction and artifacts](reference/wecsim/README.md).

The active follow-up is ONE canonical case from the collected-buoy sea-state
generator, targeting about5-10deg roll and appreciable motion in the other5DOF.
Its linear operating range will be assumed a priori from wetted-surface and
buoyancy/righting curves against WEC-Sim's fixed-mean assumptions. This work
does not require a new nonlinear trajectory model. Topwater solver implementation
and performance work remain held. See [the current handoff](docs/checkpoint-handoff.md)
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

[The shared case](tests/reference_case.py) combines this boat at calm-water
equilibrium with the [observed Southern California sea](docs/sea-states.md),
using phase seed 42. These are physical inputs, not an initialized or validated
nonlinear fluid simulation. Linear reference cases use the same boat with
explicitly declared small waves and physical constants.

The intended range is convex 18-30 ft boats, 7-11 ft beam, in ordinary ocean
conditions. Incident wavelengths must be at least twice beam and frequencies
at most 0.6 Hz. Breaking, slamming, flooding and loss of contact are outside the
supported simulation scope. Changing supported boats must require input data
only; the reference mesh is not a solver requirement.

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

## Development

```powershell
uv sync --locked
uv run --locked pytest
```

Capytaine is a validation dependency, not a runtime load model:

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

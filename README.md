# Topwater

Topwater is developing a time-domain potential-flow solver for a freely drifting
boat with six rigid-body motions and a changing wetted surface.

**There is no validated fluid-body simulation baseline yet.** The selected
implementation is the existing **P0 Rankine-panel core**, extracted from
`4e8f511`: physical inputs, calm-water hydrostatics, Laplace boundary solves and
instantaneous added mass. It has **no coupled time-domain wave-response path**.
The archived nonlinear diagnostic runtime is a different, unqualified experiment.

A validated, practical **time-domain** baseline is the immediate goal: the
canonical boat in a frozen mild sea, at least one simulated second, matching a
Capytaine-derived reference within 0.1 degree and 1 cm in every axis, including
uncertainty, and costing at most 200 wall seconds per simulated second including
normal Topwater preparation. Major optimization follows a passing physical
comparison; early cost screens reject impractical routes. See the
[current goal and checkpoints](docs/time-domain-baseline-goal.md).
[Issue #30](https://github.com/amosborne/topwater/issues/30) owns promotion
requirements and the concise investigation register.

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
and [AGENTS.md](AGENTS.md) define the validation-first work.

The later requirement is 0.1 degree in each rotation and 1 cm in each
translation throughout 60 seconds, at no more than 20 wall seconds per simulated
second on the development machine. Real time is preferred. Neither this longer
target nor the first one-second/200x milestone has been demonstrated.

Optimization history is preserved at Git tag `archive/optimization-2026-10-23`.
The retired canoe and historical experimental paths do not qualify this boat.

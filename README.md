# Topwater

**Nonlinear wave-body dynamics in Python.**

Topwater is a planned open-source, three-dimensional time-domain potential-flow solver for freely floating bodies in directional waves. The initial target is one rigid boat with six free motions, an evolving free surface and changing wetted hull geometry, explicit viscous damping, and an external force and torque interface.

**Status: planning and numerical verification infrastructure. No solver is implemented yet.**

## Approach

Solve incompressible, inviscid, irrotational water flow using boundary elements on the wetted hull and surrounding free surface. Couple pressure and body acceleration consistently, so buoyancy, wave excitation, diffraction, radiation and added-mass effects emerge from the fluid solution rather than being added from fixed coefficient tables.

Use Python with NumPy/SciPy, then optimize measured bottlenecks with compiled kernels. Existing Capytaine results can supply independent linear benchmarks.

The initial scope is a single rigid body in deep water, with reproducible directional wave realizations and no anchor or mooring forces. Viscous damping remains an explicit model requiring sensitivity analysis or calibration. Breaking/overturning waves, resolved slamming jets, flooding, propulsion and multiple bodies are outside the initial scope. Free yaw does not guarantee a settled heading.

## Roadmap

Requirements, implementation stages and acceptance criteria live in [GitHub issues](https://github.com/amosborne/topwater/issues/1). Each stage requires verification evidence and measured execution time before later stages rely on it.

- [ ] [Stage 0: Define solver contracts, physical scope and verification gates](https://github.com/amosborne/topwater/issues/2)
- [ ] [Stage 1: Implement and verify the boundary-element numerical core](https://github.com/amosborne/topwater/issues/3)
- [ ] [Stage 2: Propagate nonlinear directional waves without a boat](https://github.com/amosborne/topwater/issues/4)
- [ ] [Stage 3: Track the moving waterline and instantaneous wetted hull](https://github.com/amosborne/topwater/issues/5)
- [ ] [Stage 4: Recover consistent pressure and verify prescribed-body hydrodynamics](https://github.com/amosborne/topwater/issues/6)
- [ ] [Stage 5: Couple the fluid to six freely moving rigid-body degrees of freedom](https://github.com/amosborne/topwater/issues/7)
- [ ] [Stage 6: Add explicit viscous damping and external force and torque interfaces](https://github.com/amosborne/topwater/issues/8)
- [ ] [Stage 7: Qualify one boat in one observed directional sea state](https://github.com/amosborne/topwater/issues/9)
- [ ] [Stage 8: Profile, optimize and publish a reproducible release](https://github.com/amosborne/topwater/issues/10)

The first end-to-end demonstration will cover one boat and one observed directional sea state. Broader boat catalogs and sea-state populations follow only after the numerical and runtime checks justify them.

## Contributing

Start with an open stage issue and discuss changes that affect its physical assumptions or interfaces. Contributions should distinguish numerical verification from physical validation, preserve wave/input provenance, and include timing and convergence evidence appropriate to the change. Do not commit generated simulation output or third-party code without checking its license.

## License

[MIT](LICENSE), copyright 2026 Alex Osborne.

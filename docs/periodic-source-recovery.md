# Periodic source recovery 003

The bounded source-recovery stage passed at source checkpoint `785c3d1`:
all 84 Capytaine solves and all 384 checks completed. The elapsed 329.246 s is
offline reference generation, not a Topwater simulation or its runtime ratio.
No Topwater runtime Python changed; no coupled time-domain run exists yet.

The [prospective protocol](periodic-source-protocol.md) and
[controls](../validation/periodic-source-controls.json) were committed before
launch. The four unchanged 1,232/2,772/4,928/7,700-panel mean-hull meshes use
the frozen canonical boat, world CG, full mass/restoring matrices, density 1025,
gravity 9.81, infinite depth and zero speed. Three frequencies each retain all
six unit-displacement radiation modes and the selected unit-elevation diffraction
mode. Capytaine 3.0's indirect Delhommeau source density cannot be interchanged
with a density for Topwater's Rankine operator.

Reproduce with the locked correlation environment and a fresh output stem:

```powershell
uv run --locked --no-sync --cache-dir artifacts/uv-cache --group correlation python -m tests.validation.periodic_sources --raw artifacts/validation/reference-study-003.json --raw-bundle validation/evidence/reference-study-003.zip --frozen artifacts/validation/frozen-fixed-001.json --output artifacts/validation/periodic-sources-next.json
```

The raw/frozen artifacts can be recovered from the nested complete evidence
bundles. Do not overwrite old stems or run alongside another numerical job.

Each basis retains complex source density, hull potential/pressure, prescribed
normal velocity, and the potential/normal velocity reproduced by the original
S/K operators. Capture occurs before Capytaine replaces its cached K with an LU
factorization; the non-overwriting default LU is required. Independent pressure
integration checks body-force signs and world-CG moments. All 36 finite A/B
entries and the selected complex diffraction/Froude-Krylov loads reproduce the
pinned study 003 raw data at the prospectively declared tolerances. All original
raw Ainf/A/B/excitation/RAOs are retained; study 003's failed decision is unchanged.

The complex periodic composition retains incident amplitude and phase and the
six coupled body-response phasors. Its total hull normal velocity includes the
analytic incident velocity, and is checked against independently derived rigid
normal modes. Total hull pressure includes incident, diffraction and radiation
contributions and reproduces their combined loads. The largest check/tolerance
ratio was 0.0405471. These are reproduction checks, not field convergence or a
measured physical-error bound.

[Complete evidence bundle](../validation/evidence/periodic-sources-003.zip),
[summary/SHA256](../validation/evidence/periodic-sources-003-summary.json).
SHA256: `4a317cbb91cc51348b0059b40b4a6437315e3656624ad2903743dec4721214e5`.
The bundle includes the exact successful producers/inputs/data, the complete
nested raw-reference/freeze bundle, and both unsuccessful attempts:

- 001: interrupted during native initialization to correct review findings;
  zero completed physics calls, exact session 50564 terminal exit 1.
- 002: symbolic Capytaine rigid-motion objects were treated as sampled arrays;
  zero completed physics calls, exact session 69350 terminal exit 1. Corrected
  geometry-derived modes match Capytaine's evaluated modes in a no-Green API test.
- 003: exact session 32728 terminal exit 0, all calls/checks complete.

The full local suite passed 238 tests after completion (session 2470, exit 0).
The optional installed-Capytaine API test skips when that dependency is absent
in generic CI; the actual recovery command still fails without it. This test-only
post-run guard is the sole Python change after the archived numerical producer.

Independent adversarial review verified all 84 calls and 384 unique checks,
archived/current numerical identities, pinned case/protocol/raw bundle, exact
four meshes and every finite basis/composed array. The reviewer independently
recomputed A/B, diffraction/FK and combined loads from pressure, analytic
incident fields, and phase/normal-velocity compositions. Bundle SHA and the
false qualification/readiness flags were verified. No blocking finding remains
for source recovery; this approval does not extend to field/restart qualification.

This stage exports no nonlinear-grid eta0/phi0 packet. Full reference qualification
and nonlinear restart readiness remain explicitly false. Field spatial/mapping
uncertainty, nonlinear boundary compatibility, projection transients and approach
to the linear limit remain unresolved. Future source/pressure/trajectory histories
must never drive Topwater; only physical initial fields and prescribed incident
waves can enter its own time-domain model.

Next bounded work: freeze a minimum fluid grid/domain and admissible mean-fluid
to nonlinear-domain map; reject invalid reference points and quantify the field
and mapping uncertainty. Then build and check the compatible initial packet.
Screen actual coupled preparation/RHS cost before a one-second run; reject a
route projected above 200 wall s per simulated s. Major optimization remains held
until a physical comparison passes, and neither accuracy nor cost is qualified.

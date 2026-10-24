# Canonical Capytaine / BEMIO / WEC-Sim reference

Read ../../docs/capytaine-wecsim-workflow.md first. The active workflow uses
Capytaine native outputs, official BEMIO conversion and WEC-Sim for the
six-DOF10s time-domain simulation. No Topwater imports or solver edits.

canonical-boat.json is neutral shared geometry/mass/sea input. Its spatial
arrival input remains explicitly prospective, to be frozen before simulation.
The boat starts at rest in locally calm water; waves roll in from upstream.
No temporal ramp/global force step. Startup need not settle within10s.

history/ preserves the abandoned generic Python time-domain prototype. It
passed four mathematical checks but is NOT used for the approved workflow.
Minimum-mesh end-to-end motion/video comes first, then mesh variation is judged
on those full SIX motion histories. Intermediate coefficients are diagnostics.

Pinned WEC-Sim v7.1.0: a23d2993807410fe78cc63546c1ba9d81912024f.
Download: git clone --depth 1 --branch v7.1.0 https://github.com/WEC-Sim/WEC-Sim.git artifacts/dependencies/wec-sim-v7.1.0
Verified setup is recorded in ../../docs/capytaine-wecsim-environment.json.

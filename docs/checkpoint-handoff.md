# Checkpoint handoff

- **Active user outcome:** Topwater is TIME-DOMAIN ONLY. Deliver a small-motion
  nonlinear six-DOF simulation of the canonical 24 ft boat in a frozen mild sea,
  at least 1 s, matching a Capytaine-based time history within 0.1 degree per
  rotation and 1 cm per translation including uncertainty, and practical at
  <=200 wall seconds per simulated second including normal Topwater preparation.
  Read docs/time-domain-baseline-goal.md. The reference-only goal is a prerequisite,
  not the endpoint. These latest user instructions supersede older scope holds.
- **Priority and speed:** validation before major optimization. Measure setup/a
  few steps early and reject obviously impractical routes; do not spend indefinitely
  on a 5000x model. Bounded feasibility fixes are allowed, but major architecture
  or speed campaigns wait for an actual passing coupled comparison. Both accuracy
  and practical runtime are required for final completion. No result exists yet.
- **Current implementation:** same unqualified P0 Rankine-panel core extracted
  from 4e8f511, with hydrostatics and instantaneous inertia. It still has no coupled
  TD wave-response path. The archived P1 diagnostic runtime is a different
  implementation. Reuse sound math/tests selectively, never substitute its results.
- **Branch/source:** baseline/fluid-body. Study 003 producer commit 2de87a8, numerical
  family source 9692544; full commits/hashes in manifests. Inspect live git status/log.
  This checkpoint changes goal documentation/evidence only, not numerical Python
  producers. No PR was created; main was not pushed. Preserve other worktrees.
- **Study 003 completed NEGATIVE:** all four meshes/all 132 calls completed; every
  original coefficient/trend family failed. All six Ainf diagonal latest increments
  also fail (1.93-5.81 times budget), so failure is not solely tiny couplings.
  Full coverage, meshes, source/input identity, signs, excitation/RAO construction
  and every gate/uncertainty were independently audited. No sign/unit defect found.
  Original decision stays unqualified; do not retroactively accept or overwrite it.
- **Evidence:** docs/reference-study-003.md, validation/evidence/reference-study-003.zip
  and companion summary. Includes raw matrices, producer archive, readiness/memory
  observations and nested bound-freeze.zip. Its SHA256 is explicitly bound and
  independently verified. Prior studies 001/002 and diagnostics remain preserved.
- **Selected spatial data:** 308-panel fixed graded seed, affine transfinite N=2/3/4/5,
  counts 1,232/2,772/4,928/7,700, farthest-first vertex order, exact canonical surface.
  Frozen-fixed-001 is still source-current; complete data/proof are under
  validation/evidence/frozen-fixed-001.*. Geometry passes, hydrodynamics unqualified.
  Old reference policy 0.1%+1e-6 normalized and all 36 entries remain preserved.
  The prospective time-history protocol must be explicit/reviewed before judging
  a new mild case; no silent gate relaxation or automatic blanket mesh refinement.
- **Interpretation limit:** existing nine steady single-wave RAOs at 0.01 m amplitude
  show empirical mesh increments <=1.2883e-5 m / 0.001703 degree. They do not bound
  continuum, frequency band/tails/spacing, memory, startup or transient error, and
  are not a new mild-sea validation. Use output-budget analysis to select necessary
  reference work, not to pretend the old failed coefficient gates passed.
- **Recovered prior conversion work:** git show b568edd:tests/capytaine_time_reference.py;
  analytic tests at 135c430:tests/test_capytaine_time_reference.py; old exporter at
  b568edd:tests/experiments/capytaine_time_data.py. Cummins kernel, causal convolution
  and back-transform signs/endpoints are sound, but it is radiation math for supplied
  V/a histories, not a coupled motion/excitation solver. Its depth 4 m/1000 kg/m3 canoe
  pin is inapplicable to this deep-water/1025 kg/m3 boat. Preserve generic tests;
  do not restore old application runners or merge entire historical branches.
- **Next bounded action:** choose/freeze a mild, nontrivial sea and >=1 s comparison
  window, exact boat/body/fluid initial state and cost accounting; declare a new
  prospective time-history/uncertainty contract. No phase/amplitude fitting after
  comparison or quiet startup as validation. A scaled catalog case must say so.
  Current catalog's lowest Hm0 is about 0.96 m; that label alone does not establish
  small motion. Use amplitude/representation checks; chine submergence is only
  about 2.8 mm at equilibrium. Report which DOFs are meaningfully excited.
- **Reference choice:** steady harmonic histories need matching periodic body AND
  fluid/history initialization. Rest-start transient histories need validated
  broadband radiation/excitation and band/tail/spacing/memory/time controls. The
  current three frequencies cannot automatically supply them. Reuse prior math
  with current physical identity; Capytaine data remains outside Topwater runtime.
  Complete minimum missing TD physics and early measured cost screens toward the
  first real comparison. No frequency-domain Topwater substitute, imported Capytaine
  loads on top of pressure, or new alternative-formulation/optimization campaign.
- **Live jobs:** none. Study session 73213 terminal exit 1; memory observer 45767 terminal
  exit 0, exact PID 25868 observed through exit. Reference generation 448.064 s is
  offline validation work, not a Topwater speed result. Observed peak committed
  memory 4,453,543,936 bytes. No own-model TD timing has been measured yet.
- **Stopping checkpoints:** resolve case/reference uncertainty and physical failures,
  retain complete evidence, run jobs serially with immutable producers, confirm
  exact handles terminal before edits, obtain independent review, commit/push and
  update issue #30. A negative checkpoint leaves the goal active. Do not claim the
  one-second case establishes the later 60 s/20x target or general nonlinear accuracy.
- **CI:** PRs targeting main only; continue local tests. 51 gate/geometry tests passed
  at prior source checkpoint. No code changes in this documentation/evidence turn
  require rerunning unchanged tests. Never merge/squash-merge without all required
  current-head CI passing.

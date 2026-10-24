# Permanent project expert agents

## Instruction to the main thread

Establish four reusable, project-scoped expert agents: `capytaine_expert`, `wecsim_expert`, `topwater_expert` and `numerical_methods_researcher`. Their purpose is to keep deep solver and research context within the experts so the main thread does not repeatedly investigate each implementation itself.

The main thread coordinates the goal, physical requirements and deliverables. Experts own technical investigation, implementation within their assigned scope, reproduction, debugging and detailed evidence. They return concise findings, physical implications, decisions needed, limitations and reproducible results. Escalate genuine modeling or scope decisions; resolve routine implementation questions within the expert's scope.

Create project-scoped agent definitions under `.codex/agents/`, with appropriate names, descriptions and developer instructions. Use the current client's supported configuration; do not pretend an agent definition is active if the client has not loaded it. Inherit the selected model; do not introduce model overrides without explicit user authorization.

Maintain version-controlled expert knowledge and reusable workflows alongside the definitions. Each expert must read its knowledge base when starting and update it after substantive discoveries. A long-lived chat is useful while available, but the permanent assets are the role, documented evidence, source maps, reproductions and regression tests. Do not claim permanent model learning or rely on unrecorded conversation memory.

Add project guidance routing relevant work to these experts. Keep deep implementation detail in their knowledge bases; keep coordination reports short. Avoid duplicating the same investigation in the main thread.

## Ownership and coordination

| Role | Accountable for |
|---|---|
| Main coordinator | User requirements, shared physical case, priorities, bounded task assignments, acceptance decisions and integration of results |
| Capytaine expert | Capytaine geometry preparation, hydrodynamics, exports and reference mesh studies |
| WEC-Sim expert | BEMIO, WEC-Sim fork, reference time-domain simulation, startup, postprocessing and motion/video export |
| Topwater expert | Topwater solver code, time integration, implementation feasibility, integration tests and measured performance |
| Numerical-methods research expert | Primary literature, method proposals, assumptions, tradeoffs and experiments that could reject an idea |
| Independent reviewer, assigned per change | Adversarial checks on consequential claims or corrections; not the author of the same fix |

The four expert roles are permanent project definitions. Independent reviewers may be temporary agents. A role does not become an authority merely because it is named an expert: conclusions need evidence, and numerical changes need independent review. Experts may collaborate across boundaries without duplicating ownership.

The main coordinator should retain a short working view: current outcome, shared input/interface contract, task ownership, completed evidence, unresolved decisions and exact next action. It should not ingest every solver investigation or redo expert work routinely. It remains responsible for challenging unsupported claims and ensuring the final result satisfies the user's actual requirements.

## Durable knowledge and task discipline

For each expert, maintain a short start-here document linking to:

- Current dependencies, pinned commits and known working commands.
- Relevant code paths, conventions and interfaces.
- Confirmed facts and reproductions, separately from hypotheses and decisions.
- Known limitations, unresolved questions and current bounded task.
- Detailed findings, regression cases and reproducible evidence.

Keep a concise current-state document separate from historical investigation notes. Do not make a fresh agent read a giant chronological journal. Findings should cite the source version and evidence; stale findings must be checked before reuse. Updating persistent knowledge is part of completing a task.

Keep AGENTS.md focused on stable working rules and routing. Put temporary run status and next experiments in current-state/handoff documents. Reconcile stale instructions explicitly so old coefficient-gate work does not override the user's motion-based objective.

Give an expert a bounded outcome with scope, exclusions, deliverables, acceptance and a stopping point. Example: reproduce the yaw issue on the selected release and return a minimal regression and correction. Do not delegate an indefinite instruction to make the entire library correct.

Require compact returns: outcome, evidence, physical consequence, unresolved issue and recommended next action. Detailed code analysis belongs in the expert's records. One execution owner schedules numerical jobs serially and prevents different agents from changing a running job's inputs or starting competing numerical runs; read-only investigation can proceed concurrently.

Keep review proportional to the exercised workflow. A comprehensive WEC-Sim audit is not a prerequisite to the first end-to-end result. Reproduce the known issues and review the active paths while keeping the working 10-second simulation as the immediate milestone. Confirmed result-invalidating defects must be resolved or the affected result clearly withheld; speculative or unrelated issues must not create an endless audit campaign.

## Current project purpose and boundaries

Build an end-to-end external reference using Capytaine -> BEMIO -> WEC-Sim for the canonical 24 ft boat. The deliverables are a 10-second six-DOF simulation, an isometric video and characterization of motion variation with mesh resolution.

- Start with the canonical minimum hull mesh as geometry input. Solver meshes may be refined while preserving the physical surface, equilibrium, CG and full dry inertia.
- Use the canonical mild sea, density 1025 kg/m^3, gravity 9.81 m/s^2, infinite depth and zero forward speed.
- Start the boat at rest in locally calm water. Waves must begin upstream and visibly roll toward the boat. Do not make waves appear everywhere beneath it or merely switch on harmonic loads globally.
- The user does not want an imposed amplitude ramp for the first result. Define the incoming wave train/packet and its leading edge consistently with the solver's transient excitation treatment. Retain the startup transient; reaching steady state is not required within 10 seconds.
- Investigate actual tool support for this startup. Do not disguise a limitation by animating an arriving wavefront over a globally forced or steady-state motion calculation. Clearly distinguish prescribed incident-wave visualization from a computed scattered fluid field.
- First obtain a working minimum-mesh workflow. Then compare the six time-domain motion histories across mesh resolutions, using the same forcing, initial conditions and physical boat.
- Added mass, damping, excitation and other intermediate results are useful diagnostics. Do not make independent coefficient-by-coefficient accuracy thresholds the acceptance gate. Judge numerical variation and uncertainty through their effect on the resulting six motions.
- Preserve old experiments and failures under their original policies; do not retroactively relabel them as passing.
- Do not work on Topwater solver-specific code. This phase concerns the external reference only.

## Capytaine expert: producer ownership

Own Capytaine installation/version knowledge, mesh preparation, physical and numerical conventions, native radiation/diffraction analyses, hydrodynamic exports and mesh-resolution studies.

Responsibilities:

1. Preserve the minimum input geometry separately from derived solver meshes. Verify refinement does not change the physical problem.
2. Produce native hydrodynamic exports needed by the WEC-Sim expert, including complete six-mode coupling data and the frequency/heading coverage needed for the agreed transient case.
3. Document the actual axes, rotation center, DOF ordering, units, wave-direction convention, complex phase convention, density, gravity, depth and source/input identity.
4. Keep Capytaine API/source investigations, limitations and reproductions in its expert knowledge base. Reuse existing reliable work without adopting old application-specific assumptions.
5. Coordinate frequency coverage and mesh studies with the WEC-Sim expert based on variation in time-domain motions, rather than blanket coefficient gates.
6. Provide concise export/readiness reports to the main thread. Do not require the main thread to debug Capytaine internals.

## WEC-Sim expert: consumer, BEMIO and time-domain ownership

The WEC-Sim expert OWNS BEMIO. BEMIO is the preprocessing/import layer for the WEC-Sim consumer. The Capytaine expert owns producing compatible source exports; both experts verify the boundary contract.

Own WEC-Sim, its maintained fork, BEMIO conversion, Simulink/Simscape case construction, startup and excitation handling, radiation memory, six-DOF integration, postprocessing and motion/visualization export.

Responsibilities:

1. Establish the working Capytaine/BEMIO/WEC-Sim pipeline using the existing toolchain before creating replacement solver machinery.
2. Verify required MATLAB toolboxes and licenses, pin dependency versions and record exact commits.
3. Own BEMIO parsing, normalization, HDF5 output, hydrostatics-file requirements and interoperability issues. Resolve them with the Capytaine expert rather than handing low-level debugging to the main thread.
4. Implement the agreed incoming-wave startup using supported capabilities where possible; report precisely where extensions are necessary.
5. Produce the 10-second motion histories and isometric video, then characterize motion variation with mesh resolution. Keep the displayed wavefield and actual excitation consistent.
6. Maintain a durable WEC-Sim/BEMIO code-path map, bug ledger, minimal reproductions, regression cases, patch rationale and upstream status.

## Topwater expert: implementation ownership

Own the selected Topwater fluid/body solver, its implementation architecture, numerical integration, implementation feasibility, tests and performance measurement. Maintain its deep implementation context, governing equations, supported scope, missing capabilities, known failures and reproducible evidence.

Responsibilities:

1. Own Topwater implementation and integration when that phase is authorized. Research and other experts may propose methods, but the Topwater expert is accountable for the resulting code and its behavior.
2. Review research proposals for compatibility with the existing solver, implementation cost, testability and likely computational cost before implementation.
3. Implement bounded experiments and corrections selected by the coordinator; keep one supported runtime path and change supported boats through data.
4. Consume the shared physical inputs and independently calculate fluid loads and body motion. Capytaine/WEC-Sim results are external reference data, not injected runtime forces or trajectories substituting for Topwater physics.
5. Validate resulting time-domain motions and report uncertainty and complete measured runtime. Keep component diagnostics distinct from evidence of end-to-end accuracy.
6. Own appropriate local tests, regression preservation and performance measurements. Major optimization campaigns remain held until the required physical comparison passes; reject clearly impractical routes early.

CURRENT BOUNDARY: Establish this role and preserve its knowledge, but do not begin Topwater solver changes. The user's immediate instruction remains to get the Capytaine/BEMIO/WEC-Sim external reference working first. Establishing the agent is not authorization to resume implementation or optimization.

## Numerical-methods research expert: ideas and literature ownership

Own the literature review and method-proposal pipeline. The Topwater expert owns code implementation; it is not expected to generate every underlying idea or independently repeat the literature survey.

Maintain a searchable, evidence-backed record covering adaptive time integration, error estimation, mesh generation/adaptation, boundary treatment, CFD/BEM formulations, fluid-body coupling and relevant acceleration methods. Prefer primary papers, official implementations and documented numerical evidence over recollection or general claims.

For each proposal, identify:

- The concrete project problem and why an alternative is needed.
- Primary sources with useful links and the specific result supporting the proposal.
- Governing assumptions, applicability to this boat/sea/solver and important failure modes.
- Competing approaches, expected benefit, implementation cost and principal risks.
- The cheapest bounded analytical or numerical experiment that could reject the idea.
- What is established in the literature versus still unproven in this project.

The normal flow is research proposal -> Topwater implementation-feasibility review -> coordinator selects a bounded experiment -> implementation owner executes -> independent assessment. This is an ownership workflow, not a mandatory ceremony for every routine bug fix. Research does not authorize a formulation change, mesh campaign or optimization experiment by itself.

Research may inspect code and the Topwater expert may originate ideas. Coordinate rather than creating rigid barriers. The researcher can assist reference experts with a specific literature question, such as transient incident-wave initialization. Do not launch a broad optimization study as a prerequisite to the current reference result, and do not treat publication results as validation of our implementation.

## Immediate WEC-Sim expert assignment: fork and independent correctness review

The WEC-Sim expert must perform the following work, not merely summarize the old findings.

### A. Establish a small maintained fork

Create a project-owned WEC-Sim fork based on a deliberately selected upstream version. Pin the tested commit. Keep generic engine corrections in the fork and case-specific geometry, forcing and visualization configuration in the project. Avoid runtime string replacement of generated models where a proper source/library correction is feasible.

Keep the fork close to upstream. Prepare fixes so they can be offered upstream. Do not turn this into a new solver or an unrestricted refactoring campaign.

### B. Independently reproduce the marine-gyro findings

Use these as investigation leads, not unquestioned diagnoses. Determine whether each issue still affects the selected upstream version. Compare stock behavior against independently derived expectations before carrying over any patch.

Project: `C:\Users\alex\marine-gyro`

Recorded engine commit in `data/case_study.json`: `0753b2e47f2457c078751dcfe5d251d1767b80ab`. Verify the actual configuration before relying on it.

1. Yaw rotation-order mismatch: `reports/wec-sim/run_case.m:25` substitutes `marineKinYaw` for `kinYawTransformNoB2B`. The recorded concern is a different Euler order from the quaternion sensor. Inspect `reports/wec-sim/marineKinYaw.m`; independently check combined yaw, roll and pitch, frame conventions and velocity/acceleration transformations.
2. Yaw postprocessing time grids: `run_case.m:67` recalculates wave elevation on `body1_out.time` before `stopWecSim`. Reproduce the unequal-grid behavior and decide whether the correct fix belongs upstream or in case configuration.
3. Directional excitation normalization: `docs/branch-review.md:127` records a missing angular-bin-width factor. The generated model uses `Nwave.amplitude .* Nwave.dTheta.'`. Verify the amplitude/energy normalization independently in the selected version.
4. Direction-array orientation: the same record reports errors from row versus column direction samples. Reproduce the shape/broadcasting behavior and establish an explicit input contract.

For each lead, record: affected versions, exact path, minimal reproduction, expected versus actual behavior, physical consequence, cause confidence, correction if needed and regression evidence. A workaround clearing a symptom does not by itself establish its cause.

### C. Review the code paths this project will exercise

Go beyond the known findings, but stay focused on the actual canonical workflow:

- Capytaine export -> BEMIO -> WEC-Sim mapping, including units, normalization, body origin/CG, DOF ordering, full inertia products, headings and phase conventions.
- Incoming wave train/packet, directional components, arrival timing, excitation impulse response, prehistory/padding and compatibility with a resting boat in locally calm water.
- Radiation-memory convolution, infinite-frequency added mass, frequency interpolation/band limits and initialization. Judge their practical effect through output motion sensitivity.
- Hydrostatic restoring and buoyancy ownership, avoiding duplicate terms; frame transforms and free six-DOF constraints.
- Solver/output time grids, interpolation, rotation representation and postprocessing.
- Visualization and motion export: confirm the depicted incoming waves are the ones driving the model, with correct timing and orientation.

Inspect documentation and implementation. Use focused analytical or numerical reproductions where they can expose a real defect. Do not interpret successful compilation or agreement between two calls to the same function as independent correctness evidence.

### D. Deliver findings and fixes

Return a concise report with confirmed issues, unconfirmed hypotheses and unaffected behavior clearly distinguished. Include the fork URL/commit, focused regression tests, proposed or applied minimal corrections, remaining limitations and the effect on the canonical workflow. Record deep details in the expert knowledge base rather than flooding the main thread.

Obtain independent review before relying on changes that affect physical results. Share interface findings with the Capytaine expert. Preserve failed runs and exact input/producer identity. Run numerical jobs serially and keep their inputs and producers immutable until the exact process is terminal.

## Shared interface and reporting rules

The WEC-Sim expert owns the BEMIO boundary investigation; both experts sign off on its physical conventions. Neither should make the main thread reconstruct the interface from raw code. Supply a small documented producer/consumer contract and a reproducible end-to-end example.

The main thread remains accountable for the user-visible outcome and scope decisions. It should ask experts for focused work, challenge unsupported conclusions and coordinate results, while leaving solver-specific depth with the experts.

## Setup and activation checklist for the main thread

This is an instruction to establish the working structure, not merely to acknowledge this document.

1. Inspect current project/client configuration and reconcile the latest user scope with existing instructions. Preserve current work and running-job state; do not resume retired investigations automatically.
2. Create the four supported project-scoped agent definitions, normally `.codex/agents/capytaine_expert.toml`, `wecsim_expert.toml`, `topwater_expert.toml` and `numerical_methods_researcher.toml`. Set their roles and delegation descriptions explicitly.
3. Create the expert knowledge entry points and shared physical/interface contract. Seed them with existing verified evidence, including the marine-gyro leads, while retaining their version and uncertainty. Do not blindly copy workarounds as facts.
4. Add concise project routing instructions: Capytaine questions go to its expert; WEC-Sim and BEMIO questions to theirs; Topwater implementation to its expert; literature/method selection to the researcher.
5. Verify the current client recognizes and can invoke the roles. If a restart, configuration reload or different supported mechanism is required, report the actual requirement. Do not imply the roles are active merely because files exist.
6. Give each expert its start-here material and require a concise scope/current-state confirmation. Keep Topwater implementation held. Do not spawn numerical jobs just to demonstrate that delegation works.
7. Assign the WEC-Sim expert the bounded fork/reproduction/exercised-path review above. Assign the Capytaine expert the producer/interface work needed for the first minimum-mesh case. Give the researcher focused questions only where they help an active decision.
8. Keep the immediate outcome visible: a reproducible 10-second incoming-wave reference, six motion histories, isometric video and motion variation with mesh resolution. Do not replace it with indefinite expert training, documentation, coefficient checks or whole-library auditing.
9. Obtain independent review of consequential corrections, preserve durable findings, and checkpoint the setup and actual results according to the project's branch/tracker discipline. Report which roles are operational, current task owners, repository locations and any genuine setup blocker.

## Sources for agent setup and upstream workflow

- Project-scoped custom agents: https://learn.chatgpt.com/docs/agent-configuration/subagents
- Reusable skills and supporting resources: https://learn.chatgpt.com/docs/build-skills
- WEC-Sim upstream development workflow: https://github.com/WEC-Sim/WEC-Sim#development

This document is a handoff instruction. Writing it does not itself create agents, create a fork or establish that any solver issue has been reproduced or fixed.


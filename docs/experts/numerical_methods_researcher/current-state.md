# Current state: numerical-methods research expert

Reference support is queued, not authorized; follow the
[handoff](../../checkpoint-handoff.md). Setup is complete. The downstream
project need is the Capytaine/BEMIO/WEC-Sim incoming-wave reference.
No broad formulation, optimization, adaptive-mesh or literature campaign is
authorized. No numerical jobs or solver implementation on role activation.

Primary starting references: Capytaine3.0.0 official documentation; selected
WEC-Sim7.1.0 source/docs and their theory references (radiation memory and
excitation IRF). Current source-backed issue: imported elevation uses noncausal
excitation convolution relative to a wave reference point, so spatial arrival,
initial state and actual timeline padding must be consistent. This is an active
interface question, not permission to invent a new solver.

Native activation confirmed 2026-10-24 as `/root/researcher_setup` using exact
`numerical_methods_researcher` registry role, with native developer instructions
and inherited parent settings. Setup is complete and stopped; no numerical work
or solver changes were performed and no execution lease exists.

Reusable [proposal template](proposal-template.md) and
[focused initialization/padding question](native-setup-2026-10-24.md) are ready.
The latter records exact read paths, selected official source commit
`a23d2993807410fe78cc63546c1ba9d81912024f`, source observations versus unproven
project assumptions, and a later independent-index rejection plan. Source
inspection is not a reproduction or physical validation.

Next action: await user-authorized reference work and a bounded coordinator
assignment for a focused question from reference owners on consistent wave
reference, excitation-IRF alignment and retained timeline padding. No broad
literature search or numerical experiment is assigned. Topwater implementation
and performance work remain held; durable evidence lives in these files.

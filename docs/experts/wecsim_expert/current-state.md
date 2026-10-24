# Current state: WEC-Sim / BEMIO expert

Updated2026-10-24. Scope/authorization: [handoff](../../checkpoint-handoff.md).
Owner: consumer/BEMIO/maintained fork/dynamics/startup/output/video. Coordinator
owns selection, acceptance and numerical allocation; no automatic follow-up.

REF-001 minimum reference is COMPLETE and independently audited. Frozen10s/501
native samples, actualrest, six signed applied-force histories and full mass/frame
checks passed. Final centered12m isometric video/plots passed literal encoded
frame review. [Run/evidence](minimum-reference-2026-10-24.md).
Maintained fork f962103a348d63101ce669c2611f7d7c109755c3 on user-owned
`topwater/reference-v7.1.0`, based on officialv7.1.0. Four reviewed corrections:
import history, component sum, actual Solid inertia order, selected rigid imported
force dispatch once. [Pin/limits](../wecsim-fork.json). No upstream/main publication.

REF-002 goal is ONE fresh buoy-generator case, roll5-10deg plus appreciable other
five motions, with an a-priori geometric linear-operating assumption. No bespoke
spectrum/phase/amplitude fit, canonical suite, new nonlinear adapter or custom
integrator. Topwater solver remains held. Generated whole banks/MEM spread/seed42
are immutable; internal heading selection does not itself select the canonical case.

Native17panel3-default-bank exact-pair H and .05:.05:12rad/s radiation completed.
Official BEMIO40sIRF/exact H5 metadata/full K/Ainf/all36kernel readback passed.
Native-vs40sK harmonic comparison completed1.859s, independently reviewed integrity
PASS. April globalroll7.551/7.802deg, max history difference.396deg roll/.05393m
heave; fixedmean relative-roll subtraction SCREEN proxy peaks5.929/6.069deg.
March globalroll10.201/10.228deg; January3.810deg native. These stationary
particular histories do NOT qualify initial pose, radiation history, startup/free
surge-sway-yaw drift, actual wet set or linearity. Flatwater5deg righting is already
12.3% below mean-K prediction; large globalroll has not rescued that issue.

Signed A33native-minus-finiteband ranges-3573.04kg to+375.79kg over incidentomega;
314/384negative. Raw sampled B33 is positive. Passive upper-tail-only explanation
is sign-inconsistent with dominant negative mismatch. Dominant error is native
versus truncatedbandPV, rather than40s/time quadrature alone; missing tails versus
coarse mesh/native consistency unresolved. No Ainf/damping repair or fit accepted.
[Executed evidence, arrays, normalizations, limits](generated-bank-consumer-execution-2026-10-24.md).

Jan05offset90 source/math review passed and ONE allocated stored screen completed
1.828s/exit0. Native/CIC roll5.119/5.221deg, exactfixedmeanplane relativeroll
2.105/2.133deg and pitch2.484/2.508deg; sampledfootprintdeparture.307m, initial.183m.
Provisional global-roll hypothesis retained, NOT canonical/linearity/IC acceptance.
Independent executed review requested. No body/ODE/BEM/MATLAB run follows.
[Executed screen, source/manifest and practical next step](jan90-screen-proposal-2026-10-24.md).
All owned numerical jobs terminal/leasefree; consult shared execution ledger
for CURRENT other-owner liveness, never infer it from this document.

Four old marine-gyro leads remain unconfirmed/unexercised as applicable; no old
workaround carried over. Historical setup/open questions live in linked records,
not current blockers. Before engine replacement, query authoritative outside-
sandbox IPC/exact recorded process identity; restricted visibility is insufficient.

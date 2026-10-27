# Simple-boat native spectral motion qualification

The corrected native WEC-Sim spectral path passed the fixed simple-boat
compatibility test. This establishes compatibility for this case, not general
physical accuracy or permission to promote the held dependency.

The case retains the original boat, sea, full-band operator, radiation memory,
inertia mapping, integration settings and initial calm/rest. WEC-Sim computes
excitation independently from physical spectral components. Reference forces and
motion are read only afterward for comparison; no time/phase/gain fitting occurs.

| Check | Measured result | Prospective limit |
|---|---|---|
| Full spectral force replay, all3001x6 samples | Worst relative maximum3.414e-15; relative RMS4.360e-15 | Each <=1e-10 of its original channel scale; zero channels exact |
| Actual native signal bus | All six channels pass | Same unchanged per-channel force gate |
| Actual outer body excitation | All six channels pass | Once-only replay, <=1e-10 of prepared peak |
| Maximum world-position vector difference | 0.05975798mm at58.80s | 1mm |
| Maximum relative orientation angle | 0.001569854deg at51.46s | 0.01deg |
| Initial body velocity/radiation output | Zero | Existing rest limits |
| First radiation force-history rows | Zero, times0,-.02,-.04s | Native zero-memory startup |

The fresh-process/no-import initialization certificate and frozen native zero
initializer support the negative velocity-prehistory inference. Persistent
negative velocity arrays were not measured directly. Full retained CIC replay,
force balance and original logging-only graph checks passed. Actual motion query
clock equals the prepared integer-multiply clock exactly.

Producer: WEC topic e24c1a73372baf37aebe37a8b6bf0cd0d2368b0c plus ONLY the
independently qualified PR9 inertia library7b0c063c8c4db2dced7e42257002989d69cdbb7a789a21f0a044ca26c793d828
in a detached private worktree. Neither topic/reference nor production pin moved.

Evidence under artifacts/workflow:

- native-spectral-full-force-review-v1/acceptance.json: f460dc37.
- native-spectral-actual-bus-review-v1/acceptance.json: 547becf2.
- native-spectral-motion-allocation1/report.json: a9cadd6e.
- native-spectral-motion-review-v1/acceptance.json: d934f305; binds raw motion,
  60-channel/CIC logs, clocks, zero-history certificate and original public result.
- native-spectral-motion-watch1/watchdog-receipt.json: 0f1d391a.
- native-spectral-motion-attempt1/terminal.json: exit0, frozen inputs unchanged.

Controlled motion allocation elapsed1382.281s, with kernel peak commit6.590GB.
Sampled aggregate working set peaked6.429GB with13 disclosed sampling gaps;
this is not a continuous RSS guarantee. Exact owned processes were terminal,
all7903 frozen inputs matched, and root released the numerical lease only after
independent host-CIM and idle-ledger readback.

Remaining delivery: reconstruct the total incident+diffraction+radiation surface
from this saved new trajectory, then export private standard result.npz/result.png
and render the complete video with normal logs. The existing public video proves
the visualization fix on the old saved result, not corrected-path end-to-end output.
General physical field/mesh/applicability gaps remain in issues31-34. Dependency
promotion, upstream submissions and merges remain held.

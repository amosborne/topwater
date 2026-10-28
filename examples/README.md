# Canonical examples

`simple-boat/` is the canonical 24 ft, 2,500 kg open-hull example. Its explicit
sea realization comes from observation `ndbc-46086-20231223T0040Z`, seed 42.

| File | Purpose |
|---|---|
| boat.npz | Physical hull and mass properties |
| sea.npz | Explicit wave components |
| result.npz | Computed motion, physical hull snapshot and total free surface |
| result.png | Full six-DOF histories with common translation/rotation row scales |
| video.mp4 | Saved motion and total surface, 800-square at 12.5 fps |
| sea.log, result.log, video.log | Actual command logs |

```powershell
uv run topwater generate-sea examples/simple-boat --record ndbc-46086-20231223T0040Z --seed 42
uv run topwater solve examples/simple-boat 100 --rest 10 --ramp 5 --solver=wecsim
uv run topwater render examples/simple-boat
```

Each directory argument is optional and defaults to cwd. Solve now automatically
publishes the NPZ and PNG together; render reads only result.npz. It preserves
the producer's water mesh and computes waterline intersections separately.

The positional time is the full-amplitude dwell. This case uses 10 s rest,
5 s ramp, 100 s dwell, 5 s ramp and 10 s rest: 130 s total. Carrier phases
refer directly to the dwell midpoint, 65 s. Packet construction and WEC-Sim
use a fixed 0.01 s timestep; water/video output uses 0.08 s spacing.

The saved result contains all 13,001 native motion samples and 1,626 validated
surface/video frames, including both endpoints. The encoded movie lasts
130.08 s because its final state occupies one frame interval. The original
boat and sea files are byte-for-byte unchanged.

The native calculation completed before a MATLAB cache-restoration failure.
The saved motion and completed total field were audited and continued through
the existing geometry replay, standard Result/PNG publication and ordinary
render command. The logs record that continuation explicitly. The final video
passed complete decoding with zero warnings. This is not a claim that the
original solve invocation completed successfully. The [run record](../docs/experts/packet-integration-run-2026-10-27.json)
and [warning triage](../docs/cli/warning-triage.md) retain the evidence and limits.

![Six motion histories](simple-boat/result.png)

This is a working **provisional linear reference**, not a physical-accuracy
qualification. Read the [command documentation](../docs/cli/solve.md) and
[actual warning triage](../docs/cli/warning-triage.md). Native coupled dynamics
are not yet implemented. WEC-Sim uses the maintained release-based reference
at `6b774228e1f6b31dee473e626204b0ef65e5ec38`. The earlier selected 60-second
surface-error qualification does not automatically qualify this changed packet.

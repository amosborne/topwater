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
uv run topwater solve examples/simple-boat 60 --solver=wecsim
uv run topwater render examples/simple-boat
```

Each directory argument is optional and defaults to cwd. Solve now automatically
publishes the NPZ and PNG together; render reads only result.npz. It preserves
the producer's water mesh and computes waterline intersections separately.

The recorded full run contains 3,001 motion samples and751 surface/video frames,
including both endpoints. The encoded movie is 60.08 s because the last state
occupies one frame interval. Public solve took 2393.74 s and render 226.93 s.
The PNG feature was added afterward and the image was generated from the same
completed result without rerunning physics; those original command logs remain
unchanged. Plot histories agree with the retained native histories to rounding.

![Six motion histories](simple-boat/result.png)

This is a working **provisional linear reference**, not a physical-accuracy
qualification. Read the [command documentation](../docs/cli/solve.md) and
[actual warning triage](../docs/cli/warning-triage.md). Native coupled dynamics
are not yet implemented. The current data use the old pinned WEC-Sim fork;
the new clean reference branch has not been promoted into Topwater.

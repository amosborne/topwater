# Examples

`simple-boat-native/` is the default Python solver example.
`simple-boat-wecsim/` preserves the MATLAB/WEC-Sim reference example.
They have identical `boat.npz` and `sea.npz` inputs and separate results, plots,
videos and command logs. Both use the standard 24 ft boat and 130 s wave case.

```sh
uv run topwater solve examples/simple-boat-native --solver native
uv run topwater render examples/simple-boat-native
uv run topwater solve examples/simple-boat-wecsim --solver wecsim
uv run topwater render examples/simple-boat-wecsim
```

The reference solve requires MATLAB/Simulink and the WEC-Sim fork; the native
solve does not. See [the workflow guide](../docs/workflow.md) for input commands.

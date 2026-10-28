# MATLAB workflow tests

These tests cover WEC-Sim coefficient exchange, spectral excitation, hydrodynamic
input guards and MATLAB engine setting restoration. They exercise implementation
behavior; passing them does not establish physical accuracy.

Use MATLAB R2025b and the pinned WEC-Sim checkout described in
[setup](../../docs/setup.md). Run the suites explicitly when relevant:

```matlab
results = runtests('tests/matlab');
assertSuccess(results);
```

The BEMIO handoff tests require `TOPWATER_WECSIM_PATH` to identify the maintained
checkout. The fixtures use temporary files and synthetic coefficients rather
than running a boat simulation.

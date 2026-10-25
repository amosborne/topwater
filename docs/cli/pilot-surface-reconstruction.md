# Retained-field surface pilot

Developer diagnostic: `scripts/pilot_surface_reconstruction.py`. This is not
an alternative public solve command and does not write canonical example files.
Run only under the project's serial execution claim with frozen inputs.

```powershell
uv run --group wecsim python scripts/pilot_surface_reconstruction.py --physical-case reference/capytaine/canonical-boat.json --native-mat PATH/native.mat --output NEW_DIRECTORY
```

## Inputs and outputs

Inputs are the archived canonical physical JSON and accepted exact native MAT
containing M/K/CG/CB/volume/water constants. The pilot deliberately requires the
unchanged17-panel mean wetted hull. Output is a fresh diagnostic directory with
native coefficient exports, pilot.npz fields, and a JSON execution report.
It never integrates a boat or changes a stored trajectory.

## Assumptions and hard-coded choices

Capytaine3.0.0, infinite depth, zero forward speed and one numerical thread.
Thirteen sparse radiation frequencies span0.025..96rad/s, six modes, three
exterior XY points, two synthetic incident frequencies and two directional
carriers with both signed spatial branches. Expected137solves includes infinity
and independent directional-linearity controls. The checked wall budget is180s;
a single native call can return after that budget before the next check occurs.

## Fatal errors and warnings

Before solving: wrong physical mesh, inconsistent native mass/CG, unsupported
Capytaine version or existing output directory are fatal. During the pilot:
failed solve, nonfinite fields, failed linearity assertions or exceeded checked
budget stop it and record failure. Prior source artifacts remain unchanged.

Capytaine mesh-resolution and irregular-frequency warnings remain diagnostics.
The intentionally sparse, high-frequency probes do not establish accurate
high-frequency hydrodynamics. This command never certifies total-surface
readiness, causal-kernel accuracy, moving contact or body-motion convergence.

## Detection limits

The three points are away from the mean waterline. Directional linearity checks
compare combined versus separate solves of the same synthetic field; they do
not establish actual packet/old force-table agreement, total initial calm or
near-hull continuation accuracy. Those are separate acceptance steps.

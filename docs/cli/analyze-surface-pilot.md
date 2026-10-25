# Surface-transfer analysis and refinement

Developer diagnostic: `scripts/analyze_surface_pilot.py`. Run under the serial
execution procedure; it does not integrate motion or publish result.npz.

```powershell
uv run --group wecsim python scripts/analyze_surface_pilot.py --input PATH/pilot.npz --output NEW_DIRECTORY
```

Add `--refine --physical-case reference/capytaine/canonical-boat.json
--native-mat PATH/native.mat` to perform the bounded direct-field refinement.

## Inputs and outputs

The input is the sparse pilot's numeric field archive. Output is a fresh
diagnostic directory containing a JSON analysis report, candidate kernel data,
and, in refinement mode, retained-field/native products and refined-transfer.npz.
No example input, stored motion, or applied body load is changed.

## Assumptions and hard-coded choices

The exp(-i*omega*t) convention gives velocity transfer `i*eta_q/omega`.
The candidate is a finite-band causal cosine inversion, evaluated over60s
at0.02s cadence. Cutoffs12/24rad/s and spacings0.1/0.05rad/s are compared.
The diagnostic low-frequency endpoint is zero and the direct term is zero;
these are measured hypotheses, not automatic physical acceptance.

Default mode interpolates sparse stored data and cannot qualify the quadrature.
Refinement directly solves approximately3000same-mesh problems at the three
existing exterior points, adds independent held-out frequencies, and uses a
300s between-calls budget. Scratch Green-function caching avoids repeated
precomputation within the job. Finite-step endpoint/tail accuracy is assessed
through full complex recovery, not merely the reconstructed real part.

## Fatal errors, warnings and detection limits

Missing/malformed data, reused output directory, missing refinement inputs,
failed solves, nonfinite transfers or exceeded checked budget are fatal.
Mesh-resolution/irregular-frequency warnings do not become accuracy proofs.
Reports retain `total_surface_ready=false`: neither a fitted endpoint diagnostic
nor a small three-point error establishes full-patch, near-waterline, actual
packet, initial-total-calm or nonlinear moving-contact accuracy. No parameter
is selected to force a small total-minus-incident surface difference.

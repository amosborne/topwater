# Actual radiation work: independent stored-data review

Coordinator record of the numerical-methods reviewer's 2026-10-24 read-only
audit. No new solver run, model correction or coefficient repair. Latest user
requests controlled changes to the passing official WEC-Sim example.

## Definition and scope

Use all six native world linear/angular velocity channels and the corresponding
resisting radiation-memory force/moment channels. Native applied radiation force
has the opposite sign. At every saved major time sample,

`P_R = sum(v_world .* Frad_resisting, 2)` and `W_R = cumtrapz(time, P_R)`.

These are full six-channel sums, not selected diagonal/component powers.
Positive work means transfer out of the body's motion into the radiation-memory
subsystem. Integrate positive/negative instantaneous parts separately using the
same 0.02s saved grid. Negative instantaneous power can represent stored-energy
return. This is not a complete fluid-energy budget or a general passivity proof.

## Results

| Zero-memory case | Final W_R | Minimum cumulative W_R | Integrated negative-power magnitude |
| --- | ---: | ---: | ---: |
| Original incoming field extended to60s | +7.390385095J | 0J | 0.193731967J |
| Earlier Jan90 established-sea held release | +237893.6506J | 0J | 5236.6721J |

Original60s: 601/3001 saved power samples negative, minimum power -0.368689W.
Its10s prefix gives +1.253275844960606J, independently matching the original10s
run's +1.253275844960607J. At first5deg roll in the Jan90 case (34.64s), cumulative
CIC work is already +18042.95J. Neither case has negative cumulative radiation
work at any saved sample. These observations do not exclude mode-specific
underdamping or every coefficient error, but do not support attributing the
large response to demonstrated net CIC energy generation from startup.

The Jan90 solver completed; its approximately107deg roll/11.5m sway response was
rejected as outside the intended linear validation regime. That is not an
established numerical failure. It started from rest in an established incident
sea, not the requested locally calm incoming-wave initial-value problem.

## Evidence and limitations

- Original60s: `artifacts/reference/reference-extension-20261024/case-attempt1`;
  motion SHA256 `d052bdf0831bed728f255a4a928e9f9411b804cc4e2f5abc0d1611f5abf5ce8a`;
  raw native output `722759da50e958b748e386e6a5024443c6a5f2ea5b742b3a12993f66c1cef6f0`.
- Original10s: `artifacts/reference/minimum-20261024/case-attempt1`.
- Jan90: `artifacts/reference/generated-bank-consumer-20261024/case-jan90-attempt1`;
  immutable raw trajectory and original rejected-case controls remain preserved.
- Original60s independent full3001 state/force/CIC/wrench review PASS: all142
  frozen inputs and recorded output hashes match, initial v/radiation zero,
  empirical viscous channels zero, applied Fe error<=2.27e-12, rotated Solid
  Newton-Euler<=3.41e-13, native Ar recurrence<=1.45e-12, actual RT/three-row
  PCHIP consumption from independent zero-past CIC<=1.71e-13. First501 q/v match
  original trajectory within2.78e-17/7.81e-18. Exact job terminal/released.

The separate prescribed-motion native pressure-work witness still exposes a
negative harmonic-power direction in the17-panel coefficients. Its practical
consequence for these trajectories has not been established. Changing17 to68
panels showed sensitivity to spatial discretization; panel density as root cause
is not established. No WEC-Sim/Capytaine radiation bug is established.

Official examples are wave-driven with PTO/constraints. Their net radiation loss
is expected and does not demonstrate an energy-conservation defect. A complete
wave/body/PTO/fluid-storage energy balance has not been closed. Added mass in the
actual Solid and the native residual-Ar feedback must not be counted twice.

Next: controlled modifications of the passing official regularCIC case, with
explicit work accounting and one feature changed at a time. No further custom
sea, mesh solve, damping fit or coefficient correction follows from this review.

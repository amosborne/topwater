# Native spectral force route: read-only review, 2026-10-26

WEC-Sim already computes linear excitation by summing harmonic wave components.
For our known packet, an excitation IRF is not mathematically necessary. Prefer
one vectorized native spectral input: Topwater supplies physical components and
BEM coefficients; WEC-Sim calculates forces and motion. This is a recommendation
for investigation, not an implemented or tested replacement.

## What already works, and what does not

| Route | Selected-source finding |
| --- | --- |
| `spectrumImport` | A three-column file supplies frequency Hz, spectrum m2/Hz and explicit phase radians. Imported frequencies are retained, including nonuniform spacing; three columns avoid random phases. One exact direction is directly representable. |
| Ordinary `irregular` | Traditional/EqualEnergy grids and generated phases change our exact components. |
| One ordinary direction array | A common spectrum times fixed directional spread, with imported phase only one column, cannot represent arbitrary complex frequency-by-direction components. |
| `spectrumImportFullDir` | Supports frequency-dependent directional amplitudes, but pinned v7.1 always randomizes its phase matrix. Additional normalization/interpolation issues below prevent claiming a drop-in replacement. |
| Regular groups | No public arbitrary phase input; regular elevation has zero phase. Plain `regular` also changes radiation; `regularCIC` is needed to retain memory. |

## Exact physical mapping

Our packet coefficients `C[j,sign,m]` already contain amplitude, phase, spatial
window, Jacobian and frequency quadrature weight. The reverse spatial branch
uses positive frequency at `theta+pi`, without conjugating C. For one direction,
`C=a*exp(i*psi)`, native imported spacing `dOmega`, and `spread=1`, use:

```
frequency_Hz = omega/(2*pi)
spectrum_Hz = pi*abs(C).^2./dOmega
phase_rad = -angle(C)
```

Native amplitude is `sqrt(2*(spectrum_Hz/(2*pi))*dOmega)=abs(C)`. Do not apply
another packet weight or factor two. Preserve phase pi from negative window
coefficients. Use the actual imported spacing, including its endpoint rule.

With our verified H5 convention `RE=real(H)/(rho*g)`, `IM=-imag(H)/(rho*g)`,
`irregExcF` computes `a*(real(H)*cos(omega*t-psi)+imag(H)*sin(omega*t-psi))`,
which equals the old `real(C*H*exp(-i*omega*t))`. Native spatial elevation uses
`cos(omega*t-k*d.x+phase)`, consistent with the same origin. Retain absolute time
zero, exact physical directions, fixed mean-CG basis and disabled displacement
phase. Only identical frequency AND direction components may be combined by
complex addition; PSD averaging would destroy the calm-start interference.
The25 BEM headings are interpolation nodes, not replacement physical directions.

## Why it is not a configuration-only switch

- **Ordinary groups share first-wave coefficients.** `hydroForcePre` selects
  `waveObj(1)` for spectral H preparation; ordinary cloned branches receive the
  same `fExt` bus. Different directions cannot be assumed correct. The active
  Rigid mask creates one variant branch and Add port per group: G groups means
  G branches, potentially736 for the existing records. This is a source count,
  not a compile/runtime measurement. Imported once-dispatch does not fix it.
  Pinned `docs/user/advanced_features.rst:333-356` explicitly promises unique
  group directions and gives a0/90-degree example: this is an apparent defect
  against the documented contract, requiring a minimal native-force reproducer.
- **Single FullDir needs focused repairs or a small alternative input.** Its
  pinned phase reader ignores a supplied phase matrix. Its elevation amplitude
  includes `sqrt(dTheta)`, while its native force amplitude omits that factor;
  active Rigid wiring passes `Nwave.amplitude` directly. This source-level
  inconsistency needs a tiny numerical witness. FullDir also prepares H with
  spline interpolation, unlike the old piecewise-linear real/imaginary operator.
  Pinned `docs/user/code_structure.rst:353-355` already promises an optional
  frequency-by-direction phase matrix; the mismatch is not merely newer docs.
  Ignored phase and direction-width inconsistency are apparent defects pending
  minimal execution, while choosing spline interpolation is a method choice.
- **Low positive frequencies must survive.** Raw H5 begins at0.1rad/s; imported
  setup filters outside its range. Our packet has positive sidebands below0.1.
  Preserve their existing `H0+(omega/.1)^2*(H(.1,theta)-H0)`, `H0=K(:,Heave)`.
  Only C(0) is exactly zero and can be omitted. A native handoff/extension must
  carry this bridge without inventing low-frequency radiation A/B or changing
  the accepted radiation kernel/Ainf. Frequency/phase roundtrips need checks.

The smallest promising architecture is one vectorized native spectral payload,
through a repaired FullDir input or a small generic discrete-component input to
WEC-Sim's harmonic evaluator. Choosing between them remains a proposal. Neither
requires Topwater to calculate applied forces or hundreds of model clones.

## Cheapest falsifier and stopping point

Propose one native force-branch example: two headings, three nonuniform positive
frequencies, unequal amplitudes, explicit phases, distinct complex six-channel H,
zero mean drift; include an off-grid H query and one sub0.1 bridge sample. Check
origin/nonzero-position elevation, analytic forces and actual downstream bus.
Two unequal direction widths expose the FullDir normalization issue cheaply.
Do not run body dynamics first.

If that passes, verify the unchanged full packet at all3001 times against the
frozen spectral baseline, then motion under the shared contract. Keep ramp0,
initial rest and zero radiation memory. Native irregular/regularCIC retain CIC
radiation, so harmonic excitation does not require steady-state body motion or
warmup. No imports, numerical jobs, implementation, fork/pin/issue/PR changes
were made. Halted FFT candidate2f039e5 and failed replay remain preserved.

## Exact source references

Inspected held9952634f6f2e2746a1626484e27c3003d1e1d5a9 and combined reference
25b58bda96ab0272cbfb4b0755539d5bb8ca227e. Relevant waveClass and harmonic
functions are identical between them; bodyClass differences concern imported
elevation only. Sources below are pinned to held995:

- [waveClass](https://github.com/amosborne/WEC-Sim/blob/9952634f6f2e2746a1626484e27c3003d1e1d5a9/source/objects/waveClass.m#L333):333-364 imported grids/widths;596-623 phases;741-792 density conversion;820-830 ordinary field;859-867 FullDir field. SHA2563a1bf7378dd1e2e9bc9aa2ef27b4476a1a0ab8a352541b63387304183205fb09.
- [bodyClass](https://github.com/amosborne/WEC-Sim/blob/9952634f6f2e2746a1626484e27c3003d1e1d5a9/source/objects/bodyClass.m#L440):440-453 first-wave selection;1030-1055 ordinary H;1088-1093 FullDir spline. SHA256697552816c164c34207a7cdb32002d43a139871203864571764ef9b08f0a59bd.
- [irregExcF](https://github.com/amosborne/WEC-Sim/blob/9952634f6f2e2746a1626484e27c3003d1e1d5a9/source/functions/simulink/model/irregExcF.m#L5):5-19 harmonic sum; [irregExcFullDirF](https://github.com/amosborne/WEC-Sim/blob/9952634f6f2e2746a1626484e27c3003d1e1d5a9/source/functions/simulink/model/irregExcFullDirF.m#L87):87-98 amplitude/sum.
- Library `source/lib/WEC-Sim/WECSim_Lib_Body_Elements.slx`: archive `system_8_4071.xml` mask replication; `system_8_4071_7114.xml` FullDir amplitude wiring. Read as XML only, no model load.
- `initializeWecSim.m:250-262` frequency range;404-415 radiation variants. Local `src/topwater/_wecsim_setup.py:131-194` packet/offline operator; `data/wecsim/prepare_workflow_hydrodynamics.m:118-125` H5 sign normalization.
- [Pinned FullDir documentation](https://github.com/amosborne/WEC-Sim/blob/25b58bda96ab0272cbfb4b0755539d5bb8ca227e/docs/user/code_structure.rst#L353) already promises optional phases, and [pinned multiple-wave documentation](https://github.com/amosborne/WEC-Sim/blob/25b58bda96ab0272cbfb4b0755539d5bb8ca227e/docs/user/advanced_features.rst#L333) promises unique directions. These contracts strengthen the apparent-defect findings; they are not numerical reproductions.

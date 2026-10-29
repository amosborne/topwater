# Native solver and external models



`topwater solve` defaults to the Python solver in `src/topwater/solve/native/`.

Surface reconstruction accepts `--device auto|cpu|gpu|"GPU name"` (default auto), also

available as `simulate(..., device='auto')`. GPU installation

and compatibility are described in [GPU reconstruction](gpu.md).

Capytaine supplies the existing mean-hull coefficients. No MATLAB imports,

engine discovery, WEC-Sim checkout, or fallback backend is used by native solve.

One coordinator in `solve/workflow.py` owns case preparation, resource lifetime,

surface reconstruction and Result assembly for both native and WEC-Sim motion.

`solve/native/` and `solve/wecsim/` implement motion operations; the independent

`solve/cpu/` and `solve/gpu/` directories implement the same field-operation

interface. Shared calculation code lives directly under `solve/`.



CPU and GPU surface reconstruction run through one coordinator in

`solve/fields.py`. Source caching, batching, allocation backoff, Fourier

preparation, and output assembly are shared. Device-specific code implements

the numerical kernels and memory operations. Capytaine supplies body BEM source

solutions and tables; Topwater's specialized potential-only kernels evaluate

the surface fields. Both paths retain the same documented numerical policy.



## Dynamics and reference



The solver supports one freely floating six-DOF rigid body with coherent

directional excitation, linear restoring, and full-history radiation convolution.

The standard timestep is 0.01 s. Fixed-step RK4 integrates CG position, quaternion

attitude, translational/angular velocity, and optional external continuous states.

Angular velocity is a physical vector, not an Euler-angle derivative.



The instantaneous acceleration matrix is derived from the reference's compiled

Newton-Euler balance: the rotated augmented solid inertia and gyroscopic moment

are retained, and the residual added-mass acceleration term moves to the left

side. At mean attitude the matrix is exactly `M + Ainf`. This removes the

reference's nominal 1e-7 s delay and its acceleration extrapolation. Genuine BEM

asymmetry stays in the residual, rather than being silently symmetrized. This

preserves the reference's continuous-time equations; it is not a new large-angle

added-mass derivation.



Radiation-kernel preparation analytically integrates the piecewise-linear

`B/omega` interpolant against `(2/pi) * omega * cos(omega*t)`. Centered spherical

Bessel moments keep zero and short lags stable. This preserves the reference's

normalization and interpolation order without a frequency quadrature grid. Convolution uses causal

trapezoidal quadrature, zero prehistory, provisional stage velocities, and

stage-local endpoints over the full available history. Internal RK stages

never commit history. Native excitation solves directly at the finite packet frequencies and actual

directions. Incident boundary conditions and Froude-Krylov pressure are combined

coherently before each BEM solve, yielding one diffraction solution per active

packet frequency. Motion and reconstruction share those sources; there is no

native frequency/heading interpolation table. Native preparation stays in memory

and creates no WEC-Sim NetCDF/MAT or hydrostatics exchange files. The optional

reference backend retains its original table and export path.



`tests/data/simple_boat_wecsim.npz` freezes the current 130 s simple-boat reference,

including input hashes, fork revision, force histories, and consumed operators.

It can be recreated with `scripts/capture_native_reference.py` using the optional

WEC-Sim environment and shared engine; the script refuses to overwrite a reference.

Ordinary tests use only numeric arrays and need no MATLAB.



The measured WEC-Sim comparison is saved in

`tests/data/simple_boat_native_comparison.json`. It records historical numerical

agreement on this case. There is no production translation/rotation acceptance

check or angular-accuracy input. Linear-model validity is left to users inspecting

the results; solve does not assess departures from it.



## Calculation options



```python

from topwater.solve import SolveOptions, simulate



options = SolveOptions(time_step_s=.005, render_every=20, packet_tolerance_m=.005,

                       surface_padding_beams=1.)

result = simulate(boat, sea, 100.02, options=options)

```



`SolveOptions` in `solve/settings.py` owns calculation defaults. The

[workflow option table](workflow.md) lists the everyday CLI controls. Render stride,

view padding, spectral, radiation and resource overrides remain Python API controls.

The 1 cm packet tolerance determines Fourier packet construction; it is not a

motion error tolerance. No 0.1° comparison-only input is exposed.



The default motion step is 10 ms and one water/video frame is saved every 10 steps

(`render_every`), producing 10 fps. Frame rate follows the selected timestep and

stride. Surface samples

are native motion samples, plus the exact endpoint when it is off cadence. A

movie skips that short final interval. The solver determines a world view enclosing

the full hull trajectory with configurable padding, constructs a water grid

covering it (ten video pixels between samples), and saves the bounds and interval in the result. Rendering consumes

those values without selecting the reconstruction domain or resampling water.



## Private package integration



Your private project declares Topwater as a dependency and installs both packages

in one Python environment. Its driver imports its own model and passes it to

Topwater. Topwater never imports a private package by name or discovers plugins.



```python

from topwater.solve import simulate

from topwater.workflow_data import Boat, Sea

from private_gyro import GyroSystem



result = simulate(

    Boat.load('case/boat.npz'), Sea.load('case/sea.npz'), 100.,

    rest_s=10., ramp_s=5., influence=GyroSystem(...),

)

result.save('case/result.npz')

```



The public types are `BoatState`, `LoadEvaluation`, and the structural `Influence`

protocol. Inheritance and registration are unnecessary. A model provides:



```python

import numpy as np

from topwater.solve import LoadEvaluation



class OutsideModel:

    def initial_state(self, boat):

        return np.zeros(0)  # optional continuous device states



    def on_step(self, boat, state):

        pass  # optional: update sampled controller memory and held commands



    def evaluate(self, boat, state):

        return LoadEvaluation(

            force_world_N=np.zeros(3),

            torque_world_Nm=np.zeros(3),

            state_derivative=np.zeros_like(state),

        )

```



The methods receive read-only copies. `BoatState` contains `time_s`, `position_m`

at CG, `orientation_wxyz` (active body-to-world quaternion), `velocity_m_s` at CG,

and `angular_velocity_rad_s`. All vector quantities are resolved in world

east/north/up coordinates. Body axes are bow/port/up. Use `topwater.geometry.transforms`

helpers to rotate loads or shift moments to CG.



`evaluate` receives the provisional boat and continuous device states, and

returns world force [N], world torque about the current CG [N m], and the numeric

device-state derivative. It must be pure: no advancing controller memory,

recording an accepted sample, or performing irreversible actions on RK calls.

Continuous states participate in the same RK4 stages as the boat.



`on_step`, when supplied, runs at initialization and after every completed step,

including the last. It may change the model's held command or discrete memory;

it must not advance the continuous state independently. Its return is ignored.

The initial control cadence is 100 Hz; slower controllers can update at integer

multiples of the selected timestep. Custom timesteps are available through

`SolveOptions.time_step_s` and in the lower-level

integrator for convergence work; arbitrary control-event scheduling and separate

process transport are not implemented. The caller should create or reset its

model for each simulation.



External failures include model identity and simulation time. Nonfinite or

incorrectly shaped loads/derivatives fail the solve; no result is substituted.

The regular CLI has no external-model option. Private applications call Python.



## Ownership and boundaries



No damping coefficients, gyro plant, controller, or actuator model ships here.

A future model can implement roll damping or a closed loop of boat measurements,

controller commands, gyro dynamics, and hull reaction torque through this boundary.

Water-relative drag can query a caller-owned water-velocity model from `evaluate`;

the current interface does not supply local fluid velocity.



Boat mass properties include carried hardware exactly once. Influences add loads

and continuous states, not implicit mass. A reduced gyro model may return hull

reaction torque, but a full device formulation with unknown-acceleration coupling

needs a future mass-coupling interface or an explicitly derived reduction. The

present interface does not promise arbitrary multibody dynamics or DAEs.



The lower-level `topwater.solve.native.dynamics.integrate_motion` returns motion,

continuous-state histories, and named load/balance components for diagnostics.

`simulate` returns the existing self-contained `Result` archive format; it does

not add proprietary device fields. Private models can record accepted states in

`on_step` using their own storage.



## Choosing numerical discretization



The radiation defaults retain the established policy; excitation now uses the

packet frequencies and directions directly.

They are Python API controls, not CLI options. Their replacements should follow

the computed response and sampling requirements rather than boat length alone.

Exact excitation, analytical kernel integration and direct FFT-bin water

transfers are implemented. Motion coefficient spacing and coverage remain

chosen numerical policies.



| Current choice | Better basis |

|---|---|

| Radiation memory | Native uses the full case history; there is no truncation setting |

| Kernel quadrature grid | Removed: exact analytical integration of the coefficient interpolant |

| 0.0159155 Hz motion coefficient spacing | Resolve coefficient curvature and the retained kernel time horizon; refine where interpolation changes the response |

| Excitation frequency/heading grid | Removed from native solve; use the active packet frequencies and directions directly |

| Minimum 1.90986 Hz motion radiation coverage | Resolve the radiation kernel and its high-frequency tail/asymptotes; incident bandwidth alone is insufficient |

| Water radiation frequencies | Every native velocity FFT bin, spacing `1/(N*dt)` Hz through the highest bin; no separate cutoff or interpolation grid |



The solver's internal output catalogue owns `VIDEO_SIZE_PX = (800, 800)` and

`SURFACE_PIXELS_PER_SAMPLE = 10`. The physical water-grid spacing follows the

world view's metres per pixel. `Result.video_size_px` records the intended video

dimensions, and rendering consumes them. Archives predating this member retain

800 by 800 output compatibility.



The actual native FFT bins are separated by `1/(N*dt)` Hz, where `N` includes

the saved motion endpoint. Hydrodynamic frequency resolution also affects how

far the cosine-transformed kernel can be resolved in time. The coefficient interpolation error needs a numerical error criterion;

scaling by `sqrt(L/g)` gives a characteristic time but does not choose a justified

cutoff multiplier. The packet-envelope budget is a separate construction target.



The radiation convolution and cosine-transform definitions are documented in

[WEC-Sim theory](https://wec-sim.github.io/WEC-Sim/dev/theory/theory.html#convolution-integral-formulation).



Use the canonical `topwater.solve` public API for external controller integrations.


Reference MATLAB assets and exchange tables belong to `solve/wecsim/`;
static flotation equilibrium belongs to `boat/hydrostatics.py`.

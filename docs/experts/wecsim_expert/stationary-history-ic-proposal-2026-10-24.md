# Same Jan90 sea: proposed harmonic state/history initialization

Current status: the reduced case-only helpers were authored and ONE native
CIC/PCHIP/joint fixture passed independent executed review. Actual60s initialized
case sources are prepared for independent review ONLY; no body execution grant.
See [fixture receipt](seed-initialization-fixture-proposal-2026-10-24.md).

PROPOSAL ONLY. No implementation, engine call or numerical launch. The failed
held-release case, successful calm control, generated bank and17panel operator
remain immutable. Initializing a state does not repair/prove passive radiation
or guarantee a credible trajectory. Source/math/API review precedes any authoring.

Choose ONE concrete implementation direction: case-level missing-prehistory
correction BEFORE native radiation extrapolation, plus native joint state targets.
Coordinator narrowed the FIRST implementation proposal to leave acceleration
delay/default initial output/extrapolation UNCHANGED. Do not seed undocumented private
persistent arrays, replace CIC, impose motion, change coefficients, add damping,
change phases or make a custom time-domain integrator.

## Deterministic seed and its physical limitation

Use the existing Jan90 finite40sK component particular amplitudes Q_kernel40,
same384H/a/omega/theta/phi and fullM/K as the candidate seed. Freeze those arrays.
For s<=0 evaluate q_h(s)=Re(sum Q_j exp(-i*omega_j*s)), its first/second
derivatives analytically. No phase selection, fitting, window shift or modal
removal. Propose q0=q_h(0), nativeCG velocity and world angular velocity from
that history, instead of zero velocities which created large free-mode offsets.

This is a MEAN-LINEAR harmonic initialization approximation, not an exact
stationary solution of native rotating-Solid/gyro/discrete-CIC/delay equations.
The old harmonic angular coordinates approximate rotations in the linear solve.
Finite orientation must be handled consistently when assigning native kinematics.
Report pre0 and initial native equation residuals; do not call the state a
free-floating native equilibrium or initialized nonlinear fluid field. No
prescribed trajectory/forcing is applied for t>=0: all six native DOFs evolve.
If this approximation is materially inconsistent, reject before body simulation
rather than conceal its homogeneous transient or impose q_h as motion.

The sea's incident/diffraction field remains the existing mean-body linear
frequency-domain assumption. Actual initial immersion/curved-surface geometry and
mesh uncertainty still need their separate audits. A future nonlinear state
would require fluid potential/velocity/scattered-field and body histories on the
actual geometry, not merely these linear body/memory variables; none is proposed
as a new nonlinear model here.

## Native pose and CG/pivot velocity map

Body sensor Euler convention is R=Rz(psi)Ry(theta)Rx(phi). Set initial body pose
through body.initial.displacement and axis-angle for R_h(0), exactly once;
meanCG/STL/inertia are not reoriented by the external global heading again.

For Euler rates e_dot=[phi_dot,theta_dot,psi_dot], use WORLD omega=E(e)e_dot,

```
E = [cos(psi)cos(theta), -sin(psi), 0;
     sin(psi)cos(theta),  cos(psi), 0;
     -sin(theta),        0,        1].
```

WORLD alpha=E*e_ddot+E_dot*e_dot; translational v/a are q_h translation
derivatives. This same world velocity/acceleration convention must be used in
prehistory and native readbacks, not angular Euler rates substituted for omega.

Actual native constraint source is WECSim_Lib_Constraints.slx/system_3.xml,
CONSTRAINT (Bushing Joint). It exposes Px/Py/Pz and Rx/Ry/Rz PositionTargetSpecify,
VelocityTargetSpecify, Priority, Value and ValueUnits through proper set_param.
Keep motion ComputedMotion, force/torque actuation NoTorque, zero stiffness/
damping/limits unchanged. Set initial joint positions zero/high priority and
initial velocity targets only (not motion actuation).

Because meanCG and initial displacement are DOWNSTREAM of that joint, its pivot
is not the body's CG. At the proposed zero joint orientation, let
c=meanCG+initial.displacement (relative to constraint.location=0). Set rotational
primitive rates to WORLD omega0 (joint axes coincide with world at joint0), and
prismatic rates to v_CG0-omega0 cross c. The body axis-angle offset R_h(0) is
downstream, so assigning body Euler rates directly to primitive rates is wrong.
Acceptance must read actual CG v/worldomega and actual pose at0 against this
independent construction; no first-sample rebasing. Joint target API frame/order
must be independently traced before reliance.

## Exact finite CIC prehistory without private-state mutation

Selected fork f962103a: convolutionIntegralInterp.m creates velocityHistory=0,
shifts once per main CIC call, inserts CURRENT velocity in row1, multiplies
irkb and trapz(cicTime). It retains three force/time samples; extrapFrad.m uses
PCHIP extrapolation at minor solver stages. PCHIP is not additive in separate
histories, so adding a prehistory force AFTER extrapFrad is NOT equivalent.

Let dt=.02,T=40,N=2000,K_j=K(j*dt), and trapezoid weights
w_0=w_N=dt/2, w_j=dt otherwise. At main t_n=n*dt, native zero-history CIC
already contains actual post0 velocities for lag j<=n. Missing correction is

```
C(n) = sum_{j=n+1..N} w_j*K_j*v_h((n-j)*dt), 0<=n<N;
C(n) = 0, n>=N;
C(n) = sum_{j=0..N} w_j*K_j*v_h((n-j)*dt), n<0.
```

v_h includes worldomega, not Euler rates. At n=0 the native term already includes
HALF dt*K0*v0: C(0) starts at j=1, avoiding double counting. At negative times,
native initial FradHistory is zero, so C supplies the complete past force.
Precompute v_h down to -T-3dt, not just -T: native paired RateTransitions also
initially expose[-dt,-2dt,-3dt]. Their force InitialCondition is the fullpast
forces at exactly those times, then the combined updated three rows take over.

Proper case-only Simulink API inserts an Add on the native FradExtrap(3x6)
signal BEFORE its input to extrapFrad. For native timeExtrap rows[n,n-1,n-2],
add [C(n);C(n-1);C(n-2)]. Use the actual timeExtrap samples as lookup indices,
not an independently interpolated clock/table. Thus each native PCHIP call sees
the combined physical past+actual post0 force history; its sample timing and
minor-stage extrapolation remain native. Missing term becomes identicallyzero
after the last two history rows pass T. Excitation bus remains ONLY exact H;
the seeded contribution belongs to radiation, with its native resisting sign.
All post0 velocity samples remain measured native output, not q_h forcing.

No generic fork correction is claimed. This is an explicit initial-history case
adapter, with archived proper-API model changes and unchanged original library.

## Native residual acceleration delay: unchanged numerical startup

Ar remains the existing world-resolved residual, Solid I rotates natively.
TransportDelay SID8:4071:6251 uses delta=1e-7s/default initial0. Executed main
samples show a_D(0)=0, a_D(dt)=a_actual(0), then
a_D(n*dt)=(1+r)*a_actual((n-1)dt)-r*a_actual((n-2)dt), r=1-delta/dt.
This computational zero-history startup does not continue the harmonic seed's
past accelerations. It is not the physical radiation memory: it is the stock
numerical treatment of residual added mass. The coordinator explicitly excludes
an InitialOutput/prefix/feedback modification from the first bounded adapter.
Keep this exact native startup, record actual a0/net wrench and its discrepancy
from the approximate harmonic history/meanphysical mass equation, and retain the
entire transient. Initial actual acceleration is solved by native mechanics,
not assigned from q_h or cancelled by a correction force.

Therefore the proposal is not a fully stationary native solver state. Potential
small numerical startup impulse/drift from this unseeded computational state is
an explicit uncertainty. Reviewer must distinguish any physical requirement for
true fluid prehistory from an unasked exact native-discrete-periodicity goal.
No new acceleration-loop algebra, cancellation, delay removal, feedback fitting,
delay fixture or hidden startup trajectory belongs to this first proposal.

## Cheap rejection sequence before one possible body allocation

1. Source/API review of actual primitive frames, radiation force-history line,
   sample time, signs and delay parameters. No physical correction or new solve.
2. Bounded fixture-only regression: odd/asymmetric kernel with nonzero past,
   v0 endpoint and first three force rows; independently indexed direct sum must
   equal the COMBINED native CIC/PCHIP result at main/minor times. Include
   native joint pose/CGvelocity/worldomega initialization readbacks in this ONE
   bounded fixture. No delay-state feature, body trajectory/ODE replacement or
   second control is proposed.
3. Stored-data preflight: actual initial body pose/v, radiation and delayed-Ar
   histories and initial wrench residual under rotated Solid/gyro versus the
   harmonic seed, clearly separating interpolation/finite-angle/discrete effects.
   Preserve nonzero residuals. Coordinator decides whether the approximation is
   credible; no fitted compensation force or acceptance threshold invented here.
4. ONLY if separately allocated, ONE same full60s native free trajectory with
   all six actual force/history/initial state readbacks, full transient/drift,
   geometry review and unchanged coefficient warnings. No automatic second run.

This proposal does not ask for another frequency/mesh campaign, canonical sea
change, nonlinear adapter, custom ODE, prescribed body motion or damping tune.
Stopping point is this readable math/API proposal and independent review.

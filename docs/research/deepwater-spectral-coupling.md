# Deep-water near/far coupling: primary evidence and proposed first falsifier

Research record, 2026-10-22. This is an architecture proposal, not an implemented
solver or a status backlog. GitHub issue #30 owns the speed/accuracy work.

The target remains the canonical 24 ft, 2500 kg boat in the observed SoCal sea,
with nonlinear moving wetted geometry and 6DOF. The current 12 x 6 x 2 m
impermeable diagnostic tank does not represent that deep-water/open-water case.
The requested accuracy applies over 60 s, not the existing short screens.

## Primary evidence

Shi, Zhu and Xu's SCBEM separates the current surface potential between a broad
HOS field and a compact BEM correction. Their sequential Spectral-Green method
assumes the correction's exterior trace can be neglected. Green-Spectral instead
neglects the HOS contribution to body-normal velocity. The paper explicitly
identifies iterating the two components as a route to satisfying the complete
boundary conditions. Its reported wavelength/transition-width restriction is a
test-dependent limitation of the sequential approximation, not a universal
minimum pocket size. The primary article was accessible as indexed full text;
direct publisher opening returned HTTP403.

- [A spectral coupled boundary element method for the simulation of nonlinear
  surface gravity waves](https://doi.org/10.1016/j.joes.2023.07.003), online2023,
  Journal of Ocean Engineering and Science10(1),2025,123-135.

A subsequent moving-ship study uses SCBEM, acceleration-potential6DOF coupling,
overset geometry and instantaneous mesh cutting. Its reported KCS heave/pitch
and parametric-roll results establish relevant precedent, not qualification of
our boat or our discretization.

- [Simulation of ship motion response in head seas using a fully nonlinear
  potential flow method](https://doi.org/10.1016/j.oceaneng.2024.118309),2024.

The2024 HBEM paper is a different method: perturbations of free-surface and
wavemaker motions produce fixed-reference BVPs; quadratic BEM and pFFT accelerate
their solution. Its abstract describes special treatment of high derivatives
and C1 potential continuity at the free-surface/wavemaker intersection. It does
not itself establish incident-only HOS feeding as two-way moving-body coupling.
The author preprint record/abstract was accessible; full manuscript details were
not verified here.

- [HBEM publisher record](https://doi.org/10.1016/j.cma.2024.117396).
- [Author manuscript record](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4847720).

## Reuse the existing component

`tests/experiments/free_surface/hos_graph.py` already implements a periodic,
deep-water harmonic graph Dirichlet-to-normal-velocity approximation. Its
surface Taylor recursion reconstructs the harmonic trace order by order, with
padded products. It has manufactured finite-amplitude harmonic-field tests.
It has no hull, free-surface hole, radiation coupling or open-boundary treatment.
Reuse that owner; do not recreate a second HOS implementation.

## Proposed cheapest physical checkpoint (derived here)

Use the canonical flat-wetted hull, existing odd-mirror halfspace BEM and evolving
deep-water FFT surface states. At each stage write

\[
\phi=H\psi+R\sigma,\qquad
\partial_nR\sigma=V_{body}\cdot n-\partial_nH\psi.
\]

Here H is the harmonic extension of the current surface trace. R uses the
odd-mirror Green function, so its Dirichlet trace on z=0 is zero. Then advance

\[
\eta_t=|D|\psi+\partial_zR\sigma,\qquad \psi_t=-g\eta.
\]

The hull's radiated flux changes the evolving surface state; that state returns
through H to the hull. Body pressure/acceleration must use the total field and
its time derivative. This is two-way coupling, but deliberately only a flat
surface/small-motion checkpoint. It neither supplies nonlinear instantaneous
contact geometry nor establishes the final60s accuracy.

First falsifier: prescribe small heave and roll impulses of the actual hull.
Check the existing infinite-frequency added-mass limit, force, radiated flux and
returning pressure while independently refining far-grid spacing, domain size
and residual absorption. Quantify reflected returning load. An outgoing wave
picture or a body-free wave packet alone is insufficient.

The halfspace source correction must be evaluated on the exterior free surface,
not arbitrarily set to zero beyond the local patch. Filling the hidden plane
inside the hull is a numerical extension, not an extra physical free surface;
the total exterior solution and body condition must be insensitive to that
extension under resolution refinement.

## Nonlinear successor: proposed algebra, not a claimed published implementation

A simultaneous stage can be organized as

\[
A\sigma+B a=b_N,\qquad a+C\sigma=b_F,
\]

where a owns a far spectral trace with a defined extension into the near region,
C evaluates the local correction on the far free surface, and B returns the far
harmonic field to the near/free/body equations. Elimination gives

\[
(A-BC)\sigma=b_N-Bb_F.
\]

The precise trace extension, independent unknown counts and overlap constraints
must be specified and checked for rank before implementation. These schematic
blocks are not permission to prescribe both potential and normal flux
independently or to assume a windowed potential is harmonic. A local correction
has a global tail. Keeping C and B is what distinguishes this proposal from
discarding one direction of the coupling.

The practical first question is the cost and accuracy of applying C and B,
including spectral extension below the graph. Dense source-to-far-grid storage
may be prohibitive; price it on the flat checkpoint before building nonlinear
geometry coupling. No speed estimate is established here.

## Scope and acceptance cautions

- Analytic incident lifting avoids fitting the incoming long wavelengths on a
  small periodic box. It does not remove long-range hull radiation, generated
  short waves or nonlinear incident/scattered cross terms. The measured384-mode
  linear realization is not by itself an exact nonlinear free-surface solution.
- The final near field must keep actual nonlinear wet geometry and6DOF. A linear
  exterior or finite HOS order is an explicit approximation requiring order,
  interface and force/motion refinement checks.
- Deep-water |D| eliminates the fictitious2m bed. Periodic FFT wraparound remains;
  use a sufficiently large residual domain and measured absorption/reflection
  control. Do not damp the analytic incoming lift twice.
- The incident wavelength contract does not justify discarding body-generated
  shorter waves. Their effect must be assessed through loads and motion, not
  incident spectral energy alone.
- Across any explicit interface, match potential and opposite normal flux;
  pressure checks require the correct time derivative and mesh/frame chain rule.
  Transfer operators should be dual with respect to their discrete work pairing.
  A new exact discrete Hamiltonian is not a user requirement, but neither are
  unmeasured force/flux discontinuities acceptable.

No source implementation or numerical experiment was performed for this record.

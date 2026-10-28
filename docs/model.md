# Model assumptions and limitations

## Forces, motion and reconstruction

Capytaine computes hydrodynamic responses about one fixed mean hull, including
its equilibrium pose, draft and waterline. The maintained WEC-Sim fork computes
wave excitation through its public spectral inputs and integrates all six body
motions. Topwater supplies wave components and coefficients, not applied force
histories. Full dry inertia, restoring coupling, added mass and radiation memory
remain part of the model; no damping is fitted to force settling.

Motion does not rebuild the mean-hull hydrodynamics. After motion, Topwater
reconstructs incident, diffracted and radiated water elevation without feedback
to the trajectory. The rendered moving waterline does not change the forces.
This is a linear mean-geometry approximation, not a nonlinear moving-boundary
fluid solve. Breaking, slamming and flooding are outside its intended scope.

The current hull input supports convex open boats with one inferred opening.
The example is a 24 ft by 9 ft open V-hull, mass 2,500 kg, with a full inertia
tensor incorporating hull/structure and amidships/stern loads. Other supported
boats are supplied through data, not special solver branches.

## Sea and finite packets

The packaged catalog contains nine observed Southern California spectra from
NOAA NDBC station 46086. Its JSON retains source URLs/hashes and observed record
identities. Synthetic frequencies, directions and uniform random phases realize
the measured spectrum; they do not replay the original measured elevation history.
Observed support is 0.03-0.495 Hz. When a beam is supplied, the deep-water generator
rejects positive-energy bands above 0.6 Hz or with wavelengths shorter than twice
that beam; it does not crop or renormalize an incompatible observation. CLI sea
generation does not supply a boat beam. Solve checks the incoming carriers against
the actual boat's operating envelope.

Solve constructs finite quintic packets at 10 ms sampling. It retains DC and
conjugate Fourier pairs and searches the first certified eligible cutoff at or
above the carrier. The only packet search criterion is absolute exceedance of
the target envelope, with a shared 0.01 m origin budget allocated by carrier
amplitude. This is not a dwell-waveform-fit or whole-hull calmness certificate.
Equal frequency/heading components combine coherently. Each sideband keeps its
heading and uses its own dispersion. Generated sidebands can exceed the original
carrier operating band; warnings and coefficient coverage do not prove mesh accuracy.

Packet DC is a constant sea-level offset. Computational BEM geometry and stiffness
stay fixed. Native initial heave supplies its restoring load, and the offset is
added once to exported physical motion and water elevation.

## Surface policy

Radiation uses the original saved 10 ms velocities and their complete complex
transfer response, with one unpadded FFT and no oversampling or invented future.
This finite periodic construction has no general causality or wraparound guarantee.

For a target inside the mean-waterline footprint, only diffraction/radiation are
continued from its closest point ON the mean-waterline boundary. Incident waves
remain evaluated at the actual target. Strictly interior Green-function queries
are prohibited. This continuation is an engineering approximation; it does not
solve the newly exposed moving fluid region or guarantee smooth nearest-edge changes.

Normal solve retains finite/shape/clock/connectivity checks but does not perform
moving-waterline, wetted-area, disturbance-statistics or opening-clearance checks.
Render computes the moving waterline for display. A successful solve does not
establish safe opening clearance or physical applicability.

## Interpreting results

Mesh/coefficient resolution, irregular frequencies, field/load consistency and
fixed-mean applicability remain limitations to assess for the intended use.
Startup, applied-force, full-inertia, equation-of-motion and radiation-memory
audits check implementation consistency; they are not physical-error bounds.

An earlier 60 s case had a selected-model numerical elevation difference plus
reference uncertainty of 5.784 mm under a stated high-frequency response envelope.
That historical bound does NOT qualify the current 130 s packet or other inputs.
The current example's first/last water frames are within roughly 1.1 mm of calm,
while full rest intervals reach about 14 cm as incident water crosses the patch.
Endpoints alone do not establish calm intervals or numerical accuracy.

Physical validity is left to user interpretation of these assumptions and the
computed results. No automated physical qualification or CI merge gate is imposed.

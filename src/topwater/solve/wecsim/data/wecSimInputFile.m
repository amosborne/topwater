% Generic case inputs; all physical arrays and paths come from prepared-case.mat.
prepared = load(fullfile(pwd,'prepared-case.mat'));
n = prepared.n;
simu = simulationClass();
simu.simMechanicsFile = [prepared.model_name '.slx'];
simu.mode = 'normal';
simu.explorer = 'off';
simu.startTime = 0;
simu.endTime = prepared.time(end);
simu.rampTime = 0;
simu.solver = 'ode4';
simu.dt = n.time_step_s;
simu.cicDt = n.time_step_s;
simu.cicEndTime = 10;
simu.rho = n.rho;
simu.gravity = n.g;
simu.saveWorkspace = 0;
simu.paraview.option = 0;
% One coherent FullDir object; native WEC-Sim owns the harmonic force sum.
waves = waveClass('spectrumImportFullDir');
waves.spectrumFile = prepared.spectral_file;
assert(exist('nativeHydro','var') == 1,'Prepared native hydro structure is missing');
body(1) = bodyClass(nativeHydro);
body(1).geometryFile = fullfile(pwd,'boat_cg.stl');
body(1).mass = n.M(1,1);
body(1).inertia = diag(n.M(4:6,4:6))';
body(1).inertiaProducts = [n.M(4,5) n.M(4,6) n.M(5,6)];
% Computational mean water is zero. The static packet mean is represented by
% an equal/opposite initial heave; Python restores the physical world offset.
body(1).initial.displacement = n.initial_q(1:3);
body(1).initial.axis = [0 1 0];
body(1).initial.angle = 0;
body(1).nonlinearHydro = 0;
body(1).yaw.option = 0;
body(1).linearDamping = zeros(6);
body(1).quadDrag.drag = zeros(6);
body(1).quadDrag.cd = zeros(1,6);
body(1).quadDrag.area = zeros(1,6);
body(1).morisonElement.option = 0;
constraint(1) = constraintClass('Free6DOF');
constraint(1).location = [0 0 0];

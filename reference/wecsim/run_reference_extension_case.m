clear body1_out constraint1_out tout mildCicTimes mildCicForces mildCicRadiation
root = 'C:/Users/alex/topwater';
maxNumCompThreads(1);
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
caseDir = fullfile(root,'artifacts/reference/reference-extension-20261024/case-attempt1');
assert(~isfolder(caseDir),'Use a fresh case output directory');
mkdir(caseDir);
mkdir(fullfile(caseDir,'geometry'));
prepare_reference_extension_inputs(caseDir);
copyfile(fullfile(root,'reference/wecsim/wecSimInputFile_reference_extension.m'),fullfile(caseDir,'wecSimInputFile.m'));
n = load(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/native.mat'));
% native STL is already in mean world orientation; subtract CG once for the
% File Solid local frame, whose setup transform adds the mean world CG once.
tr = stlread(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/boat.stl'));
stlwrite(triangulation(tr.ConnectivityList,tr.Points-n.cg(:)'), ...
    fullfile(caseDir,'geometry/boat_cg.stl'),'text');
cd(caseDir);
preprocess_clock = tic;
wecSimInputFile;
build_reference_model(caseDir,true);
initializeWecSim;
build_reference_model(caseDir,false);
load_system(fullfile(caseDir,'minimum_reference.slx'));
preparedInput = load(fullfile(caseDir,'prepared-case.mat'));
dispatch = inspect_reference_extension_dispatch('minimum_reference');
% Make only this case copy editable for sink-only diagnostic logging.
set_param('minimum_reference/Boat','LinkStatus','inactive');
log_mild_packet_cic('minimum_reference');
preprocess_elapsed_s = toc(preprocess_clock);
% Actual H5 -> body import and preprocessing, independently reviewed metadata.
assert(max(abs(body(1).centerGravity(:)-n.cg(:))) < 1e-12);
assert(max(abs(body(1).centerBuoyancy(:)-n.cb(:))) < 1e-12);
assert(abs(body(1).volume-n.volume) < 1e-12);
assert(max(abs(body(1).hydroData.hydro_coeffs.added_mass.inf_freq*n.rho-n.Ainf),[],'all') < 1e-9);
assert(max(abs(body(1).hydroData.hydro_coeffs.linear_restoring_stiffness*n.rho*n.g-n.K),[],'all') < 1e-9);
assert(isequal(body(1).initial.displacement,preparedInput.initial_translation));
assert(isequal(body(1).initial.axis,preparedInput.initial_axis));
assert(body(1).initial.angle == preparedInput.initial_angle);
assert(simu.rampTime == 0 && simu.startTime == 0 && simu.endTime == 60 && simu.cicEndTime == 10);
expected = struct('time',preparedInput.time,'force',preparedInput.force);
native_prepared = body(1).hydroForce.hf1.userDefinedFe;
prepared = native_prepared;
prepared_error = max(abs(prepared-expected.force),[],1);
force_tolerance = 1e-9+1e-10*max(abs(expected.force),[],1);
assert(all(prepared_error <= force_tolerance),'Case force table differs from retained original threepacket excitation');
compile_clock = tic;
set_param('minimum_reference','SimulationCommand','update');
compile_elapsed_s = toc(compile_clock);
block = 'minimum_reference/Boat/Hydrodynamic Body/Structure/Body Properties';
expression = get_param(block,'ProductsOfInertia');
products = eval(expression);
hf = body(1).hydroForce.hf1;
solid = zeros(6);
solid(1:3,1:3) = hf.adjustedMass*eye(3);
solid(4:6,4:6) = diag(hf.adjustedInertia);
solid(5,6) = products(1); solid(6,5) = products(1);
solid(4,6) = products(2); solid(6,4) = products(2);
solid(4,5) = products(3); solid(5,4) = products(3);
mass_error = max(abs(solid+hf.fAddedMass-n.M-n.Ainf),[],'all');
assert(mass_error < 1e-8,'Compiled selected block dry/full added tensor mapping differs');
solid_world_initial = solid;
solid_world_initial(4:6,4:6) = preparedInput.initial_rotation*solid(4:6,4:6)*preparedInput.initial_rotation';
effective_world_initial = solid_world_initial+hf.fAddedMass;
preflight = struct('compiled',true,'exact_h5_metadata_readback',true, ...
    'parameter_mass_added_roundtrip_error',mass_error,'actual_products_expression',expression, ...
    'solid_local_matrix',solid,'solid_world_initial_matrix',solid_world_initial, ...
    'residual_added_mass_world_matrix',hf.fAddedMass,'effective_world_initial_matrix',effective_world_initial, ...
    'mass_scope','parameter roundtrip only; Solid rotates with initial pose, native residual acceleration feedback is world resolved', ...
    'products_yz_zx_xy',products,'mean_world_cg',n.cg(:)', ...
    'dispatch',dispatch, ...
    'local_stl','native mean-world STL minus mean CG once', ...
    'initial_rotation','identity perturbation; no repeated equilibrium quaternion');
fid = fopen(fullfile(caseDir,'preflight.json'),'w');
fprintf(fid,'%s\n',jsonencode(preflight,PrettyPrint=true));
fclose(fid);
try
    dynamics_clock = tic;
    sim(simu.simMechanicsFile,[],simset('SrcWorkspace','base'));
    dynamics_elapsed_s = toc(dynamics_clock);
catch simulationFailure
    if exist('body1_out','var')
        save(fullfile(caseDir,'native-failure.mat'),'body1_out','-v7');
    end
    fid = fopen(fullfile(caseDir,'simulation-failure.txt'),'w');
    fprintf(fid,'%s\n',getReport(simulationFailure,'extended','hyperlinks','off'));
    fclose(fid);
    diary off
    rethrow(simulationFailure)
end
% Preserve native bus and grid directly. The scalar-wave stock postprocessor
% is not required for this original threepacket extension's native output.
native_body = body1_out;
save(fullfile(caseDir,'native-output.mat'),'native_body','-v7');
save(fullfile(caseDir,'native-cic-output.mat'),'mildCicTimes','mildCicForces','mildCicRadiation','-v7');
assert(all(isfinite(native_body.time),'all') && all(isfinite(native_body.signals.values),'all'));
assert(all(isfinite(mildCicTimes.time),'all') && all(isfinite(mildCicTimes.signals.values),'all'));
assert(all(isfinite(mildCicForces.time),'all') && all(isfinite(mildCicForces.signals.values),'all'));
assert(all(isfinite(mildCicRadiation.time),'all') && all(isfinite(mildCicRadiation.signals.values),'all'));
time = native_body.time;
position_world = native_body.signals.values(:,1:6);
velocity_world = native_body.signals.values(:,7:12);
acceleration_world = native_body.signals.values(:,13:18);
q = position_world;
q(:,1:3) = q(:,1:3)-n.cg(:)';
assert(size(q,2) == 6 && all(isfinite(q),'all'));
assert(abs(time(1)) < 1e-12 && abs(time(end)-60) < 1e-10);
expectedInitial = [preparedInput.initial_translation preparedInput.initial_euler];
assert(max(abs(q(1,:)-expectedInitial)) < 1e-8,'Native initial pose differs from independently known neutral pose');
assert(max(abs(velocity_world(1,:))) < 1e-9);
assert(max(abs(native_body.signals.values(1,31:36))) < 1e-9,'Initial radiation memory force is not zero');
assert(max(abs(diff(time)-simu.dt)) < 1e-9);
assert(numel(time) == 3001 && max(abs(time(:)-expected.time(:))) < 1e-9);
applied_excitation = native_body.signals.values(:,25:30);
expected_on_grid = interp1(expected.time(:),expected.force,time(:));
applied_error = max(abs(applied_excitation-expected_on_grid),[],1);
assert(all(applied_error <= force_tolerance),'Applied native excitation differs from accepted total history');
cic_error = verify_reference_extension_cic(native_body,hf,simu.dt,caseDir,mildCicTimes,mildCicForces,mildCicRadiation);
assert(cic_error < 1e-6,'Native zero-history CIC rows/RT/PCHIP/bus differ');
x = native_body.signals.values;
wrench_error = max(abs(x(:,19:24)-(x(:,25:30)-x(:,31:36)-x(:,37:42) - ...
    x(:,43:48)-x(:,49:54)-x(:,55:60))),[],'all');
assert(wrench_error < 1e-9);
assert(all(x(:,49:60) == 0,'all'),'Viscous/linear damping channels must remain zero');
balance_error = zeros(numel(time),6);
delay_error = zeros(numel(time),6);
for j = 1:numel(time)
    roll = q(j,4); pitch = q(j,5); yaw = q(j,6);
    Rx = [1 0 0;0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
    Ry = [cos(pitch) 0 sin(pitch);0 1 0;-sin(pitch) 0 cos(pitch)];
    Rz = [cos(yaw) -sin(yaw) 0;sin(yaw) cos(yaw) 0;0 0 1];
    R = Rz*Ry*Rx;
    I = R*solid(4:6,4:6)*R';
    omega = velocity_world(j,4:6)';
    expected_force = [(solid(1:3,1:3)*acceleration_world(j,1:3)')' ...
        (I*acceleration_world(j,4:6)'+cross(omega,I*omega))'];
    balance_error(j,:) = expected_force-native_body.signals.values(j,19:24);
    if j == 1
        delayed = zeros(6,1);
    elseif j == 2
        delayed = acceleration_world(1,:)';
    else
        r = 1-1e-7/simu.dt;
        delayed = ((1+r)*acceleration_world(j-1,:)-r*acceleration_world(j-2,:))';
    end
    delay_error(j,:) = (hf.fAddedMass*delayed)'-native_body.signals.values(j,37:42);
end
save(fullfile(caseDir,'native-balance-check.mat'),'balance_error','delay_error','-v7');
assert(max(abs(balance_error),[],'all') < 1e-6,'Native Solid Newton-Euler balance failed');
assert(max(abs(delay_error),[],'all') < 1e-6,'Native residual acceleration feedback differs');
assert(any(max(abs(q-expectedInitial),[],1) > 1e-8),'No resolved motion after release relative to known initial pose');
excitation_prepared = prepared;
mean_world_cg = n.cg(:)';
save(fullfile(caseDir,'motion.mat'),'native_body','time','q','position_world', ...
    'velocity_world','acceleration_world','excitation_prepared','mean_world_cg', ...
    'applied_excitation','expected_on_grid','force_tolerance','-v7');
save(fullfile(caseDir,'native-prepared-excitation.mat'),'native_prepared','dispatch','-v7');
writematrix([time q],fullfile(caseDir,'motion.csv'));
oldMotion = load(fullfile(root,'artifacts/reference/minimum-20261024/case-attempt1/motion.mat'));
repeatability_q_error = q(1:501,:)-oldMotion.q;
repeatability_v_error = velocity_world(1:501,:)-oldMotion.velocity_world;
save(fullfile(caseDir,'first10-repeatability.mat'),'repeatability_q_error','repeatability_v_error','-v7');
summary = struct('first10_q_max_error',max(abs(repeatability_q_error),[],1), ...
    'first10_v_max_error',max(abs(repeatability_v_error),[],1),'complete',true,'duration_s',time(end)-time(1), ...
    'preprocess_elapsed_s',preprocess_elapsed_s,'compile_elapsed_s',compile_elapsed_s, ...
    'dynamics_elapsed_s',dynamics_elapsed_s, ...
    'samples',numel(time),'dt_s',simu.dt,'initial_q',q(1,:), ...
    'prepared_excitation_error_by_axis',prepared_error,'applied_excitation_error_by_axis',applied_error, ...
    'excitation_tolerance_by_axis',force_tolerance, ...
    'initial_velocity',velocity_world(1,:),'maximum_abs_motion',max(abs(q),[],1), ...
    'newton_euler_balance_max_by_axis',max(abs(balance_error),[],1), ...
    'residual_added_mass_delay_error_by_axis',max(abs(delay_error),[],1), ...
    'peak_to_peak_motion',max(q,[],1)-min(q,[],1), ...
    'maximum_abs_change_from_known_initial',max(abs(q-expectedInitial),[],1), ...
    'initial_acceleration',acceleration_world(1,:), ...
    'initial_total_force',native_body.signals.values(1,19:24), ...
    'initial_excitation_force',applied_excitation(1,:), ...
    'initial_radiation_force',native_body.signals.values(1,31:36), ...
    'initial_restoring_force',native_body.signals.values(1,43:48), ...
    'initial_mean_matrix_acceleration_diagnostic',((n.M+n.Ainf)*acceleration_world(1,:)')', ...
    'mean_world_cg',mean_world_cg,'translation','world CG minus mean world CG, m', ...
    'angular_velocity','world angular velocity omega; not Euler angle rates', ...
    'rotation','native quaternion2EulXYZ output [roll pitch yaw], rad; R=Rz(yaw)*Ry(pitch)*Rx(roll)', ...
    'initial_orientation','q0=0 in initially numerical-local-calm field; equilibrium quaternion already baked into STL and inertia', ...
    'forcing','unchanged accepted REF001 three incoming packets; futurepadding/duration only; diagnostic', ...
    'radiation_initial_memory','zero from idealized stationary held-body past; no prerun', ...
    'native_cic_error',cic_error,'initial_global_heading_rad',0, ...
    'heading_offset_deg',0,'status','ORIGINAL REF001 60s duration extension; unqualified17panel/CIC10; not finalbuoysea', ...
    'free_joint','all six DOFs free; no PTO/mooring/empirical damping');
fid = fopen(fullfile(caseDir,'summary.json'),'w');
fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));
fclose(fid);
save_system('minimum_reference');
close_system('minimum_reference',0);
diary off
cd(root);
fprintf('Native original REF001 packet extension completed for60s.\n');

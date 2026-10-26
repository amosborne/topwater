% Called in the existing shared engine's BASE workspace; never launches MATLAB.
% topwaterCaseDir and topwaterFork are supplied by the Python adapter.
clear body1_out topwaterCicTimes topwaterCicForces topwaterCicRadiation
assert(exist('topwaterEngineGuard','var') == 1, ...
    'Run through the shared-engine guard before creating a model');
cd(topwaterCaseDir);
maxNumCompThreads(1);
addpath(genpath(fullfile(topwaterFork,'source')));
prepared = load('prepared-case.mat');
n = prepared.n;
[h,nativeHydro] = prepare_workflow_hydrodynamics(prepared);
wecSimInputFile;
build_workflow_model(topwaterCaseDir,prepared.model_name,true);
initializeWecSim;
prepared = load('prepared-case.mat');
n = prepared.n;
native_wave_time = waves(1).waveAmpTime(:,1);
native_excitation_time = prepared.time(:);
assert(max(abs(native_wave_time-native_excitation_time)) < 1e-9);
native_origin_elevation = waves.waveAmpTime(:,2);
expected_origin_elevation = prepared.incident_origin_elevation_m(:);
phase_max_abs_difference_m = max(abs(native_origin_elevation-expected_origin_elevation));
save('native-phase.mat','native_excitation_time','native_wave_time','native_origin_elevation', ...
    'expected_origin_elevation','phase_max_abs_difference_m','-v7');
assert(phase_max_abs_difference_m <= 1e-12, ...
    'Native spectral incident elevation differs from the unchanged physical packet');
native_excitation = native_spectral_excitation(body(1),waves,prepared.time);
assert(isequal(size(native_excitation),[numel(prepared.time) 6]) && ...
    all(isfinite(native_excitation),'all'),'Native spectral excitation is incomplete');
save('native-excitation.mat','native_excitation_time','native_wave_time','native_excitation','-v7');
startup = native_excitation_startup(prepared,native_excitation);
save('native-excitation.mat','startup','-append');
if prepared.force_only
    % Preprocessing-only comparison: no body model or trajectory is created.
    return
end
build_workflow_model(topwaterCaseDir,prepared.model_name,false);
load_system(fullfile(topwaterCaseDir,[prepared.model_name '.slx']));
inspect_native_excitation(prepared.model_name);
configure_cic_logging(prepared.model_name);
assert(max(abs(body(1).centerGravity(:)-n.cg(:))) < 1e-12);
assert(max(abs(body(1).centerBuoyancy(:)-n.cb(:))) < 1e-12);
assert(abs(body(1).volume-n.volume) < 1e-12);
assert(max(abs(body(1).hydroData.hydro_coeffs.added_mass.inf_freq*n.rho-n.Ainf),[],'all') < 1e-9);
assert(max(abs(body(1).hydroData.hydro_coeffs.linear_restoring_stiffness*n.rho*n.g-n.K),[],'all') < 1e-9);
set_param(prepared.model_name,'SimulationCommand','update');
block = [prepared.model_name '/Boat/Hydrodynamic Body/Structure/Body Properties'];
products = eval(get_param(block,'ProductsOfInertia'));
hf = body(1).hydroForce.hf1;
solid = zeros(6);
solid(1:3,1:3) = hf.adjustedMass*eye(3);
solid(4:6,4:6) = diag(hf.adjustedInertia);
solid(5,6) = products(1); solid(6,5) = products(1);
solid(4,6) = products(2); solid(6,4) = products(2);
solid(4,5) = products(3); solid(5,4) = products(3);
assert(max(abs(solid+hf.fAddedMass-n.M-n.Ainf),[],'all') < 1e-8, ...
    'Compiled dry/full added-mass mapping differs from physical input');
sim(simu.simMechanicsFile,[],simset('SrcWorkspace','base'));
native_body = body1_out;
time = native_body.time(:);
x = native_body.signals.values;
assert(size(x,2) == 60 && all(isfinite(x),'all') && all(isfinite(time)));
assert(numel(time) == numel(prepared.time) && max(abs(time-prepared.time)) < 1e-9);
assert(abs(time(1)) < 1e-12 && abs(time(end)-prepared.time(end)) < 1e-9);
q = x(:,1:6); q(:,1:3) = q(:,1:3)-n.cg(:)';
assert(max(abs(q(1,:))) < 1e-8 && max(abs(x(1,7:12))) < 1e-9, ...
    'Actual native initial state is not the prepared equilibrium at rest');
assert(max(abs(x(1,31:36))) < 1e-9,'Initial radiation memory is not zero');
tolerance = 1e-10*max(abs(native_excitation),[],1);
assert(all(max(abs(x(:,25:30)-native_excitation),[],1) <= tolerance), ...
    'Actual applied excitation differs from WEC-Sim native spectral preparation');
assert(max(abs(x(:,19:24)-(x(:,25:30)-x(:,31:36)-x(:,37:42)-x(:,43:48)-x(:,49:54)-x(:,55:60))),[],'all') < 1e-9);
assert(all(x(:,49:60) == 0,'all'),'Unexpected viscous or linear damping');
verify_balance(x,q,solid,hf.fAddedMass,simu.dt);
verify_cic(time,x(:,7:12),x(:,31:36),hf.irkb,simu.dt, ...
    topwaterCicTimes,topwaterCicForces,topwaterCicRadiation);
velocity = x(:,7:12);
acceleration = x(:,13:18);
excitation = x(:,25:30);
radiation = x(:,31:36);
save('motion.mat','time','q','velocity','acceleration','excitation','radiation','-v7');
close_owned_model(prepared.model_name);

function inspect_native_excitation(model)
% Read-only audit of one native FullDir harmonic branch, without rebinding.
set_param(model,'SimulationCommand','update');
parent = [model '/Boat/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation'];
children = find_system(parent,'SearchDepth',1,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.allVariants,'Type','block');
branches = {};
for k = 1:numel(children)
    if ~isempty(regexp(get_param(children{k},'Name'), ...
            '^Linear Wave Excitation Force Variant Subsystem[0-9]*$','once'))
        branches{end+1} = children{k};
    end
end
assert(numel(branches) == 1,'Expected one native FullDir excitation branch');
ports = get_param(branches{1},'PortHandles');
destination = get_param(get_param(ports.Outport(1),'Line'),'DstPortHandle');
assert(numel(destination) == 1 && strcmp(get_param(destination,'Parent'),[parent '/Add']));
signs = get_param([parent '/Add'],'Inputs');
assert(signs(get_param(destination,'PortNumber')) == '+');
active = find_system(parent,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'Type','block');
assert(any(endsWith(active,'/Irregular Wave Full Direction')), ...
    'Native FullDir excitation variant is not active');
sources = find_system(parent,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'BlockType','FromWorkspace');
for k = 1:numel(sources)
    assert(~contains(get_param(sources{k},'VariableName'),'userDefinedFe'), ...
        'An imported applied-force source is unexpectedly active');
end
end

function configure_cic_logging(model)
% Only sink logging requires an editable case copy; excitation stays native.
set_param([model '/Boat'],'LinkStatus','inactive');
parent = [model '/Boat/Hydrodynamic Body/Wave Radiation Forces Calculation'];
choices = find_system(parent,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'Name','Convolution Integral Calculation');
selected = {};
for j = 1:numel(choices)
    if getSimulinkBlockHandle([choices{j} '/MATLAB Function1']) > 0
        selected{end+1} = choices{j};
    end
end
assert(numel(selected) == 1);
cic = selected{1};
if strcmp(get_param(cic,'LinkStatus'),'resolved'), set_param(cic,'LinkStatus','inactive'); end
sources = {'Rate Transition/1','Rate Transition1/1','MATLAB Function/1'};
variables = {'topwaterCicTimes','topwaterCicForces','topwaterCicRadiation'};
for j = 1:3
    name = ['Topwater log ' num2str(j)];
    add_block('simulink/Sinks/To Workspace',[cic '/' name], ...
        'VariableName',variables{j},'SaveFormat','Structure With Time',position=[1000 50 1130 90]);
    add_line(cic,sources{j},[name '/1'],autorouting='on');
end
end

function verify_balance(x,q,solid,added,dt)
for j = 1:size(x,1)
    r = q(j,4); p = q(j,5); y = q(j,6);
    R = [cos(y) -sin(y) 0;sin(y) cos(y) 0;0 0 1]* ...
        [cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)]* ...
        [1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)];
    inertia = R*solid(4:6,4:6)*R';
    omega = x(j,10:12)';
    actual = [(solid(1:3,1:3)*x(j,13:15)')' ...
        (inertia*x(j,16:18)'+cross(omega,inertia*omega))'];
    assert(max(abs(actual-x(j,19:24))) < 1e-6,'Native Newton-Euler balance failed');
    if j == 1
        delayed = zeros(6,1);
    elseif j == 2
        delayed = x(1,13:18)';
    else
        ratio = 1-1e-7/dt;
        delayed = ((1+ratio)*x(j-1,13:18)-ratio*x(j-2,13:18))';
    end
    assert(max(abs((added*delayed)'-x(j,37:42))) < 1e-6, ...
        'Native residual added-mass delay differs from its declared discretization');
end
end

function verify_cic(time,velocity,bus,kernel,dt,timeLog,forceLog,radLog)
N = size(kernel,1)-1;
assert(N == 500 && dt == .02);
weights = dt*ones(N+1,1); weights([1 end]) = dt/2;
fullVelocity = [zeros(N+3,6);velocity];
fullForce = zeros(size(fullVelocity));
for i = 1:6
    for j = 1:6
        values = conv(fullVelocity(:,j),weights.*kernel(:,i,j));
        fullForce(:,i) = fullForce(:,i)+values(1:size(fullVelocity,1));
    end
end
assert(numel(timeLog.time) == numel(time) && numel(forceLog.time) == numel(time) && numel(radLog.time) == numel(time));
assert(max(abs(timeLog.time(:)-time)) < 1e-9 && max(abs(forceLog.time(:)-time)) < 1e-9 && max(abs(radLog.time(:)-time)) < 1e-9);
for k = 1:numel(time)
    times = read_sample(timeLog,k,[3 1]);
    rows = read_sample(forceLog,k,[3 6]);
    index = round(times/dt)+N+4;
    assert(all(index >= N+1 & index <= size(fullForce,1)));
    expected = fullForce(index,:);
    tolerance = 1e-6+1e-10*max(abs(expected),[],1);
    assert(all(max(abs(rows-expected),[],1) <= tolerance),'Native CIC history differs');
    radiation = read_sample(radLog,k,[1 6]);
    expectedRadiation = interp1(times,expected,time(k),'pchip','extrap');
    assert(all(abs(radiation-expectedRadiation) <= tolerance),'Native CIC interpolation differs');
    assert(all(abs(bus(k,:)-radiation) <= tolerance),'Native radiation output bus differs');
end
end

function value = read_sample(log,k,shape)
data = log.signals.values;
if size(data,1) == numel(log.time)
    value = reshape(data(k,:,:),shape);
else
    assert(size(data,3) == numel(log.time));
    value = reshape(data(:,:,k),shape);
end
end

function close_owned_model(name)
if bdIsLoaded(name), close_system(name,0); end
end

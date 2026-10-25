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
prepare_hydrodynamics(prepared);
wecSimInputFile;
build_workflow_model(topwaterCaseDir,prepared.model_name,true);
initializeWecSim;
prepared = load('prepared-case.mat');
n = prepared.n;
build_workflow_model(topwaterCaseDir,prepared.model_name,false);
load_system(fullfile(topwaterCaseDir,[prepared.model_name '.slx']));
topwaterForceTable = [prepared.time prepared.force];
configure_force_table(prepared.model_name);
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
tolerance = 1e-9+1e-10*max(abs(prepared.force),[],1);
assert(all(max(abs(x(:,25:30)-prepared.force),[],1) <= tolerance), ...
    'Actual native excitation differs from prepared directional incoming field');
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

function prepare_hydrodynamics(prepared)
n = prepared.n;
h = readCAPYTAINE(struct(),prepared.netcdf_path);
assert(h.Nb == 1 && all(h.dof == 6));
assert(h.rho == n.rho && h.g == n.g && isinf(h.h));
% Verify the official importer's one complex-convention conversion directly.
names = erase(lower(string(ncread(prepared.netcdf_path,'influenced_dof')')),char(0));
modes = ["surge","sway","heave","roll","pitch","yaw"];
index = zeros(1,6);
for j = 1:6, index(j) = find(names == modes(j)); end
fk = named_array(prepared.netcdf_path,'Froude_Krylov_force',{'influenced_dof','wave_direction','omega','complex'});
sc = named_array(prepared.netcdf_path,'diffraction_force',{'influenced_dof','wave_direction','omega','complex'});
components = erase(string(ncread(prepared.netcdf_path,'complex')'),char(0));
re = fk(index,:,:,components == "re")+sc(index,:,:,components == "re");
im = fk(index,:,:,components == "im")+sc(index,:,:,components == "im");
assert(max(abs(h.ex_re*n.rho*n.g-re),[],'all') <= 1e-10*max(1,max(abs(re),[],'all')));
assert(max(abs(h.ex_im*n.rho*n.g+im),[],'all') <= 1e-10*max(1,max(abs(im),[],'all')));
h.cg = n.cg(:); h.cb = n.cb(:); h.Vo = n.volume;
h.Khs = n.K/(n.rho*n.g);
h = radiationIRF(h,10,501,2401,[],[]);
h = excitationIRF(h,20,2001,2401,[],[]);
h.Ainf = n.Ainf/n.rho;
h.file = fullfile(pwd,'canonical');
writeBEMIOH5(h);
end

function value = named_array(file,name,order)
info = ncinfo(file,name);
names = {info.Dimensions.Name};
permutation = zeros(1,numel(order));
for j = 1:numel(order), permutation(j) = find(strcmp(names,order{j})); end
permutation = [permutation setdiff(1:numel(names),permutation,'stable')];
value = permute(ncread(file,name),permutation);
shape = size(value);
assert(all(shape(numel(order)+1:end) == 1));
value = reshape(value,shape(1:numel(order)));
end

function configure_force_table(model)
% Only the copied case instance receives a data-binding override. The native
% source/library and single-total imported-force dispatch remain unchanged.
set_param([model '/Boat'],'LinkStatus','inactive');
set_param(model,'SimulationCommand','update');
parent = [model '/Boat/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation'];
branches = find_system(parent,'SearchDepth',1,'LookUnderMasks','all','FollowLinks','on', ...
    'regexp','on','Name','^Linear Wave Excitation Force Variant Subsystem[0-9]*$');
assert(numel(branches) == 1,'Expected one total imported excitation branch');
ports = get_param(branches{1},'PortHandles');
destination = get_param(get_param(ports.Outport(1),'Line'),'DstPortHandle');
assert(numel(destination) == 1 && strcmp(get_param(destination,'Parent'),[parent '/Add']));
signs = get_param([parent '/Add'],'Inputs');
assert(signs(get_param(destination,'PortNumber')) == '+');
active = find_system(parent,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'BlockType','FromWorkspace');
selected = {};
for j = 1:numel(active)
    if contains(get_param(active{j},'VariableName'),'userDefinedFe')
        selected{end+1} = active{j};
    end
end
assert(numel(selected) == 1,'Expected one active positive imported-force source');
set_param(selected{1},'VariableName','topwaterForceTable');
end

function configure_cic_logging(model)
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

function run_official_pto_zero()
key = 'cic';
root = 'C:/Users/alex/topwater';
base = fullfile(root,'artifacts/reference/official-pto-zero-20261024');
assert(ismember(key,{'standard','cic'}));
config = jsondecode(fileread(fullfile(base,key,'case.json')));
fork = fullfile(root,'artifacts/dependencies/topwater-wec-sim');
addpath(genpath(fullfile(fork,'source')));
caseDir = fullfile(root,config.case_directory);
cd(caseDir);
assert(~isfile('audit.mat') && ~isfile('audit.json'));
assert(~bdIsLoaded('RM3'),'Existing RM3 model must be reconciled before this run');
assert(strcmp(strrep(which('wecSim'),'\','/'),strrep(fullfile(fork,'source/wecSim.m'),'\','/')));
set(groot,'defaultFigureVisible','off');
% Sole physical change: PTO damping0; official inputs/model/H5 remain byte-identical.
mcr.header = {'simu.explorer','pto(1).damping'};
mcr.cases = {'off',0};
imcr = 1;
clockStart = tic;
wecSim;
elapsed = toc(clockStart);
bodies = output.bodies;
ptos = output.ptos;
storage = cell(1,numel(body));
for j = 1:numel(body)
    storage{j} = body(j).hydroForce.hf1.storage;
end
settings = struct('wave_type',waves.type,'height_m',waves.height,'period_s',waves.period, ...
    'end_s',simu.endTime,'dt_s',simu.dt,'ramp_s',simu.rampTime, ...
    'cic_s',simu.cicEndTime,'cic_dt_s',simu.cicDt,'stateSpace',simu.stateSpace, ...
    'pto_damping',pto(1).damping,'pto_stiffness',pto(1).stiffness);
parameters = cell(1,numel(body));
for j = 1:numel(body)
    parameters{j} = body(j).hydroForce.hf1;
end
save('audit.mat','bodies','ptos','storage','settings','parameters','-v7');
energy = audit_official_pto_energy(bodies,ptos,storage,parameters,settings);
save('energy-budget.mat','energy','-v7');
receipt = struct('native_entry','wecSim','elapsed_s',elapsed,'settings',settings, ...
    'display_only_override','native MCR simu.explorer=off; defaultFigureVisible=off', ...
    'physical_override','pto(1).damping=0 only','golden_available',config.golden_available,'golden_max_heave_error_m',[], ...
    'energy_scope','CIC resisting work is not complete fluid energy; PTO and constraints retained');
for j = 1:numel(bodies)
    t = bodies(j).time(:);
    v = bodies(j).velocity;
    fr = bodies(j).forceRadiationDamping;
    power = sum(v.*fr,2);
    work = cumtrapz(t,power);
    receipt.body(j).samples = numel(t);
    receipt.body(j).finite = all(isfinite(bodies(j).position),'all') && all(isfinite(v),'all');
    receipt.body(j).finite_forces = all(isfinite(fr),'all') && ...
        all(isfinite(storage{j}.output_forceAddedMass),'all') && ...
        all(isfinite(storage{j}.output_forceTotal),'all');
    receipt.body(j).initial_velocity = v(1,:);
    receipt.body(j).initial_radiation = fr(1,:);
    receipt.body(j).total_resisting_work_J = work(end);
    receipt.body(j).minimum_cumulative_work_J = min(work);
    receipt.body(j).minimum_power_W = min(power);
    receipt.body(j).initial_position = bodies(j).position(1,:);
    receipt.body(j).force_sum_error = max(abs(storage{j}.output_forceTotal - ...
        (bodies(j).forceExcitation-fr-storage{j}.output_forceAddedMass - ...
        bodies(j).forceRestoring-bodies(j).forceMorisonAndViscous - ...
        bodies(j).forceLinearDamping)),[],'all');
end
receipt.combined_resisting_work_J = sum([receipt.body.total_resisting_work_J]);
receipt.force_sum_scope = 'diagnostic residual, not a force or energy acceptance tolerance';
if config.golden_available
    golden = load('regularCIC_org.mat');
    errors = zeros(1,3);
    for j = 1:2
        assert(isequal(bodies(j).time,golden.output.bodies(j).time));
        actual = bodies(j).position(:,3)-bodies(j).position(1,3);
        expected = golden.output.bodies(j).position(:,3)-golden.output.bodies(j).position(1,3);
        errors(j) = max(abs(actual-expected));
    end
    actual = bodies(1).position(:,3)-bodies(1).position(1,3) - ...
        bodies(2).position(:,3)+bodies(2).position(1,3);
    expected = golden.output.bodies(1).position(:,3)-golden.output.bodies(1).position(1,3) - ...
        golden.output.bodies(2).position(:,3)+golden.output.bodies(2).position(1,3);
    errors(3) = max(abs(actual-expected));
    receipt.golden_max_heave_error_m = errors;
end
fid = fopen('audit.json','w');
fprintf(fid,'%s\n',jsonencode(receipt,PrettyPrint=true));
fclose(fid);
assert(settings.ramp_s == 100 && settings.dt_s == .1);
assert(settings.pto_damping == 0 && settings.pto_stiffness == 0);
assert(settings.stateSpace == 0);
if strcmp(key,'cic')
    assert(strcmp(settings.wave_type,'regularCIC') && settings.cic_s == 20 && settings.cic_dt_s == .1);
else
    assert(strcmp(settings.wave_type,'regular'));
end
for j = 1:numel(bodies)
    assert(receipt.body(j).finite);
    assert(receipt.body(j).finite_forces);
    assert(max(abs(receipt.body(j).initial_velocity)) < 1e-10);
    assert(max(abs(receipt.body(j).initial_radiation)) < 1e-10);
    assert(abs(bodies(j).time(end)-config.duration_s) < 1e-9);
    expectedGrid = (0:.1:config.duration_s)';
    assert(numel(bodies(j).time) == numel(expectedGrid));
    assert(max(abs(bodies(j).time(:)-expectedGrid)) < 1e-9);
end
if config.golden_available
    assert(all(errors <= 1e-5),'Official golden heave tolerance failed; preserve all outputs');
end
close_system('RM3',0);
diary off
cd(root);
end

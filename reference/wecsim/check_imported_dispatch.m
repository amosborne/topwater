root = 'C:/Users/alex/topwater';
maxNumCompThreads(1);
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
caseDir = fullfile(root,'artifacts/reference-imported-dispatch-fixed-attempt2');
assert(~isfolder(caseDir));
mkdir(caseDir);
mkdir(fullfile(caseDir,'geometry'));
copyfile(fullfile(root,'reference/wecsim/wecSimInputFile.m'),caseDir);
n = load(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/native.mat'));
tr = stlread(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/boat.stl'));
stlwrite(triangulation(tr.ConnectivityList,tr.Points-n.cg(:)'), ...
    fullfile(caseDir,'geometry/boat_cg.stl'),'text');
cd(caseDir);
if bdIsLoaded('minimum_reference')
    assert(strcmp(get_param('minimum_reference','SimulationStatus'),'stopped'));
    close_system('minimum_reference',0);
end
clear body waves simu
wecSimInputFile;
build_reference_model(caseDir,true);
initializeWecSim;
build_reference_model(caseDir,false);
load_system(fullfile(caseDir,'minimum_reference.slx'));
set_param('minimum_reference','SimulationCommand','update');
parent = 'minimum_reference/Boat/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation';
fresh = dispatch_graph(parent);
assert(fresh.branch_count == 1 && fresh.active_imported == 1);
importedWaves = waves;
prepared = body(1).hydroForce.hf1.userDefinedFe;
% Repeat initialization with the same physical object group, changed expression.
set_param(parent,'waveGroup','waves(1:3)');
set_param('minimum_reference','SimulationCommand','update');
repeated = dispatch_graph(parent);
assert(repeated.branch_count == 1 && repeated.active_imported == 1);
% Ordinary three regular waves retain the original native replication behavior.
headings = [0 45 90];
for j = 1:3
    waves(j) = waveClass('regularCIC');
    waves(j).height = 0.02;
    waves(j).period = 4;
    waves(j).direction = headings(j);
    waves(j).setup(body(1).hydroData.simulation_parameters.w, ...
        body(1).hydroData.simulation_parameters.waterDepth,0,0.02,500,0:0.02:10,n.g,n.rho);
end
regularWaves = waves;
typeNum = 11;
set_param(parent,'waveGroup','waves');
set_param('minimum_reference','SimulationCommand','update');
regular = dispatch_graph(parent);
assert(regular.branch_count == 3 && regular.active_imported == 0);
% Returning to imports must remove persisted extra branches and their sum ports.
waves = importedWaves;
typeNum = 40;
set_param(parent,'waveGroup','waves(1:3)');
set_param('minimum_reference','SimulationCommand','update');
restored = dispatch_graph(parent);
assert(restored.branch_count == 1 && restored.active_imported == 1);
assert(strcmp(restored.sum_inputs,fresh.sum_inputs));
% Mixed imported/non-imported groups are unsupported and fail before editing.
waves = [importedWaves(1),regularWaves(2:3)];
rejected = false;
try
    set_param(parent,'waveGroup','waves');
    set_param('minimum_reference','SimulationCommand','update');
catch e
    rejected = contains(getReport(e,'extended','hyperlinks','off'),'Cannot mix elevationImport');
end
assert(rejected,'Mixed imported/non-imported group was not explicitly rejected');
waves = importedWaves;
typeNum = 40;
set_param(parent,'waveGroup','waves(1:3)');
set_param('minimum_reference','SimulationCommand','update');
final = dispatch_graph(parent);
assert(final.branch_count == 1 && final.active_imported == 1);
assert(isequal(prepared,body(1).hydroForce.hf1.userDefinedFe));
result = struct('complete',true,'fresh',fresh,'repeated',repeated, ...
    'regular',regular,'restored',restored,'mixed_rejected',rejected,'final',final, ...
    'precomputed_total_unchanged',true,'body_simulated',false, ...
    'result','Actual native imported total dispatched once, including stale-clone cleanup');
fid = fopen(fullfile(caseDir,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fclose(fid);
save_system('minimum_reference');
close_system('minimum_reference',0);
diary off
cd(root);
fprintf('Actual imported dispatch fresh/repeat/regular-transition/mixed regressions passed.\n');

function r = dispatch_graph(parent)
branches = find_system(parent,'SearchDepth',1,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.allVariants,'regexp','on', ...
    'Name','^Linear Wave Excitation Force Variant Subsystem[0-9]*$');
ports = zeros(1,numel(branches));
for i = 1:numel(branches)
    p = get_param(branches{i},'PortHandles');
    destination = get_param(get_param(p.Outport(1),'Line'),'DstPortHandle');
    assert(numel(destination) == 1 && strcmp(get_param(destination,'Parent'),[parent '/Add']));
    ports(i) = get_param(destination,'PortNumber');
end
signs = get_param([parent '/Add'],'Inputs');
assert(all(signs(ports) == '+') && numel(unique(ports)) == numel(ports));
active = find_system(parent,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'BlockType','FromWorkspace');
expressions = {};
for i = 1:numel(active)
    value = get_param(active{i},'VariableName');
    if contains(value,'userDefinedFe'), expressions{end+1} = value; end
end
r = struct('branch_count',numel(branches),'positive_sum_ports',ports, ...
    'sum_inputs',signs,'active_imported',numel(expressions),'active_total_expressions',{expressions});
end

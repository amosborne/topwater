root = 'C:/Users/alex/topwater';
maxNumCompThreads(1);
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
caseDir = fullfile(root,'artifacts/reference-imported-dispatch-stock-attempt4');
assert(~isfolder(caseDir));
mkdir(caseDir);
mkdir(fullfile(caseDir,'geometry'));
copyfile(fullfile(root,'reference/wecsim/wecSimInputFile.m'),caseDir);
n = load(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/native.mat'));
tr = stlread(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/boat.stl'));
stlwrite(triangulation(tr.ConnectivityList,tr.Points-n.cg(:)'), ...
    fullfile(caseDir,'geometry/boat_cg.stl'),'text');
cd(caseDir);
% Only this test's earlier failed partial model may be discarded after its
% recorded terminal failure; no simulation has been started in that model.
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
branches = find_system(parent,'SearchDepth',1,'regexp','on', ...
    'Name','^Linear Wave Excitation Force Variant Subsystem[0-9]*$');
assert(numel(branches) == 3,'Expected stock dispatch to duplicate each imported wave');
sumPorts = zeros(1,3);
for i = 1:3
    imported = find_system(branches{i},'SearchDepth',1,'LookUnderMasks','all', ...
        'FollowLinks','on','VariantControl','sv_udfWaves');
    assert(numel(imported) == 1 && eval(sv_udfWaves.Condition));
    ports = get_param(branches{i},'PortHandles');
    line = get_param(ports.Outport(1),'Line');
    destination = get_param(line,'DstPortHandle');
    assert(numel(destination) == 1 && strcmp(get_param(destination,'Parent'),[parent '/Add']));
    sumPorts(i) = get_param(destination,'PortNumber');
end
signs = get_param([parent '/Add'],'Inputs');
assert(numel(unique(sumPorts)) == 3 && all(signs(sumPorts) == '+'));
histories = find_system(parent,'LookUnderMasks','all','FollowLinks','on','BlockType','FromWorkspace');
expressions = {};
for i = 1:numel(histories)
    value = get_param(histories{i},'VariableName');
    if contains(value,'userDefinedFe'), expressions{end+1} = value; end
end
assert(numel(expressions) == 3 && numel(unique(expressions)) == 1);
prepared = body(1).hydroForce.hf1.userDefinedFe;
assert(any(abs(prepared) > 1e-8,'all'));
result = struct('complete',true,'actual_branch_count',numel(branches), ...
    'identical_total_force_expressions',{expressions},'sum_input_signs',signs, ...
    'positive_sum_input_ports',sumPorts,'imported_variant_condition',sv_udfWaves.Condition, ...
    'expected_multiplier',1,'actual_multiplier',3,'body_simulated',false, ...
    'prepared_initial_force',prepared(1,:),'prepared_force_peak',max(abs(prepared),[],1), ...
    'result','confirmed three branches dispatch the same precomputed total excitation');
fid = fopen(fullfile(caseDir,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fclose(fid);
save_system('minimum_reference');
close_system('minimum_reference',0);
diary off
cd(root);
fprintf('Confirmed native mask applies imported total three times.\n');

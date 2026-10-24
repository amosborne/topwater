root = 'C:/Users/alex/topwater';
parent = 'minimum_reference/Boat/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation';
assert(strcmp(get_param('minimum_reference','SimulationStatus'),'stopped'));
branches = find_system(parent,'SearchDepth',1,'LookUnderMasks','all','FollowLinks','on', ...
    'regexp','on','Name','^Linear Wave Excitation Force Variant Subsystem[0-9]*$');
sumPorts = zeros(1,numel(branches));
conditions = cell(1,numel(branches));
for i = 1:numel(branches)
    ports = get_param(branches{i},'PortHandles');
    destination = get_param(get_param(ports.Outport(1),'Line'),'DstPortHandle');
    assert(numel(destination) == 1 && strcmp(get_param(destination,'Parent'),[parent '/Add']));
    sumPorts(i) = get_param(destination,'PortNumber');
    imported = find_system(branches{i},'SearchDepth',1,'LookUnderMasks','all','FollowLinks','on', ...
        'MatchFilter',@Simulink.match.allVariants,'VariantControl','sv_udfWaves');
    assert(numel(imported) == 1);
    conditions{i} = get_param(imported{1},'VariantControl');
end
signs = get_param([parent '/Add'],'Inputs');
assert(all(signs(sumPorts) == '+'));
active = find_system(parent,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'BlockType','FromWorkspace');
activeImported = 0;
for i = 1:numel(active)
    if contains(get_param(active{i},'VariableName'),'userDefinedFe'), activeImported = activeImported+1; end
end
assert(activeImported == numel(branches) && eval(sv_udfWaves.Condition));
histories = find_system(parent,'LookUnderMasks','all','FollowLinks','on','BlockType','FromWorkspace');
expressions = {};
for i = 1:numel(histories)
    value = get_param(histories{i},'VariableName');
    if contains(value,'userDefinedFe'), expressions{end+1} = value; end
end
result = struct('actual_branches',{branches},'actual_branch_count',numel(branches), ...
    'imported_history_expressions',{expressions},'wavegroup_parameter',get_param(parent,'waveGroup'), ...
    'sum_inputs',signs,'positive_sum_ports',sumPorts,'active_imported_histories',activeImported, ...
    'imported_variant_control',{conditions},'imported_variant_condition',sv_udfWaves.Condition, ...
    'condition_value',eval(sv_udfWaves.Condition),'prepared_wave_objects',numel(waves), ...
    'body_simulated',false);
fid = fopen(fullfile(root,'artifacts/reference-imported-dispatch-stock-attempt4/wiring.json'),'w');
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fclose(fid);
fprintf('%s\n',jsonencode(result));

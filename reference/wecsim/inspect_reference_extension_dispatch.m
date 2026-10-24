function evidence = inspect_reference_extension_dispatch(model)
% Read-only actual native imported-total path; no forcing/model replacement.
set_param(model,'SimulationCommand','update');
parent = [model '/Boat/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation'];
branches = find_system(parent,'SearchDepth',1,'LookUnderMasks','all','FollowLinks','on', ...
    'regexp','on','Name','^Linear Wave Excitation Force Variant Subsystem[0-9]*$');
assert(numel(branches) == 1);
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
assert(numel(selected) == 1);
evidence = struct('active_positive_imported_branches',1,'block',selected{1}, ...
    'actual_source',get_param(selected{1},'VariableName'),'forcing_override',false);
end

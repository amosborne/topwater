function cic = log_mild_packet_cic(model)
root = [model '/Boat/Hydrodynamic Body/Wave Radiation Forces Calculation'];
choices = find_system(root,'LookUnderMasks','all','FollowLinks','on', ...
    'MatchFilter',@Simulink.match.activeVariants,'Name','Convolution Integral Calculation');
selected = {};
for j = 1:numel(choices)
    if getSimulinkBlockHandle([choices{j} '/MATLAB Function1']) > 0
        selected{end+1} = choices{j};
    end
end
assert(numel(selected) == 1);
cic = selected{1};
if strcmp(get_param(cic,'LinkStatus'),'resolved')
    set_param(cic,'LinkStatus','inactive');
end
assert(any(strcmp(get_param(cic,'LinkStatus'),{'none','inactive'})));
sources = {'Rate Transition/1','Rate Transition1/1','MATLAB Function/1'};
variables = {'mildCicTimes','mildCicForces','mildCicRadiation'};
for j = 1:3
    name = ['Mild log ' num2str(j)];
    add_block('simulink/Sinks/To Workspace',[cic '/' name], ...
        'VariableName',variables{j},'SaveFormat','Structure With Time',position=[1000 50 1130 90]);
    add_line(cic,sources{j},[name '/1'],autorouting='on');
end
end

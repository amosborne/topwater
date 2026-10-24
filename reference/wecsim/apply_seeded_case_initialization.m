function evidence = apply_seeded_case_initialization(model,seed,n,hf,cicTime)
free = [model '/Free6DOF'];
if strcmp(get_param(free,'LinkStatus'),'resolved')
    set_param(free,'LinkStatus','inactive');
end
joint = [free '/CONSTRAINT (Bushing Joint)'];
jointEvidence = set_initial_joint_state(joint,n.cg,seed.q0(1:3), ...
    seed.velocity_world0(1:3),seed.velocity_world0(4:6));
correction = cic_seed_table(hf.irkb,cicTime,seed.history_velocity_world);
assignin('base','cicSeedCorrection',correction);
assignin('base','cicSeedInitialFrad',correction(3:-1:1,:));
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
assert(any(strcmp(get_param(cic,'LinkStatus'),{'none','inactive'})),'CIC must be an editable case copy');
delay = [model '/Boat/Hydrodynamic Body/' sprintf('Transport\nDelay')];
before = struct('delay',get_param(delay,'DelayTime'),'initial',get_param(delay,'InitialOutput'));
cicEvidence = apply_cic_history_adapter(cic);
log_signal(cic,'Rate Transition/1','Seed history times','seededCicTimes');
log_signal(cic,'Rate Transition1/1','Seed history forces','seededCicForces');
log_signal(cic,'MATLAB Function/1','Seed radiation','seededCicRadiation');
after = struct('delay',get_param(delay,'DelayTime'),'initial',get_param(delay,'InitialOutput'));
assert(isequal(before,after) && str2double(after.delay) == 1e-7 && str2double(after.initial) == 0);
evidence = struct('joint',jointEvidence,'cic',cicEvidence,'Ar_delay_unchanged',after, ...
    'seed_source',seed.source_file,'scope','initialstates and missingpast radiation only; no prescribedpost0motion');
end

function log_signal(system,source,name,variable)
add_block('simulink/Sinks/To Workspace',[system '/' name], ...
    'VariableName',variable,'SaveFormat','Structure With Time',position=[1000 50 1130 90]);
add_line(system,source,[name '/1'],autorouting='on');
end

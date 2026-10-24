function evidence = apply_cic_history_adapter(cicPath)
% Modify this model copy only; native function and pchip/existing delay stay unchanged.
native = [cicPath '/MATLAB Function1'];
timeRT = [cicPath '/Rate Transition'];
forceRT = [cicPath '/Rate Transition1'];
assert(getSimulinkBlockHandle(native) > 0 && getSimulinkBlockHandle(forceRT) > 0);
assert(strcmp(get_param(timeRT,'InitialCondition'),'simu.cicDt*[-1 -2 -3]'));
ports = get_param(native,'PortHandles');
forcePorts = get_param(forceRT,'PortHandles');
line = get_param(forcePorts.Inport(1),'Line');
assert(get_param(line,'SrcPortHandle') == ports.Outport(2));
block = [cicPath '/Seeded past force rows'];
assert(getSimulinkBlockHandle(block) == -1);
add_block('simulink/User-Defined Functions/MATLAB Function',block,position=[600 260 760 340]);
r = sfroot;
chart = r.find('-isa','Stateflow.EMChart','Path',block);
assert(numel(chart) == 1);
chart.Script = sprintf(['function combined = seeded_rows(timeExtrap,FradExtrap,correction,dt)\n' ...
    'combined = cic_seed_lookup(timeExtrap,FradExtrap,correction,dt);\nend\n']);
add_block('simulink/Sources/Constant',[cicPath '/Past correction table'], ...
    'Value','cicSeedCorrection',position=[400 290 520 320]);
add_block('simulink/Sources/Constant',[cicPath '/Seed dt'], ...
    'Value','simu.cicDt',position=[400 340 520 370]);
delete_line(line);
add_line(cicPath,'MATLAB Function1/1','Seeded past force rows/1',autorouting='on');
add_line(cicPath,'MATLAB Function1/2','Seeded past force rows/2',autorouting='on');
add_line(cicPath,'Past correction table/1','Seeded past force rows/3',autorouting='on');
add_line(cicPath,'Seed dt/1','Seeded past force rows/4',autorouting='on');
add_line(cicPath,'Seeded past force rows/1','Rate Transition1/1',autorouting='on');
set_param(forceRT,'InitialCondition','cicSeedInitialFrad');
evidence = struct('path',cicPath,'insertion','before forceRateTransition and nativePCHIP', ...
    'initial_force_rows','fullpast[-dt,-2dt,-3dt] pairedwithnativetimeRateTransition', ...
    'native_CIC_function','unchanged','native_extrapolation','unchangedPCHIP', ...
    'residual_added_mass_delay','unchanged0/default/extrap');
end

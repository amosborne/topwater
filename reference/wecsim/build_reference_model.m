function build_reference_model(caseDir,templateOnly)
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
name = 'minimum_reference';
if templateOnly
    assert(~bdIsLoaded(name));
    new_system(name);
    load_system('WECSim_Lib_Frames');
    add_block('WECSim_Lib_Frames/Global Reference Frame',[name '/World'],position=[40 100 180 160]);
    save_system(name,fullfile(caseDir,[name '.slx']));
    close_system(name,0);
    return
end
assert(bdIsLoaded(name));
existing = find_system(name,'SearchDepth',1,'Type','Block');
assert(isequal(existing,{[name '/World']}),'Template must contain only native global frame');
load_system('WECSim_Lib_Frames');
load_system('WECSim_Lib_Constraints');
load_system('WECSim_Lib_Body_Elements');
add_block('WECSim_Lib_Constraints/Floating (6DOF)',[name '/Free6DOF'], ...
    position=[250 110 320 150],constraint='constraint(1)');
add_block('WECSim_Lib_Body_Elements/Rigid Body',[name '/Boat'], ...
    position=[400 100 540 160],body='body(1)');
world = get_param([name '/World'],'PortHandles');
joint = get_param([name '/Free6DOF'],'PortHandles');
boat = get_param([name '/Boat'],'PortHandles');
add_line(name,world.RConn(1),joint.LConn(1),autorouting='on');
add_line(name,joint.RConn(1),boat.RConn(1),autorouting='on');
set_param(name,'Solver','ode4','FixedStep','simu.dt','StopTime','simu.endTime', ...
    'ReturnWorkspaceOutputs','off','SimMechanicsOpenEditorOnUpdate','off');
save_system(name,fullfile(caseDir,[name '.slx']));
close_system(name,0);
fprintf('Created freely floating six-DOF model from native WEC-Sim library blocks.\n');
end

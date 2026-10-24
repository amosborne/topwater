function evidence = set_initial_joint_state(joint,meanCG,initialDisplacement,vCG,omega)
c = meanCG(:)+initialDisplacement(:);
vPivot = vCG(:)-cross(omega(:),c);
primitives = {'Px','Py','Pz','Rx','Ry','Rz'};
values = [vPivot(:);omega(:)];
for j = 1:6
    prefix = primitives{j};
    if j <= 3
        positionUnit = 'm'; velocityUnit = 'm/s';
    else
        positionUnit = 'rad'; velocityUnit = 'rad/s';
    end
    set_param(joint,[prefix 'PositionTargetSpecify'],'on', ...
        [prefix 'PositionTargetPriority'],'High',[prefix 'PositionTargetValue'],'0', ...
        [prefix 'PositionTargetValueUnits'],positionUnit, ...
        [prefix 'VelocityTargetSpecify'],'on',[prefix 'VelocityTargetPriority'],'High', ...
        [prefix 'VelocityTargetValue'],mat2str(values(j),17),[prefix 'VelocityTargetValueUnits'],velocityUnit);
    assert(strcmp(get_param(joint,[prefix 'MotionActuationMode']),'ComputedMotion'));
    assert(strcmp(get_param(joint,[prefix 'TorqueActuationMode']),'NoTorque'));
end
evidence = struct('desired_CG_velocity_world',vCG(:)', ...
    'desired_omega_world',omega(:)','pivot_arm_world_at_joint0',c', ...
    'joint_velocity_targets_m_s_rad_s',values(:)', ...
    'scope','initial targets only; post0 all6computedmotion/free');
end

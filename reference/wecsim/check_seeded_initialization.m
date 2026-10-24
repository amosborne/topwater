function check_seeded_initialization()
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
maxNumCompThreads(1);
out = fullfile(root,'artifacts/reference/generated-bank-consumer-20261024/seed-initialization-fixture-attempt3');
assert(~isfolder(out));
mkdir(out);
cicResult = check_cic_fixture(out);
kinematics = check_kinematics_fixture(out);
jointResult = check_joint_fixture(out,root);
r = struct('complete',true,'cic',cicResult,'joint',jointResult,'kinematics',kinematics, ...
    'scope','ONE interface fixture; no canonicalbodytrajectory/BEM/customODE or Ardelaychange');
fid = fopen(fullfile(out,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(r,PrettyPrint=true));
fclose(fid);
fprintf('%s\n',jsonencode(r));
end

function result = check_cic_fixture(out)
model = 'cic_history_fixture';
assert(~bdIsLoaded(model));
dt = .02;
lags = (0:2)*dt;
N = 2;
K = zeros(3,6,6);
K(1,:,:) = reshape(1:36,6,6);
K(2,:,:) = -reshape(37:72,6,6)/3;
K(3,:,:) = reshape(73:108,6,6)/7;
history = zeros(N+4,6);
for j = -N-3:0
    history(j+N+4,:) = past_velocity(j,dt);
end
correction = cic_seed_table(K,lags,history);
simu = struct('cicDt',dt,'cicTime',lags);
post = zeros(5,6);
for n = 0:4
    post(n+1,:) = post_velocity(n);
end
assignin('base','simu',simu);
assignin('base','fixtureKernel',K);
assignin('base','fixtureVelocity',[(0:4)'*dt post]);
assignin('base','cicSeedCorrection',correction);
assignin('base','cicSeedInitialFrad',correction(3:-1:1,:));
load_system('WECSim_Lib_Body_Elements');
new_system(model);
source = ['WECSim_Lib_Body_Elements/Rigid Body/Hydrodynamic Body/Wave Radiation Forces Calculation/' ...
    'SS CI and Constant-Damping-CoeVariant Subsystem/Convolution Integral Calculation/' ...
    'Convolution Variant Subsystem/Convolution Integral Calculation'];
cic = [model '/CIC'];
add_block(source,cic,position=[250 50 550 250]);
set_param(cic,'LinkStatus','inactive');
assert(strcmp(get_param([cic '/MATLAB Function1'],'SystemSampleTime'),'simu.cicDt'));
nativePorts = get_param([cic '/MATLAB Function1'],'PortHandles');
busPorts = get_param([cic '/' sprintf('Bus\nSelector')],'PortHandles');
kernelLine = get_param(nativePorts.Inport(2),'Line');
assert(get_param(kernelLine,'SrcPortHandle') == busPorts.Outport(1));
assert(isequal(get_param(kernelLine,'DstPortHandle'),nativePorts.Inport(2)));
delete_line(kernelLine);
delete_block([cic '/' sprintf('Bus\nSelector')]);
delete_block([cic '/From']);
add_block('simulink/Sources/Constant',[cic '/Fixture kernel'], ...
    'Value','fixtureKernel',position=[100 90 240 120]);
add_line(cic,'Fixture kernel/1','MATLAB Function1/2',autorouting='on');
apply_cic_history_adapter(cic);
add_block('simulink/Sources/From Workspace',[model '/Velocity'], ...
    'VariableName','fixtureVelocity','SampleTime',num2str(dt),'Interpolate','off', ...
    'OutputAfterFinalValue','Holding final value',position=[30 40 170 80]);
add_line(model,'Velocity/1','CIC/1');
add_block('simulink/Sources/Clock',[model '/Clock'],position=[30 150 60 180]);
add_block('simulink/Signal Routing/Goto',[model '/Clock tag'], ...
    'GotoTag','Clock','TagVisibility','global',position=[100 150 170 180]);
add_line(model,'Clock/1','Clock tag/1');
log_signal(cic,'Rate Transition/1','Logged time rows','fixtureTimeRows');
log_signal(cic,'Rate Transition1/1','Logged force rows','fixtureForceRows');
log_signal(cic,'MATLAB Function/1','Logged radiation','fixtureRadiation');
set_param(model,'Solver','FixedStepDiscrete','FixedStep',num2str(dt/4), ...
    'StopTime',num2str(4*dt),'ReturnWorkspaceOutputs','on');
save_system(model,fullfile(out,[model '.slx']));
clear convolutionIntegralInterp
s = sim(model,'ReturnWorkspaceOutputs','on');
timeLog = s.get('fixtureTimeRows');
forceLog = s.get('fixtureForceRows');
radLog = s.get('fixtureRadiation');
save(fullfile(out,'cic-native-logs.mat'),'timeLog','forceLog','radLog','correction','K','history','-v7');
fastTimes = (0:16)'*(dt/4);
assert(numel(timeLog.time) == 17 && max(abs(timeLog.time(:)-fastTimes)) < 1e-12);
assert(numel(forceLog.time) == 17 && max(abs(forceLog.time(:)-fastTimes)) < 1e-12);
assert(numel(radLog.time) == 17 && max(abs(radLog.time(:)-fastTimes)) < 1e-12);
maximum = 0;
minorMaximum = 0;
actualPchipMaximum = 0;
zeroTail = 0;
for sample = 1:numel(timeLog.time)
    times = sample_values(timeLog,sample,[3 1]);
    actual = sample_values(forceLog,sample,[3 6]);
    expected = zeros(3,6);
    for row = 1:3
        n = round(times(row)/dt);
        expected(row,:) = direct_full_force(n,K,dt,N);
    end
    maximum = max(maximum,max(abs(actual-expected),[],'all'));
    query = timeLog.time(sample);
    radIndex = find(abs(radLog.time-query) < 1e-12,1);
    assert(~isempty(radIndex));
    actualRadiation = sample_values(radLog,radIndex,[6 1]);
    expectedRadiation = interp1(times,expected,query,'pchip','extrap')';
    actualPchipMaximum = max(actualPchipMaximum,max(abs(actualRadiation-expectedRadiation)));
    for fraction = [0 .25 .5 .75]
        query = timeLog.time(sample)+fraction*dt;
        a = extrapFrad(times,actual,query);
        b = interp1(times,expected,query,'pchip','extrap')';
        minorMaximum = max(minorMaximum,max(abs(a-b)));
    end
    if min(round(times/dt)) >= N
        zeroTail = max(zeroTail,max(abs(cic_seed_lookup(times,zeros(3,6),correction,dt)),[],'all'));
    end
end
assert(all(cic_seed_lookup(dt*[N+2;N+1;N],zeros(3,6),correction,dt) == 0,'all'));
assert(maximum < 1e-10 && minorMaximum < 1e-10 && actualPchipMaximum < 1e-10 && zeroTail == 0);
result = struct('main_force_row_error',maximum,'native_pchip_minor_query_error',minorMaximum, ...
    'actual_native_PCHIP_faststep_error',actualPchipMaximum,'fixture_base_dt',dt/4,'CIC_dt',dt, ...
    'zero_correction_after_N_plus_2',zeroTail,'main_samples',numel(timeLog.time), ...
    'kernel','3lag oddlength/asymmetric6x6,dt.02','minor_queries','actualnativePCHIPoutput atdt.005 plus standalone quarterwindow checks');
close_system(model,0);
end

function result = check_kinematics_fixture(out)
Q = [.1+.03i .2-.04i -.03+.02i .12+.02i -.08+.03i .04-.01i; ...
    -.02+.01i .03+.02i .05-.01i .01-.03i .02+.01i -.01+.02i];
w = [.7;1.3];
epsilon = 1e-6;
[q,v,a] = harmonic_seed_state(Q,w,[-epsilon;0;epsilon]);
Rm = fixture_rotation(q(1,4:6));
R0 = fixture_rotation(q(2,4:6));
Rp = fixture_rotation(q(3,4:6));
W = ((Rp-Rm)/(2*epsilon))*R0';
omega = [W(3,2) W(1,3) W(2,1)];
alpha = (v(3,4:6)-v(1,4:6))/(2*epsilon);
omegaError = max(abs(omega-v(2,4:6)));
alphaError = max(abs(alpha-a(2,4:6)));
assert(omegaError < 1e-7 && alphaError < 1e-7);
Q_kernel40 = Q;
omega_rad_s = w;
save(fullfile(out,'synthetic-seed.mat'),'Q_kernel40','omega_rad_s','-v7');
seed = make_harmonic_seed(fullfile(out,'synthetic-seed.mat'),0:.02:.04);
assert(max(abs(seed.q0-q(2,:))) < 1e-13);
assert(max(abs(seed.velocity_world0-v(2,:))) < 1e-13);
assert(isequal(size(seed.history_velocity_world),[6 6]));
result = struct('finite_rotation_derivative_worldomega_error',omegaError, ...
    'worldalpha_finite_derivative_error',alphaError,'scope','algebra/kinematic finite difference, not dynamics');
end

function R = fixture_rotation(angles)
r = angles(1); p = angles(2); y = angles(3);
Rx = [1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)];
Ry = [cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Rz = [cos(y) -sin(y) 0;sin(y) cos(y) 0;0 0 1];
R = Rz*Ry*Rx;
end

function v = past_velocity(n,dt)
v = (1:6)*(1+.3*n*dt)+(1:6).^2*(n*dt)^2;
end

function v = post_velocity(n)
v = (1:6)+n*.07*[3 1 4 2 6 5]+n^2*.01*[2 -1 3 -2 1 -3];
end

function F = direct_full_force(n,K,dt,N)
F = zeros(1,6);
for lag = 0:N
    if n-lag < 0
        v = past_velocity(n-lag,dt);
    else
        v = post_velocity(n-lag);
    end
    weight = dt;
    if lag == 0 || lag == N, weight = dt/2; end
    for i = 1:6
        for j = 1:6
            F(i) = F(i)+weight*K(lag+1,i,j)*v(j);
        end
    end
end
end

function log_signal(system,source,name,variable)
add_block('simulink/Sinks/To Workspace',[system '/' name], ...
    'VariableName',variable,'SaveFormat','Structure With Time',position=[1000 50 1130 90]);
add_line(system,source,[name '/1'],autorouting='on');
end

function value = sample_values(log,k,shape)
data = log.signals.values;
nt = numel(log.time);
if size(data,1) == nt
    value = reshape(data(k,:,:),shape);
else
    assert(size(data,3) == nt);
    value = reshape(data(:,:,k),shape);
end
end

function result = check_joint_fixture(out,root)
model = 'joint_state_fixture';
assert(~bdIsLoaded(model));
simu = simulationClass();
simu.dt = .02;
simu.gravity = 9.81;
assignin('base','simu',simu);
% These native World variants are otherwise supplied by initializeWecSim.
assignin('base','visON',0);
assignin('base','sv_visualizationON',Simulink.Variant('visON==1'));
assignin('base','sv_visualizationOFF',Simulink.Variant('visON==0'));
meanCG = [.3 -.2 .1];
displacement = [.2 .1 -.05];
angles = [.12 -.08 .04];
vCG = [.17 -.23 .05];
omega = [.07 -.04 .03];
r = angles(1); p = angles(2); y = angles(3);
Rx = [1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)];
Ry = [cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Rz = [cos(y) -sin(y) 0;sin(y) cos(y) 0;0 0 1];
R = Rz*Ry*Rx;
q = [cos(y/2)*cos(p/2)*cos(r/2)+sin(y/2)*sin(p/2)*sin(r/2), ...
    cos(y/2)*cos(p/2)*sin(r/2)-sin(y/2)*sin(p/2)*cos(r/2), ...
    cos(y/2)*sin(p/2)*cos(r/2)+sin(y/2)*cos(p/2)*sin(r/2), ...
    sin(y/2)*cos(p/2)*cos(r/2)-cos(y/2)*sin(p/2)*sin(r/2)];
axis = q(2:4)/norm(q(2:4));
angle = 2*atan2(norm(q(2:4)),q(1));
new_system(model);
load_system('WECSim_Lib_Frames');
source = 'WECSim_Lib_Body_Elements/Rigid Body/Hydrodynamic Body/Structure';
add_block('WECSim_Lib_Frames/Global Reference Frame',[model '/World'],position=[30 50 140 100]);
add_block('sm_lib/Joints/Bushing Joint',[model '/Joint'],position=[200 50 260 130]);
for primitive = {'Px','Py','Pz','Rx','Ry','Rz'}
    pre = primitive{1};
    set_param([model '/Joint'],[pre 'MotionActuationMode'],'ComputedMotion', ...
        [pre 'TorqueActuationMode'],'NoTorque',[pre 'SpringStiffness'],'0', ...
        [pre 'DampingCoefficient'],'0',[pre 'LowerLimitSpecify'],'off',[pre 'UpperLimitSpecify'],'off');
end
targets = set_initial_joint_state([model '/Joint'],meanCG,displacement,vCG,omega);
add_block([source '/' sprintf('Rigid\nTransform\n(Setup CG)')],[model '/Mean CG'], ...
    'TranslationCartesianOffset',mat2str(meanCG,17),position=[310 60 360 110]);
add_block([source '/' sprintf('Rigid\nTransform\n(Setup Initial Disp)')],[model '/Initial pose'], ...
    'TranslationCartesianOffset',mat2str(displacement,17),'RotationArbitraryAxis',mat2str(axis,17), ...
    'RotationAngle',mat2str(angle,17),position=[410 60 460 110]);
add_block([source '/Body Properties'],[model '/Solid'], ...
    'ExtGeomFileName',fullfile(root,'artifacts/reference/generated-bank-native-20261024/attempt1/boat.stl'), ...
    'Mass','1','MomentsOfInertia','[1 1 1]','ProductsOfInertia','[0 0 0]', ...
    'GraphicDiffuseColor','[0.8 0.8 0.8]','GraphicOpacity','1',position=[520 50 600 130]);
add_block([source '/Motion Sensor'],[model '/Sensor'],position=[450 200 560 280]);
world = get_param([model '/World'],'PortHandles');
joint = get_param([model '/Joint'],'PortHandles');
cg = get_param([model '/Mean CG'],'PortHandles');
init = get_param([model '/Initial pose'],'PortHandles');
solid = get_param([model '/Solid'],'PortHandles');
sensor = get_param([model '/Sensor'],'PortHandles');
add_line(model,world.RConn(1),joint.LConn(1),autorouting='on');
add_line(model,joint.RConn(1),cg.LConn(1),autorouting='on');
add_line(model,cg.RConn(1),init.LConn(1),autorouting='on');
add_line(model,init.RConn(1),solid.RConn(1),autorouting='on');
add_line(model,world.RConn(1),sensor.LConn(1),autorouting='on');
add_line(model,init.RConn(1),sensor.RConn(1),autorouting='on');
add_block('simulink/Signal Routing/Bus Selector',[model '/Response'], ...
    'OutputSignals','disp,vel,acc',position=[650 180 690 300]);
add_line(model,'Sensor/1','Response/1');
log_signal(model,'Response/1','Pose log','fixturePose');
log_signal(model,'Response/2','Velocity log','fixtureVelocityReadback');
set_param(model,'Solver','ode4','FixedStep','.02','StopTime','0', ...
    'ReturnWorkspaceOutputs','on','SimMechanicsOpenEditorOnUpdate','off');
save_system(model,fullfile(out,[model '.slx']));
s = sim(model,'ReturnWorkspaceOutputs','on');
pose = s.get('fixturePose');
velocity = s.get('fixtureVelocityReadback');
save(fullfile(out,'joint-native-logs.mat'),'pose','velocity','targets','R','-v7');
actualPose = reshape(pose.signals.values(1,:),1,6);
actualVelocity = reshape(velocity.signals.values(1,:),1,6);
expectedPose = [meanCG+displacement angles];
expectedVelocity = [vCG omega];
assert(max(abs(actualPose-expectedPose)) < 1e-9);
assert(max(abs(actualVelocity-expectedVelocity)) < 1e-9);
result = struct('actual_initial_pose',actualPose,'expected_pose',expectedPose, ...
    'actual_CG_velocity_worldomega',actualVelocity,'expected_velocity',expectedVelocity, ...
    'targets',targets,'simulation_duration',0,'scope','nativemechanicalinitialization only; nohydrotrajectory');
close_system(model,0);
end

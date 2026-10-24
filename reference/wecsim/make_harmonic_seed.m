function seed = make_harmonic_seed(seedFile,cicTime)
s = load(seedFile,'Q_kernel40','omega_rad_s');
dt = cicTime(2)-cicTime(1);
N = numel(cicTime)-1;
historyTime = (-N-3:0)'*dt;
[q,v,a] = harmonic_seed_state(s.Q_kernel40,s.omega_rad_s,historyTime);
angles = q(end,4:6);
r = angles(1); p = angles(2); y = angles(3);
quat = [cos(y/2)*cos(p/2)*cos(r/2)+sin(y/2)*sin(p/2)*sin(r/2), ...
    cos(y/2)*cos(p/2)*sin(r/2)-sin(y/2)*sin(p/2)*cos(r/2), ...
    cos(y/2)*sin(p/2)*cos(r/2)+sin(y/2)*cos(p/2)*sin(r/2), ...
    sin(y/2)*cos(p/2)*cos(r/2)-cos(y/2)*sin(p/2)*sin(r/2)];
if norm(quat(2:4)) == 0
    axis = [0 1 0];
else
    axis = quat(2:4)/norm(quat(2:4));
end
angle = 2*atan2(norm(quat(2:4)),quat(1));
seed = struct('q0',q(end,:),'velocity_world0',v(end,:),'harmonic_acceleration_world0',a(end,:), ...
    'initial_axis',axis,'initial_angle',angle,'history_time',historyTime, ...
    'history_q',q,'history_velocity_world',v,'history_acceleration_world',a, ...
    'source_file',seedFile,'scope','meanlinear particular seed/finiteworldkinematics, not exactnative stationary solution', ...
    'Ar_delay','unchanged defaultzero/extrapolation; no acceleration state correction');
end

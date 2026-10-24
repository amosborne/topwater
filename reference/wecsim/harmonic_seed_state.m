function [q,velocity,acceleration] = harmonic_seed_state(Q,omega,time)
phase = exp(-1i*time(:)*omega(:)');
q = real(phase*Q);
dq = real(phase*(-1i*omega(:).*Q));
ddq = real(phase*(-omega(:).^2.*Q));
velocity = dq;
acceleration = ddq;
for j = 1:numel(time)
    pitch = q(j,5); yaw = q(j,6);
    dp = dq(j,5); dy = dq(j,6);
    E = [cos(yaw)*cos(pitch) -sin(yaw) 0; ...
        sin(yaw)*cos(pitch) cos(yaw) 0;-sin(pitch) 0 1];
    dE = [-sin(yaw)*dy*cos(pitch)-cos(yaw)*sin(pitch)*dp -cos(yaw)*dy 0; ...
        cos(yaw)*dy*cos(pitch)-sin(yaw)*sin(pitch)*dp -sin(yaw)*dy 0;-cos(pitch)*dp 0 0];
    velocity(j,4:6) = (E*dq(j,4:6)')';
    acceleration(j,4:6) = (E*ddq(j,4:6)'+dE*dq(j,4:6)')';
end
end

function prepare_jan90_case_inputs(caseDir)
root = 'C:/Users/alex/topwater';
native = fullfile(root,'artifacts/reference/generated-bank-native-20261024/attempt1');
scout = fullfile(root,'artifacts/reference/radiation-heading-scout-20261024/attempt1');
n = load(fullfile(native,'native.mat'));
pair = load(fullfile(scout,'jan05-offset90-pairs.mat'));
screen = jsondecode(fileread(fullfile(root,'artifacts/reference/generated-bank-consumer-20261024/jan90-attempt1/summary.json')));
assert(screen.complete);
assert(isequal(pair.component_index(:),(1:384)'));
assert(isequal(pair.excitation_complex,pair.fk_complex+pair.diffraction_complex));
w = pair.omega_rad_s(:);
a = pair.amplitude_m(:);
theta = pair.heading_native_rad(:);
phi = pair.phase_local_rad(:);
k = w.^2/n.g;
phaseCG = phi+k.*(n.cg(1)*cos(theta)+n.cg(2)*sin(theta));
eta0 = sum(a.*cos(phaseCG));
sx = -sum(a.*k.*cos(theta).*sin(phaseCG));
sy = -sum(a.*k.*sin(theta).*sin(phaseCG));
stored = screen.records.incident_fields.meanCG_t0_eta_slopes_Hessian_etat;
assert(max(abs([eta0 sx sy]-stored(1:3)')) < 1e-12);
roll0 = atan(sy/sqrt(1+sx^2));
pitch0 = -atan(sx);
initial_translation = [0 0 eta0];
initial_euler = [roll0 pitch0 0];
qw = cos(pitch0/2)*cos(roll0/2);
qxyz = [cos(pitch0/2)*sin(roll0/2),sin(pitch0/2)*cos(roll0/2),-sin(pitch0/2)*sin(roll0/2)];
initial_angle = 2*atan2(norm(qxyz),qw);
if norm(qxyz) == 0
    initial_axis = [0 1 0];
else
    initial_axis = qxyz/norm(qxyz);
end
Rx = [1 0 0;0 cos(roll0) -sin(roll0);0 sin(roll0) cos(roll0)];
Ry = [cos(pitch0) 0 sin(pitch0);0 1 0;-sin(pitch0) 0 cos(pitch0)];
initial_rotation = Ry*Rx;
u = initial_axis(:);
crossU = [0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0];
axis_rotation = eye(3)+sin(initial_angle)*crossU+(1-cos(initial_angle))*(crossU*crossU);
assert(max(abs(axis_rotation-initial_rotation),[],'all') < 1e-13);
time = (0:3000)'*.02;
force = real(exp(-1i*time*w')*(a.*exp(1i*phi).*pair.excitation_complex));
jan90ForceTable = [time force];
for index = [1 20 158 3001]
    expected = zeros(1,6);
    for j = 1:384
        phase = phi(j)-w(j)*time(index);
        expected = expected+a(j)*(real(pair.excitation_complex(j,:))*cos(phase)-imag(pair.excitation_complex(j,:))*sin(phase));
    end
    assert(max(abs(expected-force(index,:))) < 1e-8);
end
historyTime = (-1000:4000)'*.02;
historyEta = real(exp(-1i*historyTime*w')*(a.*exp(1i*phaseCG)));
writematrix([historyTime historyEta],fullfile(caseDir,'wave-full-bank.txt'),delimiter='tab');
save(fullfile(caseDir,'prepared-case.mat'),'n','pair','time','force','jan90ForceTable', ...
    'initial_translation','initial_euler','initial_axis','initial_angle','initial_rotation','eta0','sx','sy','-v7');
r = struct('initial_translation_m',initial_translation,'initial_euler_rad',initial_euler, ...
    'initial_axis',initial_axis,'initial_angle_rad',initial_angle, ...
    'initial_eta_slopes',[eta0 sx sy],'duration_s',60,'dt_s',.02,'radiation_memory_s',40, ...
    'holding_load_linear_estimate', (n.K*[initial_translation initial_euler]'-force(1,:)')', ...
    'holding_load_convention','external generalized holding load = K*q0 - Fexc0; zero v and radiation history; not exact finite-wave load', ...
    'history','body held stationary before0, removed holding constraint at0; incident sea stationary', ...
    'initial_field','native mean-body frequency-domain diffraction idealization; no initialized nonlinear field', ...
    'curvature_caveat','sampled initial footprint plane departure about0.183m; tangentplane pose is not curvedsurface wetset/physicalequilibrium guarantee', ...
    'case_status','provisional selected17panel model trajectory, not trusted mesh/geometry validation');
fid = fopen(fullfile(caseDir,'initialization.json'),'w');
fprintf(fid,'%s\n',jsonencode(r,PrettyPrint=true));
fclose(fid);
end

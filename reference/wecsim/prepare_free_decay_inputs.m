function prepare_free_decay_inputs(caseDir)
root = 'C:/Users/alex/topwater';
n = load(fullfile(root,'artifacts/reference/generated-bank-native-20261024/attempt1/native.mat'));
initial_translation = [0 0 0];
initial_euler = [pi/360 0 0];
initial_axis = [1 0 0];
initial_angle = pi/360;
r = initial_angle;
initial_rotation = [1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)];
time = (0:3000)'*.02;
force = zeros(numel(time),6);
% Keep the reviewed case force-hook workspace name; every value is exactly zero.
jan90ForceTable = [time force];
historyTime = (-1000:4000)'*.02;
writematrix([historyTime zeros(size(historyTime))],fullfile(caseDir,'wave-full-bank.txt'),delimiter='tab');
save(fullfile(caseDir,'prepared-case.mat'),'n','time','force','jan90ForceTable', ...
    'initial_translation','initial_euler','initial_axis','initial_angle','initial_rotation','-v7');
r = struct('control','ONE calm-water free-decay rejection diagnostic; not a canonicalsea', ...
    'initial_translation_m',initial_translation,'initial_euler_rad',initial_euler, ...
    'initial_axis',initial_axis,'initial_angle_rad',initial_angle, ...
    'duration_s',60,'dt_s',.02,'radiation_memory_s',40, ...
    'prior_holding_load_linear_estimate',(n.K*[initial_translation initial_euler]')', ...
    'history','held stationary at0.5degheel in calm water before0; release0; allv/radiationhistory0', ...
    'forcing','all6appliedexcitation exactlyzero; no generatedbank altered', ...
    'operator','same original17panelH5/M/K/Ainf/Ar/CIC/delay; no damping/tuning/modeconstraints', ...
    'interpretation','Sustained/severe growing response rejects bounded small-perturbation use; transient amplification alone is not asymptotic instability proof; decay does not prove passivity or certify forcedcase');
fid = fopen(fullfile(caseDir,'initialization.json'),'w');
fprintf(fid,'%s\n',jsonencode(r,PrettyPrint=true));
fclose(fid);
end

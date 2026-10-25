function prepare_calm_buoy_case_inputs(caseDir)
root = 'C:/Users/alex/topwater';
base = getenv('TOPWATER_BUOY_ROOT');
if isempty(base)
    base = fullfile(root,'artifacts/reference/calm-buoy-20261024');
end
assessment = jsondecode(fileread(fullfile(base,'assessment-attempt1/diagnostic.json')));
assert(assessment.complete && assessment.accepted && strcmp(assessment.field_status,'proved_calm'));
assert(assessment.old_rejection_preserved || ...
    (isfield(assessment,'strict_policy_report_preserved') && assessment.strict_policy_report_preserved), ...
    'Original strict-policy source report must remain preserved');
heading = jsondecode(fileread(fullfile(base,'export-attempt1/case.json')));
n = load(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/native.mat'));
trainDir = getenv('TOPWATER_BUOY_TRAIN');
if isempty(trainDir)
    trainDir = 'train-attempt2';
end
forcing = load(fullfile(base,trainDir,'body-force.mat'));
time = forcing.time(:);
force = forcing.force;
assert(isequal(size(force),[3001 6]) && all(isfinite(force),'all'));
assert(max(abs(time-(0:3000)'*.02)) < 1e-12);
assert(isequal(force(1,:),assessment.force_policy.native_force_N_Nm(:)'));
initial_translation = [0 0 0];
initial_euler = [0 0 0];
initial_axis = [0 1 0];
initial_angle = 0;
initial_rotation = eye(3);
mean_box_yaw_global_rad = heading.mean_box_yaw_global_rad;
global_from_boat_rotation = heading.global_from_boat_rotation;
global_from_boat_translation_m = heading.global_from_boat_translation_m;
% Same native model-world/meanboat axes as original M/K/Ainf/H. Global R is NOT
% applied to model inertia/geometry/force again; output/video applies it once.
jan90ForceTable = [time force];
historyTime = (-1000:4000)'*.02;
writematrix([historyTime zeros(size(historyTime))],fullfile(caseDir,'wave-full-bank.txt'),delimiter='tab');
save(fullfile(caseDir,'prepared-case.mat'),'n','time','force','jan90ForceTable', ...
    'initial_translation','initial_euler','initial_axis','initial_angle','initial_rotation', ...
    'mean_box_yaw_global_rad','global_from_boat_rotation','global_from_boat_translation_m','-v7');
copyfile(fullfile(base,'assessment-attempt1/diagnostic.json'),fullfile(caseDir,'incident-preflight.json'));
copyfile(fullfile(base,'export-attempt1/case.json'),fullfile(caseDir,'heading-case.json'));
end

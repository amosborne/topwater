function prepare_reference_extension_inputs(caseDir)
root = 'C:/Users/alex/topwater';
packetDir = fullfile(root,'artifacts/reference/reference-extension-20261024/incident-attempt1');
startup = jsondecode(fileread(fullfile(packetDir,'diagnostic.json')));
assert(startup.complete && startup.accepted && startup.unchanged_physical_field);
assert(startup.old_wave_prefix_bytes_identical && startup.first501_retained_force_bit_identical);
n = load(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/native.mat'));
forcing = load(fullfile(packetDir,'body-force.mat'));
time = forcing.time(:);
force = forcing.force;
assert(isequal(size(force),[3001 6]) && all(isfinite(force),'all'));
assert(max(abs(time-(0:3000)'*.02)) < 1e-12);
initial_translation = [0 0 0];
initial_euler = [0 0 0];
initial_axis = [0 1 0];
initial_angle = 0;
initial_rotation = eye(3);
save(fullfile(caseDir,'prepared-case.mat'),'n','time','force', ...
    'initial_translation','initial_euler','initial_axis','initial_angle','initial_rotation','-v7');
copyfile(fullfile(packetDir,'diagnostic.json'),fullfile(caseDir,'incident-preflight.json'));
end

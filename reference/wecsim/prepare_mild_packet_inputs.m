function prepare_mild_packet_inputs(caseDir)
root = 'C:/Users/alex/topwater';
packetDir = fullfile(root,'artifacts/reference/mild-packet-20261024/incident-attempt1');
startup = jsondecode(fileread(fullfile(packetDir,'diagnostic.json')));
assert(startup.complete && startup.accepted,'Incident/full-region calm preflight must PASS before any body run');
n = load(fullfile(root,'artifacts/reference/generated-bank-native-20261024/attempt1/native.mat'));
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
jan90ForceTable = [time force];
% Only the existing case hook name is reused; no Jan90 sea/history is consumed.
historyTime = (-1000:4000)'*.02;
writematrix([historyTime zeros(size(historyTime))],fullfile(caseDir,'wave-full-bank.txt'),delimiter='tab');
save(fullfile(caseDir,'prepared-case.mat'),'n','time','force','jan90ForceTable', ...
    'initial_translation','initial_euler','initial_axis','initial_angle','initial_rotation','-v7');
copyfile(fullfile(packetDir,'diagnostic.json'),fullfile(caseDir,'incident-preflight.json'));
end

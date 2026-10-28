function tests = nativeExcitationStartupTest
% Diagnostic scale and unchanged startup limits; no model/body simulation.
tests = functiontests(localfunctions);
end

function setup(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
originalPath = path;
originalDirectory = pwd;
testCase.addTeardown(@() path(originalPath));
addpath(fullfile(root,'src','topwater','data','wecsim'));
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
% Teardowns run last-in first-out: leave the folder before its fixture deletes it.
testCase.addTeardown(@() cd(originalDirectory));
cd(folder.Folder);
w = [1;2];
theta = [0;pi;2*pi];
H = zeros(2,3,6);
for j = 1:6
    H(:,:,j) = repmat(j*(1+.2*w)*(1+.5i),1,3);
end
native = fullfile(folder.Folder,'native.mat');
excitation_omega_rad_s = w;
excitation_headings_rad = theta;
excitation_native = H;
save(native,'excitation_omega_rad_s','excitation_headings_rad','excitation_native');
testCase.TestData.prepared = struct('native_mat_path',native, ...
    'n',struct('M',100*eye(6),'g',9.81),'boat_length_m',2, ...
    'carrier_omega_rad_s',[1.2;1.8],'carrier_directions_rad',[.2;1.1], ...
    'carrier_amplitudes_m',[.3;.4]);
end

function testNativeZeroForceAndDiagnosticScale(testCase)
prepared = testCase.TestData.prepared;
report = native_excitation_startup(prepared,zeros(3,6));
expected = (1:6)*sqrt(1.25)*(.3*(1+.2*1.2)+.4*(1+.2*1.8));
testCase.verifyEqual(report.carrier_scale,expected,'AbsTol',1e-12);
testCase.verifyEqual(report.ramp_time_s,0);
testCase.verifyEqual(report.force_norm_limit,.001*100*9.81,'AbsTol',1e-14);
testCase.verifyEqual(report.moment_norm_limit,.001*100*9.81*2,'AbsTol',1e-14);
end

function testNoncalmNativeForceFailsBeforeBodyRun(testCase)
prepared = testCase.TestData.prepared;
excitation = zeros(3,6);
excitation(1,2) = 1;
failed = false;
try
    native_excitation_startup(prepared,excitation);
catch failure
    failed = contains(failure.message,'unchanged engineering calm limits');
end
testCase.verifyTrue(failed,'Noncalm native force was not rejected');
testCase.verifyTrue(isfile('native-startup.mat'),'Failed startup diagnostics were not retained');
end

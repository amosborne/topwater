function tests = workflowEngineSettingTest
% Current-path regression: initializeWecSim changes this global setting.
tests = functiontests(localfunctions);
end

function setup(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
originalPath = path;
originalSetting = get_param(0,'ErrorIfLoadNewModel');
testCase.addTeardown(@() path(originalPath));
testCase.addTeardown(@() set_param(0,'ErrorIfLoadNewModel',originalSetting));
addpath(fullfile(root,'src','topwater','data','wecsim'));
end

function testRestoreModelVersionSettingAfterSuccess(testCase)
for initial = {'on','off'}
    value = initial{1};
    set_param(0,'ErrorIfLoadNewModel',value);
    originalPath = path;
    [~,model] = fileparts(tempname);
    guard = WorkflowEngineGuard(model,originalPath);
    cleanup = onCleanup(@() release_guard(guard,originalPath));
    guard.prepare();
    set_param(0,'ErrorIfLoadNewModel',opposite(value));
    report = guard.restore();
    testCase.verifyTrue(report.restored,strjoin(report.failures,'; '));
    testCase.verifyEqual(get_param(0,'ErrorIfLoadNewModel'),value);
    release_guard(guard,originalPath);
    clear cleanup guard
end
end

function testRestoreModelVersionSettingAfterFailure(testCase)
for initial = {'on','off'}
    value = initial{1};
    set_param(0,'ErrorIfLoadNewModel',value);
    testCase.verifyError(@() forced_failure(value),'Topwater:ForcedEngineSettingFailure');
    testCase.verifyEqual(get_param(0,'ErrorIfLoadNewModel'),value);
end
end

function forced_failure(initial)
originalPath = path;
[~,model] = fileparts(tempname);
guard = WorkflowEngineGuard(model,originalPath);
cleanup = onCleanup(@() release_guard(guard,originalPath));
guard.prepare();
set_param(0,'ErrorIfLoadNewModel',opposite(initial));
error('Topwater:ForcedEngineSettingFailure','Intentional failure after native setting change');
end

function value = opposite(initial)
if strcmp(initial,'on')
    value = 'off';
else
    value = 'on';
end
end

function release_guard(guard,originalPath)
if isvalid(guard)
    delete(guard);
end
path(originalPath);
end

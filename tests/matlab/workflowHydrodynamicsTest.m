function tests = workflowHydrodynamicsTest
% Finite coefficient guards without MATLAB models, BEM solves or native dynamics.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
originalPath = path;
testCase.addTeardown(@() path(originalPath));
folder = fullfile(root,'src','topwater','data','wecsim');
addpath(folder);
testCase.assertEqual(which('validate_workflow_hydrodynamics'), ...
    fullfile(folder,'validate_workflow_hydrodynamics.m'));
end

function testFiniteControlAllowsInfiniteWaterDepth(testCase)
[h,re,im] = fixture_hydro();
h.h = Inf;
h.unused_numeric_metadata = Inf;
testCase.verifyWarningFree(@() validate_workflow_hydrodynamics(h,'imported',re,im));
testCase.verifyWarningFree(@() validate_workflow_hydrodynamics(h,'prepared'));
end

function testNaNRealExcitationRejectedWithLegacyWitness(testCase)
[h,re,im] = fixture_hydro();
% One nonfinite coefficient appears in both the imported and raw record.
h.ex_re(1,1,end) = NaN;
re(1,1,end) = NaN;
legacyAccepted = max(abs(h.ex_re*h.rho*h.g-re),[],'all') ...
    <= 1e-10*max(1,max(abs(re),[],'all'));
testCase.verifyTrue(legacyAccepted,'The pre-fix phase gate should expose its NaN omission');
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:NonfiniteHydrodynamics','imported.ex_re');
end

function testNaNImaginaryExcitationRejectedWithLegacyWitness(testCase)
[h,re,im] = fixture_hydro();
h.ex_im(1,1,end) = NaN;
im(1,1,end) = NaN;
legacyAccepted = max(abs(h.ex_im*h.rho*h.g+im),[],'all') ...
    <= 1e-10*max(1,max(abs(im),[],'all'));
testCase.verifyTrue(legacyAccepted,'The pre-fix imaginary gate should omit the injected NaN');
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:NonfiniteHydrodynamics','imported.ex_im');
end

function testRawNonfiniteEntryCannotHideInComparison(testCase)
[h,re,im] = fixture_hydro();
re(2,1,end) = NaN;
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:NonfiniteHydrodynamics','netCDF real excitation');
end

function testNonfiniteRadiationCoefficientFailsBeforeTransform(testCase)
[h,re,im] = fixture_hydro();
h.B(2,5,end) = Inf;
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:NonfiniteHydrodynamics','imported.B');
end

function testPreparedKernelNonfiniteFailsBeforeH5(testCase)
h = fixture_hydro();
h.ra_K(2,5,end) = NaN;
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'prepared'), ...
    'Topwater:NonfiniteHydrodynamics','prepared.ra_K');
end

function testDirectNativeOverrideIsCheckedBeforeH5(testCase)
h = fixture_hydro();
h.Ainf(1,5) = NaN;
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'prepared'), ...
    'Topwater:NonfiniteHydrodynamics','prepared.Ainf');
end

function testDiscardedRoundedMetadataIsNotAnImportedGate(testCase)
[h,re,im] = fixture_hydro();
h.cg = NaN(3,1);
h.Khs = NaN(6);
h.Ainf = NaN(6);
testCase.verifyWarningFree(@() validate_workflow_hydrodynamics(h,'imported',re,im));
end

function testCoefficientShapeRejectedBeforeBroadcast(testCase)
[h,re,im] = fixture_hydro();
h.ex_re = h.ex_re(:,1,:);
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:InvalidHydrodynamics','imported.ex_re');
end

function testReferenceGridAndCountsMustBeCompatible(testCase)
[h,re,im] = fixture_hydro();
h.w = [.1 .1];
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:InvalidHydrodynamics','imported.w');
h = fixture_hydro();
h.Nf = 3;
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'imported',re,im), ...
    'Topwater:InvalidHydrodynamics','reference grids');
end

function testKernelShapeMustMatchTimeGrid(testCase)
h = fixture_hydro();
h.ra_t = [0 .02 .04 .06];
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'prepared'), ...
    'Topwater:InvalidHydrodynamics','prepared.ra_K');
end

function testOnlyExcitationKernelAllowsZeroFrequency(testCase)
h = fixture_hydro();
h.ex_w = [0 .5 1];
testCase.verifyWarningFree(@() validate_workflow_hydrodynamics(h,'prepared'));
h.ra_w = [0 1];
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'prepared'), ...
    'Topwater:InvalidHydrodynamics','prepared.ra_w');
end

function testExcitationKernelRejectsNegativeFrequency(testCase)
h = fixture_hydro();
h.ex_w = [-.1 .5 1];
verify_failure(testCase,@() validate_workflow_hydrodynamics(h,'prepared'), ...
    'Topwater:InvalidHydrodynamics','prepared.ex_w');
end

function [h,re,im] = fixture_hydro()
h = struct('rho',1025,'g',9.81,'h',Inf,'w',[.1 1], ...
    'theta',[0 45 90],'Nf',2,'Nh',3,'cg',zeros(3,1),'cb',zeros(3,1), ...
    'Vo',1,'Khs',eye(6),'Ainf',eye(6));
h.T = 2*pi./h.w;
h.A = reshape(1:72,6,6,2);
h.B = 2*h.A;
for name = {'fk','sc','ex'}
    h.([name{1} '_re']) = reshape(1:36,6,3,2);
    h.([name{1} '_im']) = -h.([name{1} '_re']);
    h.([name{1} '_ma']) = sqrt(2)*h.([name{1} '_re']);
    h.([name{1} '_ph']) = -pi/4*ones(6,3,2);
end
h.ra_t = [0 .02 .04]; h.ra_w = h.w;
h.ex_t = [-.02 0 .02]; h.ex_w = h.w;
h.ra_K = reshape(1:108,6,6,3);
h.ex_K = ones(6,3,3);
re = h.ex_re*h.rho*h.g;
im = -h.ex_im*h.rho*h.g;
end

function verify_failure(testCase,call,identifier,field)
caught = [];
try
    call();
catch failure
    caught = failure;
end
testCase.assertNotEmpty(caught,'Invalid hydrodynamics was accepted');
testCase.verifyEqual(caught.identifier,identifier);
testCase.verifyTrue(contains(caught.message,field),'Failure did not identify the affected field');
end

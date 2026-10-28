function tests = workflowHydrodynamicHandoffTest
% Synthetic coefficient exchange through official BEMIO; no BEM/model/body solve.
tests = functiontests(localfunctions);
end

function setup(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
fork = getenv('TOPWATER_WECSIM_PATH');
if isempty(fork)
    userHome = getenv('USERPROFILE');
    if isempty(userHome)
        userHome = getenv('HOME');
    end
    fork = fullfile(userHome,'.cache','topwater','wec-sim');
end
assert(isfolder(fork),'The maintained WEC-Sim checkout is required');
originalPath = path;
testCase.addTeardown(@() path(originalPath));
addpath(genpath(fullfile(fork,'source')));
addpath(fullfile(root,'src','topwater','data','wecsim'));
testCase.assertEqual(which('readCAPYTAINE'), ...
    fullfile(fork,'source','functions','BEMIO','readCAPYTAINE.m'));
fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
originalDirectory = pwd;
testCase.addTeardown(@() cd(originalDirectory));
cd(fixture.Folder);
[prepared,native] = coefficient_fixture(fixture.Folder);
testCase.TestData.prepared = prepared;
testCase.TestData.native = native;
end

function testGenuineRawSubsetRetainsFullRadiationKernel(testCase)
p = testCase.TestData.prepared;
reference = readCAPYTAINE(struct(),p.netcdf_path);
reference = radiationIRF(reference,10,501,2401,[],[]);
h = prepare_workflow_hydrodynamics(p);
testCase.verifyEqual(h.w,[.1 .2]);
testCase.verifyEqual(h.Nh,25);
testCase.verifyEqual(h.theta(1),0);
testCase.verifyEqual(h.theta(end),360,'AbsTol',1e-12);
testCase.verifyEqual(h.B,zeros(6,6,2));
testCase.verifyGreaterThan(max(abs(reference.ra_K),[],'all'),0);
testCase.verifyEqual(h.ra_K,reference.ra_K);
testCase.verifyEqual(h.ra_t,reference.ra_t);
testCase.verifyEqual(h.ra_w,reference.ra_w);
testCase.verifyEqual(h.ra_w(end),.3,'AbsTol',1e-14);
testCase.verifyEqual(h.Ainf,p.n.Ainf/p.n.rho);
readback = readBEMIOH5([h.file '.h5'],1,0);
testCase.verifyEqual(readback.hydro_coeffs.radiation_damping.impulse_response_fun.K,reference.ra_K);
testCase.verifyEqual(readback.hydro_coeffs.added_mass.inf_freq,h.Ainf);
testCase.verifyEqual(readback.hydro_coeffs.added_mass.all,h.A);
end

function testNativePhaseAndBridgeMatchIndependentKernelFormula(testCase)
p = testCase.TestData.prepared;
h = prepare_workflow_hydrodynamics(p);
testCase.verifyEqual(h.ex_w(1),0);
testCase.verifyEqual(h.ex_w(end),.2,'AbsTol',1e-14);
testCase.verifyEqual(h.ex_t,[-.1:.02:.1],'AbsTol',1e-14);
mode = 4;
direction = 4;
theta = h.theta(direction)*pi/180;
omega = h.ex_w;
transfer = analytic_transfer(max(omega,.1),theta,mode);
low = omega < .1;
limit = p.n.K(mode,3);
transfer(low) = limit+(omega(low)/.1).^2.*(analytic_transfer(.1,theta,mode)-limit);
for index = [1 ceil(numel(h.ex_t)/2) numel(h.ex_t)]
    lag = h.ex_t(index);
    expected = trapz(omega,real(transfer).*cos(omega*lag) ...
        +imag(transfer).*sin(omega*lag))/(pi*p.n.rho*p.n.g);
    testCase.verifyEqual(h.ex_K(mode,direction,index),expected,'AbsTol',1e-13);
end
native = testCase.TestData.native;
expected = permute(native.excitation_native,[3 2 1]);
testCase.verifyEqual(h.ex_re*p.n.rho*p.n.g,real(expected),'AbsTol',1e-12);
testCase.verifyEqual(h.ex_im*p.n.rho*p.n.g,-imag(expected),'AbsTol',1e-12);
end

function testNonfiniteDenseRecordFailsBeforeH5(testCase)
p = testCase.TestData.prepared;
native = testCase.TestData.native;
native.excitation_native(1,2,4) = NaN;
save(p.native_mat_path,'-struct','native','-v7');
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'excitation_native');
testCase.verifyFalse(isfile('canonical.h5'));
end

function testUnrecordedRawFrequencyIsNotFabricated(testCase)
p = testCase.TestData.prepared;
native = testCase.TestData.native;
native.excitation_omega_rad_s(2) = .25;
save(p.native_mat_path,'-struct','native','-v7');
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'genuine recorded radiation frequencies');
testCase.verifyFalse(isfile('canonical.h5'));
end

function testComponentIdentityAndPeriodicSeamAreRequired(testCase)
p = testCase.TestData.prepared;
native = testCase.TestData.native;
native.excitation_native(1,2,4) = native.excitation_native(1,2,4)+1;
save(p.native_mat_path,'-struct','native','-v7');
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'Froude-Krylov plus diffraction');
native = testCase.TestData.native;
native.excitation_diffraction_native(1,end,4) = native.excitation_diffraction_native(1,end,4)+1;
save(p.native_mat_path,'-struct','native','-v7');
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'periodic seam');
end

function testNativePhysicalIdentityMustMatchPreparedCase(testCase)
p = testCase.TestData.prepared;
native = testCase.TestData.native;
native.cg(1) = 1;
save(p.native_mat_path,'-struct','native','-v7');
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'Native/prepared cg');
testCase.verifyFalse(isfile('canonical.h5'));
end

function testStandardSpectralStructPreservesConsumedRadiationAndPhysics(testCase)
p = testCase.TestData.prepared;
legacy = prepare_workflow_hydrodynamics(p);
expected = readBEMIOH5([legacy.file '.h5'],1,0);
delete([legacy.file '.h5']);
p.spectral_hydro = spectral_struct_fixture(p.n);
p = rmfield(p,{'excitation_irf_half_width_s','excitation_irf_dt_s', ...
    'excitation_irf_frequency_step_rad_s'});
[h,nativeHydro] = prepare_workflow_hydrodynamics(p);
testCase.verifyFalse(isfile('canonical.h5'));
testCase.verifyFalse(isfield(h,'ex_K'));
testCase.verifyFalse(isfield(nativeHydro.hydro_coeffs.excitation,'impulse_response_fun'));
testCase.verifyFalse(isfield(nativeHydro.hydro_coeffs.added_mass,'all'));
testCase.verifyFalse(isfield(nativeHydro.hydro_coeffs.radiation_damping,'all'));
testCase.verifyFalse(isfield(nativeHydro.simulation_parameters,'T'));
testCase.verifyEqual(nativeHydro.simulation_parameters.w,p.spectral_hydro.simulation_parameters.w(:)');
testCase.verifyEqual(nativeHydro.simulation_parameters.direction,p.spectral_hydro.simulation_parameters.direction(:)');
testCase.verifyEqual(nativeHydro.hydro_coeffs.excitation,p.spectral_hydro.hydro_coeffs.excitation);
testCase.verifyEqual(nativeHydro.hydro_coeffs.mean_drift,p.spectral_hydro.hydro_coeffs.mean_drift);
testCase.verifyEqual(nativeHydro.properties.centerGravity,expected.properties.centerGravity);
testCase.verifyEqual(nativeHydro.properties.centerBuoyancy,expected.properties.centerBuoyancy);
testCase.verifyEqual(nativeHydro.properties.volume,expected.properties.volume);
testCase.verifyEqual(nativeHydro.hydro_coeffs.linear_restoring_stiffness,expected.hydro_coeffs.linear_restoring_stiffness);
testCase.verifyEqual(nativeHydro.hydro_coeffs.added_mass.inf_freq,expected.hydro_coeffs.added_mass.inf_freq);
testCase.verifyEqual(nativeHydro.hydro_coeffs.radiation_damping.impulse_response_fun, ...
    expected.hydro_coeffs.radiation_damping.impulse_response_fun);
testCase.verifyEqual(h.ra_K,legacy.ra_K);
testCase.verifyEqual(h.ra_t,legacy.ra_t);
testCase.verifyEqual(h.ra_w,legacy.ra_w);
testCase.verifyEqual(h.Ainf,legacy.Ainf);
testCase.verifyEqual(h.A,legacy.A);
testCase.verifyEqual(h.B,legacy.B);
% Exercise this exact zero-boundary transport through ordinary native preparation.
frequencies = [.05;.1;.15]/(2*pi);
directions = [0;45];
phase = [.2 -.3;.7 .1;-.4 .8];
spectrum = ones(3,1);
spread = ones(3,2);
file = fullfile(pwd,'standard-spectrum.mat');
save(file,'frequencies','directions','phase','spectrum','spread');
simu = simulationClass;
simu.dt = .02; simu.endTime = .04; simu.rampTime = 0;
simu.rho = p.n.rho; simu.gravity = p.n.g; simu.cicEndTime = 10;
simu.outputDir = fullfile(pwd,'unused-simulation-output');
simu.simMechanicsFile = fullfile(pwd,'not-loaded.slx');
simu.setup();
body = bodyClass(nativeHydro);
body.mass = p.n.M(1,1); body.inertia = diag(p.n.M(4:6,4:6))';
body.setDOF(1,0);
body.loadHydroData(body.hydroStruct,1);
waves = waveClass('spectrumImportFullDir');
waves.spectrumFile = file;
waves.setup(nativeHydro.simulation_parameters.w,'infinite',0,simu.dt, ...
    simu.maxIt,simu.time,simu.gravity,simu.rho);
body.checkInputs('off',simu.stateSpace,0,waves.typeNum);
body.hydroForcePre(waves,simu,1);
testCase.verifyEqual(waves.omega,frequencies.*2.*pi);
testCase.verifyEqual(waves.phase,phase);
testCase.verifyEqual(body.hydroForce.hf1.fAddedMass,legacy.Ainf.*p.n.rho);
expectedCIC = zeros(numel(simu.cicTime),6,6);
dimensionalKernel = legacy.ra_K.*p.n.rho;
for row = 1:6
    for column = 1:6
        expectedCIC(:,row,column) = interp1(legacy.ra_t, ...
            squeeze(dimensionalKernel(row,column,:)),simu.cicTime,'spline');
    end
end
testCase.verifyEqual(body.hydroForce.hf1.irkb,expectedCIC);
end

function testInvalidStandardSpectralInputFailsBeforeExport(testCase)
p = testCase.TestData.prepared;
standard = spectral_struct_fixture(p.n);
p.spectral_hydro = standard;
p.spectral_hydro.hydro_coeffs.excitation.re(1,2,2) = NaN;
testCase.verifyError(@() prepare_workflow_hydrodynamics(p),'Topwater:NonfiniteHydrodynamics');
testCase.verifyFalse(isfile('canonical.h5'));
p.spectral_hydro = standard;
p.spectral_hydro.hydro_coeffs.mean_drift(1,2,2) = .01;
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'mean drift');
p.spectral_hydro = standard;
p.spectral_hydro.hydro_coeffs.excitation.re(1,2,1) = ...
    p.spectral_hydro.hydro_coeffs.excitation.re(1,2,1)+.01;
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'zero limit');
testCase.verifyFalse(isfile('canonical.h5'));
end

function h = spectral_struct_fixture(n)
omega = [0 .05 .1 .15 .2];
headings = [-180 0 45 180];
h.simulation_parameters = struct('w',omega,'direction',headings);
h.hydro_coeffs.excitation = struct('re',zeros(6,numel(headings),numel(omega)), ...
    'im',zeros(6,numel(headings),numel(omega)));
h.hydro_coeffs.mean_drift = zeros(6,numel(headings),numel(omega));
for mode = 1:6
    for direction = 1:numel(headings)
        theta = mod(headings(direction)*pi/180,2*pi);
        for iw = 1:numel(omega)
            value = analytic_transfer(max(omega(iw),.1),theta,mode);
            if omega(iw) < .1
                value = n.K(mode,3)+(omega(iw)/.1)^2*(value-n.K(mode,3));
            end
            h.hydro_coeffs.excitation.re(mode,direction,iw) = real(value)/(n.rho*n.g);
            h.hydro_coeffs.excitation.im(mode,direction,iw) = -imag(value)/(n.rho*n.g);
        end
    end
end
end

function [p,native] = coefficient_fixture(folder)
radOmega = [.1 .2 .3];
importHeadings = [0 pi/4 pi/2];
omega = [.1 .2];
headings = linspace(0,2*pi,25);
n = struct('M',blkdiag(2*eye(3),diag([3 4 5])), ...
    'K',diag([0 0 100 20 30 0]),'Ainf',reshape(1:36,6,6)/100, ...
    'cg',[0 0 0],'cb',[0 0 -.1],'volume',.002,'rho',1000,'g',10);
n.K(3,4) = .5; n.K(4,3) = .5;
[fk,sc] = dense_table(omega,headings);
native = n;
native.excitation_omega_rad_s = omega;
native.excitation_headings_rad = headings;
native.excitation_froude_krylov_native = fk;
native.excitation_diffraction_native = sc;
native.excitation_native = fk+sc;
native.excitation_zero_limit_native = n.K(:,3);
p = struct('n',n,'native_mat_path',fullfile(folder,'native.mat'), ...
    'netcdf_path',fullfile(folder,'canonical.nc'),'excitation_irf_half_width_s',.1, ...
    'excitation_irf_dt_s',.02,'excitation_irf_frequency_step_rad_s',.005);
save(p.native_mat_path,'-struct','native','-v7');
write_synthetic_netcdf(p.netcdf_path,n,radOmega,importHeadings);
end

function [fk,sc] = dense_table(omega,headings)
fk = complex(zeros(numel(omega),numel(headings),6));
sc = fk;
for iw = 1:numel(omega)
    for direction = 1:numel(headings)
        for mode = 1:6
            transfer = analytic_transfer(omega(iw),headings(direction),mode);
            fk(iw,direction,mode) = .4*transfer;
            sc(iw,direction,mode) = .6*transfer;
        end
    end
end
if abs(headings(end)-2*pi) < 1e-12
    fk(:,end,:) = fk(:,1,:);
    sc(:,end,:) = sc(:,1,:);
end
end

function value = analytic_transfer(omega,theta,mode)
value = mode*(2+omega)*(1+.25*cos(theta)) ...
    +1i*mode*(.5+omega)*(1+.2*sin(theta));
end

function write_synthetic_netcdf(file,n,omega,headings)
nccreate(file,'rho','Datatype','double','Format','classic');
ncwrite(file,'rho',n.rho);
for name = {'g','water_depth','volume'}
    nccreate(file,name{1},'Datatype','double');
end
ncwrite(file,'g',n.g); ncwrite(file,'water_depth',Inf); ncwrite(file,'volume',n.volume);
write_array(file,'omega',omega(:),{'omega',numel(omega)});
write_array(file,'wave_direction',headings(:),{'wave_direction',numel(headings)});
write_array(file,'center_of_mass',n.cg(:),{'xyz',3});
write_array(file,'center_of_buoyancy',n.cb(:),{'xyz',3});
write_array(file,'hydrostatic_stiffness',n.K,{'influenced_dof',6,'radiating_dof',6});
write_labels(file,'body_name',{'boat'},'body_strlen','body_count');
write_labels(file,'influenced_dof',{'Surge','Sway','Heave','Roll','Pitch','Yaw'},'mode_strlen','influenced_dof');
write_labels(file,'radiating_dof',{'Surge','Sway','Heave','Roll','Pitch','Yaw'},'mode_strlen','radiating_dof');
write_labels(file,'complex',{'re','im'},'complex_strlen','complex');
A = zeros(6,6,numel(omega)); B = A;
for iw = 1:numel(omega)
    A(:,:,iw) = n.rho*(reshape(1:36,6,6)+iw);
    if iw == numel(omega)
        B(:,:,iw) = n.rho*omega(iw)*reshape(1:36,6,6);
    end
end
dimensions = {'influenced_dof',6,'radiating_dof',6,'omega',numel(omega)};
write_array(file,'added_mass',A,dimensions);
write_array(file,'radiation_damping',B,dimensions);
[fk,sc] = dense_table(omega,headings);
dimensions = {'influenced_dof',6,'wave_direction',numel(headings),'omega',numel(omega),'complex',2};
fk = permute(fk,[3 2 1]); sc = permute(sc,[3 2 1]);
write_array(file,'Froude_Krylov_force',cat(4,real(fk),imag(fk)),dimensions);
write_array(file,'diffraction_force',cat(4,real(sc),imag(sc)),dimensions);
write_sidecars(fileparts(file),n);
end

function write_sidecars(folder,n)
fid = fopen(fullfile(folder,'Hydrostatics.dat'),'w');
assert(fid ~= -1,'Cannot write owned hydrostatics fixture');
cleanup = onCleanup(@() fclose(fid));
for j = 1:3
    fprintf(fid,'XF = %.17g - XG = %.17g\n',n.cb(j),n.cg(j));
end
fprintf(fid,'Displacement = %.17g\n',n.volume);
clear cleanup
fid = fopen(fullfile(folder,'KH.dat'),'w');
assert(fid ~= -1,'Cannot write owned stiffness fixture');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%.17g %.17g %.17g %.17g %.17g %.17g\n',n.K');
end

function write_array(file,name,value,dimensions)
nccreate(file,name,'Dimensions',dimensions,'Datatype','double');
ncwrite(file,name,value);
end

function write_labels(file,name,labels,lengthName,countName)
width = max(cellfun(@numel,labels));
value = char(zeros(width,numel(labels)));
for j = 1:numel(labels)
    value(1:numel(labels{j}),j) = labels{j}';
end
nccreate(file,name,'Dimensions',{lengthName,width,countName,numel(labels)},'Datatype','char');
ncwrite(file,name,value);
end

function verify_rejection(testCase,operation,message)
failure = [];
try
    operation();
catch caught
    failure = caught;
end
testCase.assertNotEmpty(failure,'Invalid coefficient handoff was accepted');
testCase.verifyTrue(contains(failure.message,message),failure.message);
end

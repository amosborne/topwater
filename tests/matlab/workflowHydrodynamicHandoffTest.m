function tests = workflowHydrodynamicHandoffTest
% Synthetic coefficient exchange through official BEMIO; no BEM/model/body solve.
tests = functiontests(localfunctions);
end

function setup(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
fork = getenv('TOPWATER_WECSIM_PATH');
if isempty(fork)
    fork = fullfile(root,'artifacts','dependencies','topwater-wec-sim');
end
assert(isfolder(fork),'A coordinator-frozen WEC-Sim checkout is required');
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

function testSpectralViewOmitsUnusedIRFAndPreservesEveryRawRadiationField(testCase)
p = testCase.TestData.prepared;
legacy = prepare_workflow_hydrodynamics(p);
expected = readBEMIOH5([legacy.file '.h5'],1,0);
expected.hydro_coeffs.excitation = rmfield(expected.hydro_coeffs.excitation,'impulse_response_fun');
p.component_view = spectral_view_fixture(p.n);
p = rmfield(p,{'excitation_irf_half_width_s','excitation_irf_dt_s', ...
    'excitation_irf_frequency_step_rad_s'});
[h,nativeHydro] = prepare_workflow_hydrodynamics(p);
testCase.verifyFalse(isfield(h,'ex_K'));
testCase.verifyFalse(isfield(h,'ex_t'));
testCase.verifyFalse(isfield(h,'ex_w'));
testCase.verifyFalse(isfield(nativeHydro.hydro_coeffs.excitation,'impulse_response_fun'));
testCase.verifyEqual(nativeHydro.hydro_coeffs.excitation.component_view,p.component_view);
actual = nativeHydro;
actual.hydro_coeffs.excitation = rmfield(actual.hydro_coeffs.excitation,'component_view');
testCase.verifyEqual(actual,expected);
testCase.verifyEqual(h.ra_K,legacy.ra_K);
testCase.verifyEqual(h.ra_t,legacy.ra_t);
testCase.verifyEqual(h.ra_w,legacy.ra_w);
testCase.verifyEqual(h.Ainf,legacy.Ainf);
info = h5info([h.file '.h5'],'/body1/hydro_coeffs/excitation');
testCase.verifyFalse(any(endsWith({info.Groups.Name},'/impulse_response_fun')));
end

function testPartialWriterIRFIsFatalBeforeExistingOutputChanges(testCase)
h = prepare_workflow_hydrodynamics(testCase.TestData.prepared);
before = read_file_bytes([h.file '.h5']);
partial = rmfield(h,'ex_t');
testCase.verifyError(@() writeBEMIOH5(partial),'WECsim:IncompleteExcitationIRF');
testCase.verifyEqual(read_file_bytes([h.file '.h5']),before);
end

function testInvalidSpectralViewFailsBeforeH5(testCase)
p = testCase.TestData.prepared;
view = spectral_view_fixture(p.n);
p.component_view = view;
p.component_view.re(1,2,2) = NaN;
testCase.verifyError(@() prepare_workflow_hydrodynamics(p),'Topwater:NonfiniteHydrodynamics');
testCase.verifyFalse(isfile('canonical.h5'));
p.component_view = view;
p.component_view.md(1,2,2) = .01;
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'mean drift');
testCase.verifyFalse(isfile('canonical.h5'));
p.component_view = view;
p.component_view.re(1,2,1) = p.component_view.re(1,2,1)+.01;
verify_rejection(testCase,@() prepare_workflow_hydrodynamics(p),'zero limit');
testCase.verifyFalse(isfile('canonical.h5'));
end

function view = spectral_view_fixture(n)
omega = [0 .05 .1 .15 .2];
headings = [-180 0 45 180];
view = struct('omega_rad_s',omega(:),'headings_deg',headings(:), ...
    're',zeros(6,numel(headings),numel(omega)), ...
    'im',zeros(6,numel(headings),numel(omega)), ...
    'md',zeros(6,numel(headings),numel(omega)));
for mode = 1:6
    for direction = 1:numel(headings)
        theta = mod(headings(direction)*pi/180,2*pi);
        for iw = 1:numel(omega)
            value = analytic_transfer(max(omega(iw),.1),theta,mode);
            if omega(iw) < .1
                value = n.K(mode,3)+(omega(iw)/.1)^2*(value-n.K(mode,3));
            end
            view.re(mode,direction,iw) = real(value)/(n.rho*n.g);
            view.im(mode,direction,iw) = -imag(value)/(n.rho*n.g);
        end
    end
end
end

function value = read_file_bytes(file)
fid = fopen(file,'rb');
assert(fid ~= -1,'Cannot read owned test output');
cleanup = onCleanup(@() fclose(fid));
value = fread(fid,Inf,'*uint8');
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

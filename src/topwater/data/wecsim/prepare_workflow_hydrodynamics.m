function [h,nativeHydro] = prepare_workflow_hydrodynamics(prepared)
% Build native WEC excitation inputs without computing an applied force history.
n = prepared.n;
native = load(prepared.native_mat_path);
verify_native_identity(native,n);
if isfield(prepared,'spectral_hydro')
    transport = prepared.spectral_hydro;
    validate_workflow_hydrodynamics(transport,'spectral_input');
    zero = repmat(n.K(:,3)/(n.rho*n.g),1,numel(transport.simulation_parameters.direction));
    assert(isequal(transport.hydro_coeffs.excitation.re(:,:,1),zero) ...
        && all(transport.hydro_coeffs.excitation.im(:,:,1) == 0,'all'), ...
        'Standard coefficient input zero limit differs from physical heave stiffness');
end
h = readCAPYTAINE(struct(),prepared.netcdf_path);
assert(h.Nb == 1 && all(h.dof == 6),'The workflow requires one six-mode body');
assert(h.rho == n.rho && h.g == n.g && isinf(h.h),'BEMIO physical parameters differ');
[rawReal,rawImaginary] = raw_excitation(prepared.netcdf_path);
validate_workflow_hydrodynamics(h,'imported',rawReal,rawImaginary);
require_close(h.ex_re*n.rho*n.g,rawReal,'Imported excitation real convention');
require_close(h.ex_im*n.rho*n.g,-rawImaginary,'Imported excitation imaginary convention');
h.cg = n.cg(:); h.cb = n.cb(:); h.Vo = n.volume;
h.Khs = n.K/(n.rho*n.g);
% This unchanged full-band transform precedes any common-grid raw-table subset.
fullRadiation = radiationIRF(h,10,501,2401,[],[]);
fullRadiation.Ainf = n.Ainf/n.rho;
[omega,headings,fk,sc,total,zeroLimit] = dense_excitation(native,n);
[found,frequencyIndex] = ismember(omega(:),fullRadiation.w(:));
assert(all(found),'Excitation raw frequencies must be genuine recorded radiation frequencies');
originalHeadings = fullRadiation.theta(:);
headingIndex = zeros(size(originalHeadings));
for j = 1:numel(originalHeadings)
    matches = find(abs(headings(:)*180/pi-originalHeadings(j)) <= 1e-10);
    assert(numel(matches) == 1,'Original import heading is missing or ambiguous in dense transfer data');
    headingIndex(j) = matches;
end
% Shared points bind producer axes/units/signs to the independently read netCDF.
require_close(real(total(:,headingIndex,:)), ...
    permute(h.ex_re(:,:,frequencyIndex)*n.rho*n.g,[3 2 1]),'Dense/native shared real transfer');
require_close(imag(total(:,headingIndex,:)), ...
    -permute(h.ex_im(:,:,frequencyIndex)*n.rho*n.g,[3 2 1]),'Dense/native shared imaginary transfer');
require_close(real(fk(:,headingIndex,:)), ...
    permute(h.fk_re(:,:,frequencyIndex)*n.rho*n.g,[3 2 1]),'Dense/native shared Froude-Krylov real transfer');
require_close(imag(fk(:,headingIndex,:)), ...
    -permute(h.fk_im(:,:,frequencyIndex)*n.rho*n.g,[3 2 1]),'Dense/native shared Froude-Krylov imaginary transfer');
require_close(real(sc(:,headingIndex,:)), ...
    permute(h.sc_re(:,:,frequencyIndex)*n.rho*n.g,[3 2 1]),'Dense/native shared diffraction real transfer');
require_close(imag(sc(:,headingIndex,:)), ...
    -permute(h.sc_im(:,:,frequencyIndex)*n.rho*n.g,[3 2 1]),'Dense/native shared diffraction imaginary transfer');
h = fullRadiation;
h.w = omega(:)'; h.T = 2*pi./h.w; h.Nf = numel(h.w);
h.theta = headings(:)'*180/pi; h.Nh = numel(h.theta);
h.A = fullRadiation.A(:,:,frequencyIndex);
h.B = fullRadiation.B(:,:,frequencyIndex);
h = set_excitation_coefficients(h,'fk',fk);
h = set_excitation_coefficients(h,'sc',sc);
h = set_excitation_coefficients(h,'ex',total);
validate_workflow_hydrodynamics(h,'imported', ...
    permute(real(total),[3 2 1]),permute(imag(total),[3 2 1]));
if isfield(prepared,'spectral_hydro')
    validate_workflow_hydrodynamics(h,'spectral');
    % Official BEMIO conversion, followed by ordinary public input preparation.
    % CIC consumes Ainf and its own kernel/time, not raw frequency A/B tables.
    nativeHydro = rebuildHydroStruct(h,1,0);
    nativeHydro.simulation_parameters = rmfield(nativeHydro.simulation_parameters,'T');
    nativeHydro.simulation_parameters.w = transport.simulation_parameters.w(:)';
    nativeHydro.simulation_parameters.direction = transport.simulation_parameters.direction(:)';
    nativeHydro.hydro_coeffs.excitation = transport.hydro_coeffs.excitation;
    nativeHydro.hydro_coeffs.mean_drift = transport.hydro_coeffs.mean_drift;
    nativeHydro.hydro_coeffs.added_mass = rmfield(nativeHydro.hydro_coeffs.added_mass,'all');
    nativeHydro.hydro_coeffs.radiation_damping = rmfield(nativeHydro.hydro_coeffs.radiation_damping,'all');
    assert(isequal(nativeHydro.hydro_coeffs.added_mass.inf_freq,fullRadiation.Ainf) ...
        && isequal(nativeHydro.hydro_coeffs.radiation_damping.impulse_response_fun.K,fullRadiation.ra_K) ...
        && isequal(nativeHydro.hydro_coeffs.radiation_damping.impulse_response_fun.t(:),fullRadiation.ra_t(:)), ...
        'Standard struct preparation changed the consumed radiation operator');
    fprintf('WEC-Sim standard struct handoff: %i derived frequency nodes, %i headings; CIC/Ainf retained.\n', ...
        numel(nativeHydro.simulation_parameters.w),numel(nativeHydro.simulation_parameters.direction));
    return
end
halfWidth = positive_scalar(prepared.excitation_irf_half_width_s,'excitation IRF half width');
dt = positive_scalar(prepared.excitation_irf_dt_s,'excitation IRF time step');
dw = positive_scalar(prepared.excitation_irf_frequency_step_rad_s,'excitation IRF frequency step');
halfSteps = round(halfWidth/dt);
assert(halfSteps >= 1 && abs(halfSteps*dt-halfWidth) <= 1e-10, ...
    'Excitation IRF half width must be an integer multiple of its time step');
nFrequency = max(2,ceil(omega(end)/dw)+1);
kernelOmega = linspace(0,omega(end),nFrequency);
kernelInput = struct('dof',6,'Nh',h.Nh,'w',kernelOmega, ...
    'ex_re',zeros(6,h.Nh,nFrequency),'ex_im',zeros(6,h.Nh,nFrequency));
low = kernelOmega < omega(1);
query = max(kernelOmega,omega(1));
for mode = 1:6
    for direction = 1:h.Nh
        realPart = interp1(h.w,squeeze(h.ex_re(mode,direction,:)),query,'linear');
        imaginaryPart = interp1(h.w,squeeze(h.ex_im(mode,direction,:)),query,'linear');
        limit = zeroLimit(mode)/(n.rho*n.g);
        realPart(low) = limit+(kernelOmega(low)/omega(1)).^2.*(realPart(low)-limit);
        imaginaryPart(low) = (kernelOmega(low)/omega(1)).^2.*imaginaryPart(low);
        kernelInput.ex_re(mode,direction,:) = realPart;
        kernelInput.ex_im(mode,direction,:) = imaginaryPart;
    end
end
% BEMIO receives only the independent excitation quadrature grid, including the
% declared zero-frequency limit. No radiation coefficient is fabricated at zero.
kernel = excitationIRF(kernelInput,halfWidth,2*halfSteps+1,nFrequency,0,omega(end));
h.ex_K = kernel.ex_K; h.ex_t = kernel.ex_t; h.ex_w = kernel.ex_w;
assert(isequal(h.ra_K,fullRadiation.ra_K) && isequal(h.ra_t,fullRadiation.ra_t) ...
    && isequal(h.ra_w,fullRadiation.ra_w) && isequal(h.Ainf,fullRadiation.Ainf), ...
    'Dense excitation preparation changed the full-band radiation operator');
validate_workflow_hydrodynamics(h,'prepared');
h.file = fullfile(pwd,'canonical');
writeBEMIOH5(h);
nativeHydro = verify_h5_roundtrip(h,fullRadiation);
fprintf('WEC-Sim excitation IRF: +/-%.8g s, dt %.8g s, omega 0..%.8g rad/s, domega %.8g rad/s; %i headings.\n', ...
    halfWidth,dt,omega(end),kernelOmega(2)-kernelOmega(1),h.Nh);
fprintf('WEC-Sim radiation IRF retained: omega %.8g..%.8g rad/s, time %.8g..%.8g s.\n', ...
    fullRadiation.ra_w(1),fullRadiation.ra_w(end),fullRadiation.ra_t(1),fullRadiation.ra_t(end));
fprintf('WEC-Sim raw H5 coefficient grid: %i frequencies, omega %.8g..%.8g rad/s; no unsolved high-frequency headings.\n', ...
    h.Nf,h.w(1),h.w(end));
end

function [omega,headings,fk,sc,total,zeroLimit] = dense_excitation(native,n)
omega = native.excitation_omega_rad_s(:);
headings = native.excitation_headings_rad(:);
require_array(omega,'native excitation frequency grid',true);
require_array(headings,'native excitation heading grid',true);
assert(numel(omega) >= 2 && all(diff(omega) > 0) && abs(omega(1)-.1) <= 1e-12, ...
    'Dense excitation frequency grid must increase from the declared 0.1 rad/s bridge anchor');
assert(numel(headings) >= 2 && all(diff(headings) > 0) && abs(headings(1)) <= 1e-12 ...
    && abs(headings(end)-2*pi) <= 1e-12,'Dense excitation headings must include the periodic 0/2pi seam');
fk = native.excitation_froude_krylov_native;
sc = native.excitation_diffraction_native;
total = native.excitation_native;
shape = [numel(omega) numel(headings) 6];
for name = {'excitation_froude_krylov_native','excitation_diffraction_native','excitation_native'}
    value = native.(name{1});
    require_array(value,name{1},false);
    assert(isequal(size(value),shape),'Dense excitation arrays must be frequency/heading/mode ordered');
    require_close(value(:,1,:),value(:,end,:),[name{1} ' periodic seam']);
end
require_close(fk+sc,total,'Dense total excitation must equal Froude-Krylov plus diffraction');
zeroLimit = native.excitation_zero_limit_native(:);
require_array(zeroLimit,'native zero-frequency excitation limit',true);
assert(numel(zeroLimit) == 6,'Zero-frequency excitation limit must have six components');
require_close(zeroLimit,n.K(:,3),'Zero-frequency excitation limit differs from physical heave stiffness');
end

function h = set_excitation_coefficients(h,prefix,value)
native = permute(value,[3 2 1]);
% Capytaine uses exp(-i*omega*t). BEMIO performs this convention conversion once.
h.([prefix '_re']) = real(native)/(h.rho*h.g);
h.([prefix '_im']) = -imag(native)/(h.rho*h.g);
h.([prefix '_ma']) = abs(native)/(h.rho*h.g);
h.([prefix '_ph']) = angle(conj(native));
end

function verify_native_identity(native,n)
for name = {'M','K','Ainf','cg','cb','volume','rho','g'}
    field = name{1};
    assert(isfield(native,field) && isfield(n,field),'Native coefficient identity field is missing');
    require_array(native.(field),['native.' field],true);
    require_array(n.(field),['prepared.' field],true);
    actual = native.(field); expected = n.(field);
    require_close(actual(:),expected(:),['Native/prepared ' field]);
end
end

function [re,im] = raw_excitation(file)
names = erase(lower(string(ncread(file,'influenced_dof')')),char(0));
modes = ["surge","sway","heave","roll","pitch","yaw"];
index = zeros(1,6);
for j = 1:6, index(j) = find(names == modes(j)); end
fk = named_array(file,'Froude_Krylov_force',{'influenced_dof','wave_direction','omega','complex'});
sc = named_array(file,'diffraction_force',{'influenced_dof','wave_direction','omega','complex'});
components = erase(string(ncread(file,'complex')'),char(0));
re = fk(index,:,:,components == "re")+sc(index,:,:,components == "re");
im = fk(index,:,:,components == "im")+sc(index,:,:,components == "im");
end

function value = named_array(file,name,order)
info = ncinfo(file,name);
names = {info.Dimensions.Name};
permutation = zeros(1,numel(order));
for j = 1:numel(order), permutation(j) = find(strcmp(names,order{j})); end
permutation = [permutation setdiff(1:numel(names),permutation,'stable')];
value = permute(ncread(file,name),permutation);
shape = size(value);
assert(all(shape(numel(order)+1:end) == 1));
value = reshape(value,shape(1:numel(order)));
end

function readback = verify_h5_roundtrip(h,fullRadiation)
file = [h.file '.h5'];
readback = readBEMIOH5(file,1,0);
assert(isequal(readback.hydro_coeffs.radiation_damping.impulse_response_fun.K,fullRadiation.ra_K) ...
    && isequal(readback.hydro_coeffs.radiation_damping.impulse_response_fun.t(:),fullRadiation.ra_t(:)) ...
    && isequal(readback.hydro_coeffs.added_mass.inf_freq,fullRadiation.Ainf), ...
    'H5 roundtrip changed the full-band radiation kernel or direct infinity added mass');
radiationOmega = h5read(file,'/body1/hydro_coeffs/radiation_damping/impulse_response_fun/w');
assert(isequal(radiationOmega(:),fullRadiation.ra_w(:)), ...
    'H5 roundtrip changed radiation quadrature frequency metadata');
excitationOmega = h5read(file,'/body1/hydro_coeffs/excitation/impulse_response_fun/w');
assert(isequal(readback.hydro_coeffs.excitation.impulse_response_fun.f,h.ex_K) ...
    && isequal(readback.hydro_coeffs.excitation.impulse_response_fun.t(:),h.ex_t(:)), ...
    'H5 roundtrip changed the native excitation kernel');
assert(isequal(excitationOmega(:),h.ex_w(:)), ...
    'H5 roundtrip changed excitation quadrature frequency metadata');
assert(isequal(readback.simulation_parameters.w(:),h.w(:)) ...
    && isequal(readback.simulation_parameters.direction(:),h.theta(:)), ...
    'H5 roundtrip changed raw coefficient reference grids');
assert(isequal(readback.hydro_coeffs.excitation.re,h.ex_re) ...
    && isequal(readback.hydro_coeffs.excitation.im,h.ex_im) ...
    && isequal(readback.hydro_coeffs.added_mass.all,h.A) ...
    && isequal(readback.hydro_coeffs.radiation_damping.all,h.B), ...
    'H5 roundtrip changed raw coefficient values or axes');
for prefix = {'froude-krylov','scattering'}
    if strcmp(prefix{1},'froude-krylov'), field = 'fk'; else, field = 'sc'; end
    for part = {'re','im'}
        value = reverseDimensionOrder(h5read(file, ...
            ['/body1/hydro_coeffs/excitation/' prefix{1} '/' part{1}]));
        assert(isequal(value,h.([field '_' part{1}])),'H5 roundtrip changed raw FK/scattering coefficients');
    end
end
end

function require_array(value,name,realOnly)
if ~isnumeric(value) || isempty(value) || ~all(isfinite(value),'all') ...
        || (realOnly && ~isreal(value))
    error('Topwater:InvalidHydrodynamics','%s must contain finite compatible numeric values',name);
end
end

function require_close(actual,expected,name)
require_array(actual,name,false);
require_array(expected,name,false);
assert(isequal(size(actual),size(expected)),[name ': incompatible array shape']);
assert(max(abs(actual-expected),[],'all') <= 1e-10*max(1,max(abs(expected),[],'all')),name);
end

function value = positive_scalar(value,name)
assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value) && value > 0, ...
    [name ' must be positive and finite']);
end

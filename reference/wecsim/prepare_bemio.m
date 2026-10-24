function prepare_bemio()
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
native = fullfile(root,'artifacts/reference/minimum-20261024/native-gl2');
out = fullfile(root,'artifacts/reference/minimum-20261024/consumer');
assert(~isfolder(out),'Use a fresh output directory');
mkdir(out);
n = load(fullfile(native,'native.mat'));
nc = fullfile(native,'canonical.nc');
h = readCAPYTAINE(struct(),nc);
assert(h.Nb == 1 && all(h.dof == 6));
assert(h.rho == n.rho && h.g == n.g && isinf(h.h));
% Native legacy hydrostatics are rounded text; retain the exact producer arrays.
rounded = struct('cg',h.cg,'cb',h.cb,'volume',h.Vo,'K',h.Khs*n.rho*n.g);
h.cg = n.cg(:);
h.cb = n.cb(:);
h.Vo = n.volume;
h.Khs = n.K/(n.rho*n.g);
% Independently check native variable names, DOFs and the single phase conversion.
modes = ["surge","sway","heave","roll","pitch","yaw"];
names = erase(lower(string(ncread(nc,'influenced_dof')')),char(0));
idx = zeros(1,6);
for i = 1:6
    idx(i) = find(names == modes(i));
end
fk = named_array(nc,'Froude_Krylov_force',{'influenced_dof','wave_direction','omega','complex'});
sc = named_array(nc,'diffraction_force',{'influenced_dof','wave_direction','omega','complex'});
cx = erase(string(ncread(nc,'complex')'),char(0));
ir = find(cx == "re");
ii = find(cx == "im");
rawRe = fk(idx,:,:,ir)+sc(idx,:,:,ir);
rawIm = fk(idx,:,:,ii)+sc(idx,:,:,ii);
er = max(abs(h.ex_re*n.rho*n.g-rawRe),[],'all');
ei = max(abs(h.ex_im*n.rho*n.g+rawIm),[],'all');
assert(er <= 1e-10*max(1,max(abs(rawRe),[],'all')));
assert(ei <= 1e-10*max(1,max(abs(rawIm),[],'all')));
h = radiationIRF(h,10,501,2401,[],[]);
h = excitationIRF(h,20,2001,2401,[],[]);
% Direct native omega=infinity solve takes precedence over finite-band Ainf fit.
h.Ainf = n.Ainf/n.rho;
h.file = fullfile(out,'canonical');
writeBEMIOH5(h);
save(fullfile(out,'bemio.mat'),'h','n','rounded','-v7');
r = struct('native_directory',native,'official_importer',which('readCAPYTAINE'), ...
    'official_writer',which('writeBEMIOH5'),'phase_real_max_error',er, ...
    'phase_imaginary_sign_max_error',ei,'exact_metadata','native.mat', ...
    'infinite_added_mass','native direct omega=infinity solve', ...
    'radiation_seconds',10,'excitation_lag_seconds',[-20 20],'dt_seconds',0.02);
fid = fopen(fullfile(out,'interface-check.json'),'w');
fprintf(fid,'%s\n',jsonencode(r,PrettyPrint=true));
fclose(fid);
fprintf('Official BEMIO import/write complete; native phase and exact metadata checked.\n');
end

function a = named_array(file,var,order)
info = ncinfo(file,var);
names = {info.Dimensions.Name};
perm = zeros(1,numel(order));
for i = 1:numel(order)
    perm(i) = find(strcmp(names,order{i}));
end
perm = [perm setdiff(1:numel(names),perm,'stable')];
a = permute(ncread(file,var),perm);
sz = size(a);
assert(all(sz(numel(order)+1:end) == 1));
a = reshape(a,sz(1:numel(order)));
end

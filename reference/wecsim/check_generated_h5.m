function check_generated_h5()
root = 'C:/Users/alex/topwater';
base = fullfile(root,'artifacts/reference/generated-bank-consumer-20261024/bemio-attempt1');
s = load(fullfile(base,'bemio.mat'));
file = fullfile(base,'canonical.h5');
assert(~isfile(fullfile(base,'h5-readback.json')));
for name = {'cg','cb'}
    actual = h5read(file,['/body1/properties/' name{1}]);
    expected = s.n.(name{1});
    assert(isequal(actual(:),expected(:)));
end
assert(isequal(h5read(file,'/body1/properties/disp_vol'),s.n.volume));
k = h5read(file,'/body1/hydro_coeffs/linear_restoring_stiffness')';
a = h5read(file,'/body1/hydro_coeffs/added_mass/inf_freq')';
kernel = permute(h5read(file,'/body1/hydro_coeffs/radiation_damping/impulse_response_fun/K'),[3 2 1]);
assert(isequal(k,s.n.K/(s.n.rho*s.n.g)));
assert(isequal(a,s.n.Ainf/s.n.rho));
assert(isequal(kernel,s.h.ra_K));
r = struct('complete',true,'cg_cb_volume_exact',true,'full_K_Ainf_exact',true, ...
    'all36_radiation_kernel_exact',true,'radiation_seconds',40);
fid = fopen(fullfile(base,'h5-readback.json'),'w');
fprintf(fid,'%s\n',jsonencode(r,PrettyPrint=true));
fclose(fid);
fprintf('H5 exact metadata, full K/Ainf and all36 radiation kernel readback passed.\n');
end

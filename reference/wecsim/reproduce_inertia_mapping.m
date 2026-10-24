function reproduce_inertia_mapping()
clear classes
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
out = fullfile(root,'artifacts/reference-inertia-stock-attempt2');
assert(~isfolder(out));
mkdir(out);
native = load(fullfile(root,'artifacts/reference/minimum-20261024/native-gl2/native.mat'));
h0 = readBEMIOH5(fullfile(root,'artifacts/reference/minimum-20261024/consumer/canonical.h5'),1,0);
load_system('WECSim_Lib_Body_Elements');
block = 'WECSim_Lib_Body_Elements/Rigid Body/Hydrodynamic Body/Structure/Body Properties';
expr = get_param(block,'ProductsOfInertia');
assert(strcmp(get_param(block,'MaskType'),'File Solid'));
s = struct('cicTime',0:0.02:10,'dt',0.02,'rho',native.rho,'gravity',native.g, ...
    'stateSpace',0,'b2b',0,'rampTime',0);
Mtest = native.M;
Mtest(4,5) = 11; Mtest(5,4) = 11;
Mtest(4,6) = 13; Mtest(6,4) = 13;
Mtest(5,6) = 17; Mtest(6,5) = 17;
for j = 1:3
    h = h0;
    if j == 2, M = Mtest; else, M = native.M; end
    if j == 3, A = native.Ainf; else, A = zeros(6); end
    h.hydro_coeffs.added_mass.inf_freq = A/native.rho;
    body = bodyClass(h);
    body.setNumber(1);
    body.loadHydroData(h,1);
    body.setDOF(1,0);
    body.mass = M(1,1);
    body.inertia = diag(M(4:6,4:6))';
    body.inertiaProducts = [M(4,5) M(4,6) M(5,6)];
    w = waveClass('noWaveCIC');
    w.setup(h.simulation_parameters.w,inf,0,0.02,500,0:0.02:10,native.g,native.rho);
    body.hydroForcePre(w,s,1);
    body.adjustMassMatrix(0,native.rho);
    hf = body.hydroForce.hf1;
    % Evaluate the actual selected library expression, then use File Solid's
    % documented [Iyz Izx Ixy] convention to reconstruct the physical tensor.
    product = eval(expr);
    solid = zeros(6);
    solid(1:3,1:3) = hf.adjustedMass*eye(3);
    solid(4:6,4:6) = diag(hf.adjustedInertia);
    solid(5,6) = product(1); solid(6,5) = product(1);
    solid(4,6) = product(2); solid(6,4) = product(2);
    solid(4,5) = product(3); solid(5,4) = product(3);
    reconstructed = solid+hf.fAddedMass;
    expected = M+A;
    res(j).input_mass = M;
    res(j).input_added_mass = A;
    res(j).solid_products_yz_zx_xy = product;
    res(j).solid_mass = solid;
    res(j).residual_added_mass = hf.fAddedMass;
    res(j).actual_total = reconstructed;
    res(j).expected_total = expected;
    res(j).max_error = max(abs(reconstructed-expected),[],'all');
end
assert(res(1).max_error < 1e-9,'Canonical dry middle product should be immune');
assert(abs(res(2).max_error-6) < 1e-9,'Distinct dry product permutation expectation differs');
assert(res(3).max_error > 100,'Canonical added products should exercise the order defect');
fid = fopen(fullfile(out,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('library_block',block,'actual_expression',expr, ...
    'block_reference','https://www.mathworks.com/help/releases/R2025b/sm/ref/filesolid.html', ...
    'cases',res,'result','confirmed xy/yz permutation at File Solid boundary'),PrettyPrint=true));
fclose(fid);
close_system('WECSim_Lib_Body_Elements',0);
fprintf('Confirmed File Solid xy/yz permutation in native added-mass path.\n');
end

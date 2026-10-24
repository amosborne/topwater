function reproduce_multiple_imports()
clear classes
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
source = fullfile(root,'artifacts/reference-multiple-imports-stock-source');
addpath(source,'-begin');
assert(contains(which('bodyClass'),'reference-multiple-imports-stock-source'));
out = fullfile(root,'artifacts/reference-multiple-imports-stock-attempt5');
assert(~isfolder(out));
mkdir(out);
t = (-2:12)';
for j = 1:2
    file = fullfile(out,sprintf('wave_%d.txt',j));
    writematrix([t j*double(t == 5)],file,delimiter='tab');
    w(j) = waveClass('elevationImport');
    w(j).elevationFile = file;
    w(j).direction = 0;
    w(j).setup([0.1 12],inf,0,1,10,0:10,1,1);
end
b = paddingProbe();
s = simulationClass();
s.dt = 1;
s.rho = 1;
s.gravity = 1;
s.rampTime = 0;
actual = b.probeMultiple(w,s);
expected = zeros(11,1);
expected(4:8) = 3*[1 2 4 8 16];
assert(max(abs(actual-expected)) == 32,'Stock multidirection expectation differs');
r = struct('body_source',which('bodyClass'),'actual',actual','expected',expected', ...
    'max_error',max(abs(actual-expected)),'result','confirmed second imported component omitted');
fid = fopen(fullfile(out,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(r,PrettyPrint=true));
fclose(fid);
fprintf('Confirmed: second imported component omitted (max error32).\n');
end

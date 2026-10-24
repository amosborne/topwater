%[text] # Imported excitation history reproduction
%[text] Independent indexed convolution checks physical samples beyond the requested output interval.
%%
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
out = fullfile(root,'artifacts/reference-imported-excitation-stock');
assert(~isfolder(out),'Use a fresh evidence directory');
mkdir(out);
t = (-2:12)';
to = (0:10)';
tau = -2:2;
k = [1 2 4 8 16];
for j = 1:3
    impulse = [11 -1 5];
    eta = double(t == impulse(j));
    file = fullfile(out,sprintf('eta_%d.txt',j));
    writematrix([t eta],file,delimiter='tab');
    w = waveClass('elevationImport');
    w.elevationFile = file;
    w.direction = 0;
    w.setup([0.1 12],inf,0,1,10,to',1,1);
    b = paddingProbe();
    actual = b.probe(w);
    expected = zeros(size(to));
    for it = 1:numel(to)
        for ik = 1:numel(tau)
            expected(it) = expected(it) + k(ik)*double(to(it)-tau(ik) == impulse(j));
        end
    end
    res(j).impulse_time = impulse(j);
    res(j).actual = actual';
    res(j).expected = expected';
    res(j).max_error = max(abs(actual-expected));
end
assert(res(1).max_error == 2 && res(2).max_error == 16 && res(3).max_error == 0, ...
    'Stock behavior differs from the independently predicted cropping defect');
fid = fopen(fullfile(out,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('dependency','a23d2993807410fe78cc63546c1ba9d81912024f', ...
    'result','confirmed crop defect; interior indexing agrees','checks',res),PrettyPrint=true));
fclose(fid);
fprintf('Confirmed: imported future/past samples dropped; interior convolution agrees.\n');

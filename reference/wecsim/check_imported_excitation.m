function check_imported_excitation()
clear classes
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
addpath(fullfile(root,'reference/wecsim'));
assert(contains(which('bodyClass'),'topwater-wec-sim'));
out = fullfile(root,'artifacts/reference-imported-excitation-fixed-attempt2');
assert(~isfolder(out));
mkdir(out);
n = 0;
for dt = [1 0.5]
    t = (-3:dt:13)';
    to = (0:dt:10)';
    tau = -2:dt:2;
    k = interp1(-2:2,[1 2 4 8 16],tau);
    for ti = [11 -1 5]
        n = n+1;
        eta = double(t == ti);
        file = fullfile(out,sprintf('wave_%d.txt',n));
        writematrix([t eta],file,delimiter='tab');
        w = waveClass('elevationImport');
        w.elevationFile = file;
        w.direction = 0;
        w.setup([0.1 12],inf,0,dt,numel(to)-1,to',1,1);
        before = w.waveAmpTime;
        b = paddingProbe();
        actual = b.probeExtended(w,dt,0);
        expected = zeros(size(to));
        for it = 1:numel(to)
            for ik = 1:numel(tau)
                expected(it) = expected(it) + dt*k(ik)*double(to(it)-tau(ik) == ti);
            end
        end
        assert(max(abs(actual-expected)) < 1e-12);
        assert(isequal(before,w.waveAmpTime),'Output elevation grid changed');
        res(n).dt = dt;
        res(n).impulse_time = ti;
        res(n).max_error = max(abs(actual-expected));
    end
end
% insufficient future support must be explicit rather than silently zero filled
w.elevationFile = fullfile(out,'short.txt');
writematrix([to zeros(size(to))],w.elevationFile,delimiter='tab');
rejected = false;
try
    b.probeExtended(w,dt,0);
catch e
    rejected = contains(e.message,'past/future support');
end
assert(rejected);
% nonzero legacy ramp retains its output-grid elevation and known zero prehistory
w.elevationFile = fullfile(out,'ramp.txt');
writematrix([t ones(size(t))],w.elevationFile,delimiter='tab');
w.setup([0.1 12],inf,2,dt,numel(to)-1,to',1,1);
before = w.waveAmpTime;
actual = b.probeExtended(w,dt,2);
expected = zeros(size(to));
for it = 1:numel(to)
    for ik = 1:numel(tau)
        tq = to(it)-tau(ik);
        ramp = (1-cos(pi*tq/2))/2;
        if tq <= 0, ramp = 0; end
        if tq >= 2, ramp = 1; end
        expected(it) = expected(it)+dt*k(ik)*ramp;
    end
end
assert(max(abs(actual-expected)) < 1e-12);
assert(isequal(before,w.waveAmpTime));
% unequal headings, unequal kernels and every output DOF through production wiring
t = (-3:0.5:13)';
to = (0:0.5:10)';
tau = -2:0.5:2;
base = [1 2 4 8 16];
other = [-3 1 2 -2 5];
kernel = zeros(6,2,5);
for q = 1:6
    kernel(q,1,:) = q*base;
    kernel(q,2,:) = (7-q)*other;
end
b.configure([0 90],-2:2,kernel);
for j = 1:2
    wm(j) = waveClass('elevationImport');
    wm(j).direction = 90*(j-1);
    wm(j).elevationFile = fullfile(out,sprintf('direction_%d.txt',j));
    writematrix([t j*double(t == 5)],wm(j).elevationFile,delimiter='tab');
    wm(j).setup([0.1 12],inf,0,0.5,numel(to)-1,to',1,1);
end
s = simulationClass();
s.dt = 0.5;
s.rho = 1;
s.gravity = 1;
s.rampTime = 0;
actual = b.probeMultipleAll(wm,s);
expected = zeros(numel(to),6);
for it = 1:numel(to)
    for ik = 1:numel(tau)
        for q = 1:6
            expected(it,q) = expected(it,q)+0.5*(q*interp1(-2:2,base,tau(ik)) ...
                +2*(7-q)*interp1(-2:2,other,tau(ik)))*double(to(it)-tau(ik) == 5);
        end
    end
end
multiError = max(abs(actual-expected),[],'all');
assert(multiError < 1e-12);
% asymmetric physical lag support verifies zero lag is time based, not array centre
b.configure(0,-1:2,repmat(reshape([2 4 8 16],1,1,4),6,1,1));
actual = b.probeExtended(wm(1),0.5,0);
expected = zeros(numel(to),1);
ta = -1:0.5:2;
for it = 1:numel(to)
    for ik = 1:numel(ta)
        expected(it) = expected(it)+0.5*interp1(-1:2,[2 4 8 16],ta(ik))*double(to(it)-ta(ik) == 5);
    end
end
asymmetricError = max(abs(actual-expected));
assert(asymmetricError < 1e-12);
fid = fopen(fullfile(out,'result.json'),'w');
fprintf(fid,'%s\n',jsonencode(struct('checks',res,'short_input_rejected',rejected, ...
    'multiple_direction_all_dof_error',multiError,'asymmetric_error',asymmetricError, ...
    'source',which('bodyClass')),PrettyPrint=true));
fclose(fid);
fprintf('Imported excitation retained-history regressions passed.\n');
end

function evidence = verify_seeded_cic_output(caseDir,seed,hf,cicTime,time,velocity,forceBus,timeLog,forceLog,radLog)
dt = cicTime(2)-cicTime(1);
N = numel(cicTime)-1;
Nt = numel(time);
weights = dt*ones(N+1,1);
weights([1 end]) = dt/2;
fullVelocity = [seed.history_velocity_world(1:end-1,:);velocity];
fullForce = zeros(size(fullVelocity));
for i = 1:6
    for j = 1:6
        values = conv(fullVelocity(:,j),weights.*hf.irkb(:,i,j));
        fullForce(:,i) = fullForce(:,i)+values(1:size(fullVelocity,1));
    end
end
assert(numel(timeLog.time) == Nt && numel(forceLog.time) == Nt && numel(radLog.time) == Nt);
assert(max(abs(timeLog.time(:)-time(:))) < 1e-9);
assert(max(abs(forceLog.time(:)-time(:))) < 1e-9 && max(abs(radLog.time(:)-time(:))) < 1e-9);
rowsError = zeros(Nt,6);
radiationError = zeros(Nt,6);
busError = zeros(Nt,6);
expectedRadiation = zeros(Nt,6);
for sample = 1:Nt
    times = read_sample(timeLog,sample,[3 1]);
    actualRows = read_sample(forceLog,sample,[3 6]);
    index = round(times/dt)+N+4;
    assert(all(index >= N+1 & index <= size(fullForce,1)));
    expectedRows = fullForce(index,:);
    rowsError(sample,:) = max(abs(actualRows-expectedRows),[],1);
    expectedRadiation(sample,:) = interp1(times,expectedRows,time(sample),'pchip','extrap');
    actualRadiation = read_sample(radLog,sample,[1 6]);
    radiationError(sample,:) = actualRadiation-expectedRadiation(sample,:);
    busError(sample,:) = forceBus(sample,:)-actualRadiation;
end
tolerance = 1e-6+1e-10*max(abs(expectedRadiation),[],1);
save(fullfile(caseDir,'seeded-cic-verification.mat'),'rowsError','radiationError','busError','expectedRadiation','tolerance','-v7');
assert(all(max(rowsError,[],1) <= tolerance));
assert(all(max(abs(radiationError),[],1) <= tolerance));
assert(all(max(abs(busError),[],1) <= tolerance));
evidence = struct('force_history_row_max_error',max(rowsError,[],1), ...
    'native_pchip_radiation_max_error',max(abs(radiationError),[],1), ...
    'native_radiation_bus_max_error',max(abs(busError),[],1), ...
    'initial_expected_radiation',expectedRadiation(1,:),'tolerance',tolerance, ...
    'scope','pre0specifiedhistory plus ACTUALpost0worldvelocities; no prescribedmotion/Archange');
end

function value = read_sample(log,k,shape)
data = log.signals.values;
if size(data,1) == numel(log.time)
    value = reshape(data(k,:,:),shape);
else
    assert(size(data,3) == numel(log.time));
    value = reshape(data(:,:,k),shape);
end
end

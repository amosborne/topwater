function correction = cic_seed_table(irkb,cicTime,historyVelocity)
dt = cicTime(2)-cicTime(1);
N = numel(cicTime)-1;
assert(size(irkb,1) == N+1 && size(irkb,2) == 6 && size(irkb,3) == 6);
assert(isequal(size(historyVelocity),[N+4 6]));
assert(max(abs(cicTime(:)-(0:N)'*dt)) < 1e-10);
weights = dt*ones(N+1,1);
weights([1 end]) = dt/2;
correction = zeros(N+4,6);
% History rows cover integer times -N-3:0; correction rows cover -3:N.
for n = -3:N-1
    if n < 0
        lag = (0:N)';
    else
        lag = (n+1:N)';
    end
    velocity = historyVelocity(n-lag+N+4,:);
    integrand = sum(irkb(lag+1,:,:).*reshape(velocity,[numel(lag) 1 6]),3);
    correction(n+4,:) = weights(lag+1)'*integrand;
end
end

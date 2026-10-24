function combined = cic_seed_lookup(timeExtrap,FradExtrap,correction,dt)
%#codegen
combined = FradExtrap;
for row = 1:3
    n = round(timeExtrap(row)/dt);
    assert(abs(timeExtrap(row)-n*dt) < 1e-8);
    assert(n >= -3);
    if n+4 <= size(correction,1)
        combined(row,:) = combined(row,:)+correction(n+4,:);
    end
end
end

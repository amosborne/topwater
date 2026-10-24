function errorMax = verify_reference_extension_cic(native_body,hf,dt,caseDir,timeLog,forceLog,radLog)
N = size(hf.irkb,1)-1;
assert(N == 500 && dt == .02);
% Reuse only the independently checked verifier, with IDENTICALLY ZERO pre0 v.
% No seed initializer/history adapter is called or inserted into the model.
zeroPast.history_velocity_world = zeros(N+4,6);
evidence = verify_seeded_cic_output(caseDir,zeroPast,hf,(0:N)*dt, ...
    native_body.time,native_body.signals.values(:,7:12), ...
    native_body.signals.values(:,31:36),timeLog,forceLog,radLog);
errorMax = max([evidence.force_history_row_max_error ...
    evidence.native_pchip_radiation_max_error evidence.native_radiation_bus_max_error]);
evidence.scope = 'zero pre0 bodyvelocity; ACTUALpost0worldv; actual3rows/pairedRT/PCHIP/bus; no history seed';
save(fullfile(caseDir,'zero-history-cic-check.mat'),'evidence','errorMax','-v7');
end

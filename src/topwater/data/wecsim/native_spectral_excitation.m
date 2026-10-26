function excitation = native_spectral_excitation(body,waves,time)
% Evaluate the actual WEC-Sim FullDir function with its native coefficient cache.
assert(isscalar(waves) && strcmp(waves.type,'spectrumImportFullDir'), ...
    'Expected one native FullDir wave object');
assert(body.dof == 6 && body.yaw.option == 0 && body.largeXYDisplacement.option == 0, ...
    'Spectral workflow requires fixed linear six-DOF excitation');
data = getExcitationData(body.hydroData(1),waves.type);
assert(data.uses_component_view,'Native spectral component view is missing');
assert(all(ismember(waves.omega(:),data.w(:))) && ...
    all(ismember(waves.fullDirectionalSpectrum.directions(:),data.direction(:))), ...
    'Native spectral queries differ from the prepared component nodes');
h = body.hydroForce.hf1.fExt;
directions = waves.fullDirectionalSpectrum.directions;
nDirection = numel(directions);
nFrequency = numel(waves.omega);
previousYaw = zeros(size(directions));
previousMD = zeros(nDirection,nFrequency,6);
previousRE = zeros(nDirection,nFrequency,6);
previousIM = zeros(nDirection,nFrequency,6);
excitation = zeros(numel(time),6);
for k = 1:numel(time)
    [excitation(k,:),previousYaw,previousMD,previousRE,previousIM] = ...
        irregExcFullDirF(0,waves.amplitude,waves.omega, ...
        h.dofGrd,h.dirGrd,h.wGrd,h.qDofGrd,h.qWGrd,h.fEHRE,h.fEHIM,h.fEHMD, ...
        waves.phase,waves.dOmega,time(k),directions,zeros(1,6),body.yaw.threshold, ...
        previousYaw,previousMD,previousRE,previousIM);
end
assert(all(isfinite(excitation),'all'),'Native spectral excitation is nonfinite');
end

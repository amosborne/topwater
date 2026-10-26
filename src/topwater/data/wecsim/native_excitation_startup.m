function report = native_excitation_startup(prepared,excitation)
% Diagnostic scale only; WEC-Sim's own imported excitation is the applied load.
n = prepared.n;
data = load(prepared.native_mat_path);
w = data.excitation_omega_rad_s(:);
theta = data.excitation_headings_rad(:);
H = data.excitation_native;
omega = prepared.carrier_omega_rad_s(:);
direction = mod(prepared.carrier_directions_rad(:),2*pi);
amplitude = prepared.carrier_amplitudes_m(:);
assert(isequal(size(H),[numel(w) numel(theta) 6]) && all(isfinite(H),'all'));
assert(numel(omega) == numel(direction) && numel(omega) == numel(amplitude));
scale = zeros(1,6);
for j = 1:6
    value = interp2(theta,w,H(:,:,j),direction,omega,'linear');
    assert(all(isfinite(value)),'Carrier scale lies outside the native coefficient coverage');
    scale(j) = sum(amplitude.*abs(value));
end
first = excitation(1,:);
weight = n.M(1,1)*n.g;
length = prepared.boat_length_m;
report = struct('initial_excitation',first,'carrier_scale',scale, ...
    'force_norm_limit',.001*weight,'moment_norm_limit',.001*weight*length, ...
    'per_axis_limit',.001*scale,'ramp_time_s',0);
save('native-startup.mat','report','-v7');
assert(norm(first(1:3)) <= report.force_norm_limit && ...
    norm(first(4:6)) <= report.moment_norm_limit && all(abs(first) <= report.per_axis_limit), ...
    'WEC-Sim native initial excitation exceeds the unchanged engineering calm limits');
end

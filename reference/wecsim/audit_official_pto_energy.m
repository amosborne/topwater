function energy = audit_official_pto_energy(bodies,ptos,storage,parameters,settings)
root = 'C:/Users/alex/topwater';
baseline = load(fullfile(root,'artifacts/reference/official-examples-20261024-ready2/cic/official/tests/RegressionTests/RegularWaves/regularCIC/audit.mat'));
assert(settings.pto_damping == 0 && baseline.settings.pto_damping == 1200000);
assert(settings.dt_s == baseline.settings.dt_s && settings.cic_s == baseline.settings.cic_s);
assert(settings.ramp_s == baseline.settings.ramp_s && settings.period_s == baseline.settings.period_s);
assert(settings.height_m == baseline.settings.height_m && settings.end_s == baseline.settings.end_s);
energy.modified = budget(bodies,ptos,storage,parameters,settings.pto_damping);
energy.baseline = budget(baseline.bodies,baseline.ptos,baseline.storage,parameters,baseline.settings.pto_damping);
energy.scope = 'actual heave-only two-body native mechanical work budget; radiation work is not farfield dissipated energy';
energy.mass_accounting = 'E_S includes transferred addedmass once; W_Ar uses stored APPLIED residual force, not postprocessed full Ainf';
energy.comparison = 'same official models/H5/inertia/ramps/CIC/dt; sole physical change PTO damping1200000 to0';
end

function result = budget(bodies,ptos,storage,parameters,c)
t = bodies(1).time(:);
assert(numel(bodies) == 2 && isequal(t,bodies(2).time(:)));
Nt = numel(t);
Tdry = zeros(Nt,1);
Ttransfer = zeros(Nt,1);
U = zeros(Nt,1);
powerExc = zeros(Nt,1);
powerRad = zeros(Nt,1);
powerAr = zeros(Nt,1);
powerOther = zeros(Nt,1);
powerSolid = zeros(Nt,1);
powerHydro = zeros(Nt,1);
for j = 1:2
    hf = parameters{j};
    b = bodies(j);
    v = b.velocity;
    a = b.acceleration;
    assert(all(isfinite(v),'all') && all(isfinite(a),'all'));
    assert(max(abs(v(:,[1 2 4 5 6])),[],'all') < 1e-8,'Budget requires actual heave-only kinematics');
    cg = b.centerGravity(:)';
    z = b.position(:,3)-cg(3);
    Tdry = Tdry+.5*hf.mass*v(:,3).^2;
    Ttransfer = Ttransfer+.5*(hf.adjustedMass-hf.mass)*v(:,3).^2;
    U = U+.5*hf.linearHydroRestCoef(3,3)*z.^2;
    powerExc = powerExc+sum(v.*b.forceExcitation,2);
    powerRad = powerRad+sum(v.*b.forceRadiationDamping,2);
    powerAr = powerAr+sum(v.*storage{j}.output_forceAddedMass,2);
    powerOther = powerOther+sum(v.*(b.forceMorisonAndViscous+b.forceLinearDamping),2);
    powerSolid = powerSolid+hf.adjustedMass*v(:,3).*a(:,3);
    powerHydro = powerHydro+sum(v.*storage{j}.output_forceTotal,2);
end
relativeV = bodies(1).velocity(:,3)-bodies(2).velocity(:,3);
powerPTO = c*relativeV.^2;
nativePTO = sum(ptos(1).powerInternalMechanics,2);
Esolid = Tdry+Ttransfer+U;
Wexc = cumtrapz(t,powerExc);
Wrad = cumtrapz(t,powerRad);
War = cumtrapz(t,powerAr);
Wother = cumtrapz(t,powerOther);
Wpto = cumtrapz(t,powerPTO);
residual = Esolid-Esolid(1)-(Wexc-Wrad-War-Wother-Wpto);
pointwiseResidual = powerSolid-powerHydro+powerPTO;
result = struct('time',t,'dry_kinetic_J',Tdry,'transferred_added_kinetic_J',Ttransfer, ...
    'hydrostatic_J',U,'Solid_plus_hydrostatic_J',Esolid,'incident_work_J',Wexc, ...
    'radiation_resisting_work_J',Wrad,'residual_added_work_J',War, ...
    'other_damping_work_J',Wother,'PTO_absorbed_work_J',Wpto, ...
    'native_PTO_signed_power_W',nativePTO,'expected_PTO_absorbed_power_W',powerPTO, ...
    'native_PTO_plus_expected_max_W',max(abs(nativePTO+powerPTO)), ...
    'native_PTO_minus_expected_max_W',max(abs(nativePTO-powerPTO)), ...
    'native_PTO_magnitude_error_W',max(abs(abs(nativePTO)-powerPTO)), ...
    'mechanical_energy_residual_J',residual,'pointwise_Solid_PTO_power_residual_W',pointwiseResidual, ...
    'maximum_abs_energy_residual_J',max(abs(residual)), ...
    'maximum_abs_pointwise_power_residual_W',max(abs(pointwiseResidual)), ...
    'scale_max_abs_incident_work_J',max(abs(Wexc)),'peak_mechanical_energy_J',max(Esolid), ...
    'radiation_work_minimum_J',min(Wrad),'radiation_power_minimum_W',min(powerRad));
assert(all(isfinite(residual)) && all(isfinite(pointwiseResidual)));
if c == 0
    assert(all(powerPTO == 0) && all(nativePTO == 0),'PTO zero-damping consumption not zero');
end
end

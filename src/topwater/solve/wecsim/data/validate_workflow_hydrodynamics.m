function validate_workflow_hydrodynamics(h,stage,rawReal,rawImaginary)
% Check only arrays consumed by this single-body BEMIO workflow, never water depth.
if strcmp(stage,'spectral_input')
    assert(nargin == 2,'Standard struct validation takes no raw excitation arrays');
    validate_spectral_input(h);
    return
end
fields = {'rho','g','w','T','theta','Nf','Nh','A','B', ...
    'fk_re','fk_im','fk_ma','fk_ph','sc_re','sc_im','sc_ma','sc_ph', ...
    'ex_re','ex_im','ex_ma','ex_ph'};
if strcmp(stage,'prepared') || strcmp(stage,'spectral')
    assert(nargin == 2,[stage ' validation takes no raw excitation arrays']);
    fields = [fields {'cg','cb','Vo','Khs','Ainf','ra_K','ra_t','ra_w'}];
    if strcmp(stage,'prepared')
        fields = [fields {'ex_K','ex_t','ex_w'}];
    end
elseif strcmp(stage,'imported')
    assert(nargin == 4,'Imported validation requires raw real and imaginary excitation');
else
    error('Topwater:InvalidHydrodynamics','Unknown hydrodynamics validation stage');
end
for k = 1:numel(fields)
    name = fields{k};
    if ~isfield(h,name)
        error('Topwater:InvalidHydrodynamics','%s hydrodynamics is missing %s',stage,name);
    end
    require_finite(h.(name),[stage '.' name]);
end
require_grid(h.w,[stage '.w'],2,true);
require_grid(h.theta,[stage '.theta'],1,false);
require_shape(h.Nf,[1 1],[stage '.Nf'],false);
require_shape(h.Nh,[1 1],[stage '.Nh'],false);
if h.Nf ~= numel(h.w) || h.Nh ~= numel(h.theta)
    error('Topwater:InvalidHydrodynamics','%s frequency/heading counts do not match their reference grids',stage);
end
require_shape(h.T,[1 numel(h.w)],[stage '.T'],true);
require_shape(h.A,[6 6 numel(h.w)],[stage '.A'],false);
require_shape(h.B,[6 6 numel(h.w)],[stage '.B'],false);
excitationShape = [6 numel(h.theta) numel(h.w)];
for name = {'fk_re','fk_im','fk_ma','fk_ph','sc_re','sc_im','sc_ma','sc_ph', ...
        'ex_re','ex_im','ex_ma','ex_ph'}
    require_shape(h.(name{1}),excitationShape,[stage '.' name{1}],false);
end
require_shape(h.rho,[1 1],[stage '.rho'],false);
require_shape(h.g,[1 1],[stage '.g'],false);
if strcmp(stage,'imported')
    require_finite(rawReal,'netCDF real excitation');
    require_finite(rawImaginary,'netCDF imaginary excitation');
    require_shape(rawReal,excitationShape,'netCDF real excitation',false);
    require_shape(rawImaginary,excitationShape,'netCDF imaginary excitation',false);
else
    require_shape(h.cg,[1 3],[stage '.cg'],true);
    require_shape(h.cb,[1 3],[stage '.cb'],true);
    require_shape(h.Vo,[1 1],[stage '.Vo'],false);
    require_shape(h.Khs,[6 6],[stage '.Khs'],false);
    require_shape(h.Ainf,[6 6],[stage '.Ainf'],false);
    require_grid(h.ra_t,[stage '.ra_t'],2,false);
    require_grid(h.ra_w,[stage '.ra_w'],2,true);
    require_shape(h.ra_K,[6 6 numel(h.ra_t)],[stage '.ra_K'],false);
    if strcmp(stage,'prepared')
        require_grid(h.ex_t,'prepared.ex_t',2,false);
        require_grid(h.ex_w,'prepared.ex_w',2,false);
        if any(h.ex_w < 0)
            error('Topwater:InvalidHydrodynamics','prepared.ex_w must be nonnegative');
        end
        require_shape(h.ex_K,[6 numel(h.theta) numel(h.ex_t)],'prepared.ex_K',false);
    end
end
end

function validate_spectral_input(h)
if ~isstruct(h) || ~isscalar(h) || ~isfield(h,'simulation_parameters') ...
        || ~isfield(h,'hydro_coeffs') || ~isfield(h.hydro_coeffs,'excitation') ...
        || ~isfield(h.hydro_coeffs,'mean_drift') ...
        || ~all(isfield(h.simulation_parameters,{'w','direction'})) ...
        || ~all(isfield(h.hydro_coeffs.excitation,{'re','im'}))
    error('Topwater:InvalidHydrodynamics','Standard spectral hydro struct is incomplete');
end
w = h.simulation_parameters.w;
theta = h.simulation_parameters.direction;
require_finite(w,'spectral_input.w');
require_finite(theta,'spectral_input.direction');
require_grid(w,'spectral_input.w',3,false);
require_grid(theta,'spectral_input.direction',2,false);
if w(1) ~= 0 || theta(1) ~= -180 || theta(end) ~= 180
    error('Topwater:InvalidHydrodynamics','Standard coefficient input needs zero and exact -180/180 boundary nodes');
end
shape = [6 numel(theta) numel(w)];
values = {h.hydro_coeffs.excitation.re,h.hydro_coeffs.excitation.im,h.hydro_coeffs.mean_drift};
for k = 1:numel(values)
    value = values{k};
    require_finite(value,'spectral_input.coefficients');
    require_shape(value,shape,'spectral_input.coefficients',false);
    if ~isequal(value(:,1,:),value(:,end,:))
        error('Topwater:InvalidHydrodynamics','Standard coefficient periodic boundary values differ');
    end
end
if any(h.hydro_coeffs.mean_drift ~= 0,'all')
    error('Topwater:InvalidHydrodynamics','Standard spectral input mean drift must be explicitly disabled/zero');
end
end

function require_finite(value,name)
if ~isnumeric(value) || ~isreal(value) || isempty(value)
    error('Topwater:InvalidHydrodynamics','%s must be a nonempty real numeric array',name);
end
if ~all(isfinite(value),'all')
    error('Topwater:NonfiniteHydrodynamics','%s contains nonfinite values',name);
end
end

function require_grid(value,name,minimumCount,positive)
if ~isvector(value) || numel(value) < minimumCount || any(diff(value(:)) <= 0) ...
        || (positive && any(value(:) <= 0))
    error('Topwater:InvalidHydrodynamics','%s must be an ordered compatible reference grid',name);
end
end

function require_shape(value,expected,name,vector)
if vector
    matches = isvector(value) && numel(value) == expected(end);
else
    actual = arrayfun(@(dimension) size(value,dimension),1:numel(expected));
    matches = ndims(value) <= numel(expected) && isequal(actual,expected);
end
if ~matches
    error('Topwater:InvalidHydrodynamics','%s has a shape incompatible with this six-mode workflow',name);
end
end

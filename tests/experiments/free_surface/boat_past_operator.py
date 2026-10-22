"""Past-only two-direction linear operator for bounded offline/coupled probes.

One coefficient preserves the CURRENT odd pressure chord. Only previous
operators fit the transverse coefficient. A same-data ordinary LS fit is a
negative control. Neither model promises a pressure-error bound.
"""
import numpy as np

from .boat_operator_reuse import area_weights, _BASIS, InvalidReuse


LIMITS = dict(maximum_age_seconds=.051, maximum_normalized_coordinate=80.,
    maximum_body_residual_rms_m=.00025, maximum_free_residual_rms_m=.00025,
    maximum_direction_residual_fraction=.05, refresh_age_seconds=.04)


def coordinates(data, points):
    change = points-data['past_origin']
    directions = data['past_directions']
    weight = data['past_metric']
    alpha = np.linalg.solve(data['past_gram'],
        [np.sum(weight[:, None]*direction*change) for direction in directions])
    residual = change-np.einsum('i,ijk->jk', alpha, directions)
    diagnostics = dict(coefficients=alpha.tolist())
    for label in ('body', 'free'):
        diagnostics[label+'_residual_rms_m'] = float(np.sqrt(np.sum(
            data['past_'+label+'_weight'][:, None]*residual**2)))
    return alpha, diagnostics


def predict(data, arrays, *, ordinary=False):
    alpha, diagnostics = coordinates(data, arrays['points'])
    key = 'past_ordinary_coefficients' if ordinary else 'past_constrained_coefficients'
    return data['past_energy']+np.einsum('i,ijk->jk', alpha, data[key]), diagnostics


def direction_diagnostics(data, plus, minus):
    change = (plus['points']-minus['points'])/2
    alpha = np.linalg.solve(data['past_gram'],
        [np.sum(data['past_metric'][:, None]*g*change) for g in data['past_directions']])
    residual = change-np.einsum('i,ijk->jk', alpha, data['past_directions'])
    result = dict(coefficients=alpha.tolist())
    for label in ('body', 'free'):
        weight = data['past_'+label+'_weight'][:, None]
        result[label+'_out_of_subspace_fraction'] = float(np.sqrt(
            np.sum(weight*residual**2)/max(np.sum(weight*change**2), 1e-30)))
    return result


def fit(samples, energies, *, ordinary=True):
    """Six records: previous base/plus/minus, then current base/plus/minus."""
    previous, _, _, base, plus, minus = samples
    for arrays in samples:
        if any(not np.array_equal(arrays[key], base[key]) for key in _BASIS):
            raise ValueError('past fit requires identical chart basis and ownership')
    weights = area_weights(base)
    metric = .5*(weights['body']+weights['free'])

    def inner(a, b):
        return float(np.sum(metric[:, None]*a*b))

    first = (plus['points']-minus['points'])/2
    first_norm2 = inner(first, first)
    if first_norm2 <= 1e-24:
        raise ValueError('current pressure chord is numerically stationary')
    change = previous['points']-base['points']
    second = change-first*inner(first, change)/first_norm2
    raw_second_norm = float(np.sqrt(inner(second, second)))
    first_norm = float(np.sqrt(first_norm2))
    change_norm = float(np.sqrt(inner(change, change)))
    if raw_second_norm <= max(1e-12, 128*np.finfo(float).eps*(change_norm+first_norm)):
        raise ValueError('past displacement gives no independent second direction')
    amplification = first_norm/raw_second_norm
    second *= amplification
    directions = np.array([first, second])
    gram = np.array([[inner(a, b) for b in directions] for a in directions])
    condition = float(np.linalg.cond(gram))
    if condition*np.finfo(float).eps > 1e-8:
        raise ValueError('past geometry metric exceeds the conditioning budget')
    data = dict(past_origin=base['points'], past_metric=metric, past_directions=directions,
        past_body_weight=weights['body'], past_free_weight=weights['free'],
        past_gram=gram, past_energy=energies[3])
    features = np.array([coordinates(data, a['points'])[0] for a in samples])
    rhs = np.stack([energy-energies[3] for energy in energies])
    first_coefficient = (energies[4]-energies[5])/2
    transverse = features[:3, 1]
    denominator = transverse@transverse
    if not np.isfinite(denominator) or denominator <= 1e-24:
        raise ValueError('past transverse coefficient is unobservable')
    second_coefficient = np.einsum('i,ijk->jk', transverse,
        rhs[:3]-features[:3, 0, None, None]*first_coefficient)/denominator
    constrained = np.array([first_coefficient, second_coefficient])
    if not np.isfinite(constrained).all():
        raise ValueError('nonfinite past fit')
    data['past_constrained_coefficients'] = constrained
    description = dict(complete=True, scope='Past-only fixed model; current center and odd operator chord constrained, pressure need not match. No future geometry or operator informs fit.',
        raw_current_chord_rms_m=first_norm, raw_transverse_rms_m=raw_second_norm,
        raw_transverse_to_past_displacement_ratio=raw_second_norm/change_norm,
        transverse_normalization_factor=amplification, gram_condition=condition,
        sample_coordinates=features.tolist())
    fits = dict(constrained=constrained)
    if ordinary:
        coefficients, _, rank, singular = np.linalg.lstsq(features, rhs.reshape(6, -1), rcond=None)
        feature_condition = float(singular[0]/singular[-1])
        if rank != 2 or feature_condition*np.finfo(float).eps > 1e-8 or not np.isfinite(coefficients).all():
            raise ValueError('ordinary past fit exceeds the conditioning budget')
        data['past_ordinary_coefficients'] = coefficients.reshape(2, *energies[3].shape)
        fits['ordinary'] = data['past_ordinary_coefficients']
        description.update(ordinary_feature_singular_values=singular.tolist(), ordinary_feature_condition=feature_condition)
    for label, coefficients in fits.items():
        error = rhs-np.einsum('ij,jkl->ikl', features, coefficients)
        description[label+'_sample_relative_errors'] = [float(np.linalg.norm(e)/max(np.linalg.norm(r), 1e-30))
                                                        for e, r in zip(error, rhs)]
        chord_error = (features[4]-features[5])@coefficients.reshape(2, -1)-(energies[4]-energies[5]).ravel()
        description[label+'_current_odd_chord_relative_error'] = float(np.linalg.norm(chord_error)
            /max(np.linalg.norm(energies[4]-energies[5]), 1e-30))
    return data, description


class PastLocalOperator:
    """Bounded coupled falsifier; geometric limits still do not certify forces."""
    def __init__(self, previous, current, time, rhs_index):
        self.data, self.diagnostics = fit([*previous[0], *current[0]], [*previous[1], *current[1]],
                                          ordinary=False)
        self.basis = {key: current[0][0][key].copy() for key in _BASIS}
        self.time, self.rhs_index = float(time), rhs_index
        self.diagnostics.update(anchor_time_seconds=self.time, anchor_rhs_index=rhs_index,
            previous_time_seconds=previous[2], previous_rhs_index=previous[3])

    def _same_basis(self, arrays):
        if (arrays['points'].shape != self.data['past_origin'].shape
                or any(not np.array_equal(arrays[key], value) for key, value in self.basis.items())):
            raise InvalidReuse('chart basis or boundary ownership changed')

    def predict(self, arrays, time):
        self._same_basis(arrays)
        age = float(time-self.time)
        if abs(age) > LIMITS['maximum_age_seconds']+1e-14:
            raise InvalidReuse('operator age limit')
        alpha, diagnostics = coordinates(self.data, arrays['points'])
        if not np.isfinite(alpha).all() or np.max(abs(alpha)) > LIMITS['maximum_normalized_coordinate']:
            raise InvalidReuse('operator extrapolation limit')
        for label in ('body', 'free'):
            error = diagnostics[label+'_residual_rms_m']
            if not np.isfinite(error) or error > LIMITS['maximum_'+label+'_residual_rms_m']:
                raise InvalidReuse(label+' geometry leaves the local directions')
        energy = self.data['past_energy']+np.einsum('i,ijk->jk', alpha,
                                                  self.data['past_constrained_coefficients'])
        diagnostics.update(age_seconds=age, anchor_rhs_index=self.rhs_index)
        return energy, diagnostics

    def check_direction(self, plus, minus):
        self._same_basis(plus)
        self._same_basis(minus)
        diagnostics = direction_diagnostics(self.data, plus, minus)
        for label in ('body', 'free'):
            error = diagnostics[label+'_out_of_subspace_fraction']
            if not np.isfinite(error) or error > LIMITS['maximum_direction_residual_fraction']:
                raise InvalidReuse(label+' pressure direction leaves the local directions')
        return diagnostics

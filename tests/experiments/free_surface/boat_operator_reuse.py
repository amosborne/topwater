"""Past-only, one-direction quadratic reuse for the bounded boat experiment.

These predeclared geometry limits are falsifiable heuristics, not force-error
bounds. A caller must use one model for an entire RHS pressure stencil, or
discard that attempt and rebuild the entire stencil. No fluid state is cached.
"""
import numpy as np

from topwater.rankine import TrianglePanels


class InvalidReuse(ValueError):
    pass


LIMITS = dict(maximum_age_seconds=.01, maximum_normalized_coordinate=12.,
              maximum_body_residual_rms_m=.00002, maximum_free_residual_rms_m=.00002,
              maximum_direction_residual_fraction=.05)
_BASIS = ('faces', 'body_faces', 'free_faces', 'free_nodes', 'body_source_group',
          'free_source_parent', 'free_parent_indices')


def area_weights(arrays):
    """Separate normalized body/free lumped-area metrics at one fixed geometry."""
    result = {}
    for label in ('body', 'free'):
        faces = arrays['faces'][arrays[label+'_faces']]
        area = TrianglePanels(arrays['points'][faces], 1).areas
        weight = np.bincount(faces.ravel(), weights=np.repeat(area/3, 3),
                             minlength=len(arrays['points']))
        result[label] = weight/weight.sum()
    return result


class LocalOperator:
    def __init__(self, arrays, energies, time, rhs_index):
        """Train from fresh [base, plus, minus] of one completed physical stencil."""
        base, plus, minus = arrays
        self.basis = {key: base[key].copy() for key in _BASIS}
        self.origin = base['points'].copy()
        self.time, self.rhs_index = float(time), rhs_index
        self.weight = area_weights(base)
        # Equal body/free weights prevent a large outer surface hiding hull error.
        self.metric = .5*(self.weight['body']+self.weight['free'])
        for value in (plus, minus):
            self._same_basis(value)
        self.direction = (plus['points']-minus['points'])/2
        self.norm2 = self._inner(self.direction, self.direction)
        if not np.isfinite(self.norm2) or self.norm2 <= (1e-12)**2:
            raise InvalidReuse('training geometry motion is too small')
        coordinates, residuals = [], []
        for value in (plus, minus):
            z, diagnostics = self.coordinates(value)
            coordinates.append([z, .5*z*z])
            residuals.append(diagnostics)
        features = np.array(coordinates)
        condition = float(np.linalg.cond(features))
        if not np.isfinite(condition) or condition*np.finfo(float).eps > 1e-8:
            raise InvalidReuse('training coordinates exceed the conditioning budget')
        self.energy = energies[0]
        self.coefficients = np.linalg.solve(features,
            np.stack([energy-self.energy for energy in energies[1:]]).reshape(2, -1)
            ).reshape(2, *self.energy.shape)
        if not np.isfinite(self.coefficients).all():
            raise InvalidReuse('nonfinite local operator fit')
        self.diagnostics = dict(anchor_time_seconds=self.time, anchor_rhs_index=rhs_index,
            training_feature_condition=condition, training_coordinates=coordinates,
            training_geometry=residuals)

    def _inner(self, a, b):
        return float(np.sum(self.metric[:, None]*a*b))

    def _same_basis(self, arrays):
        if (arrays['points'].shape != self.origin.shape
                or any(not np.array_equal(arrays[key], value) for key, value in self.basis.items())):
            raise InvalidReuse('chart basis or boundary ownership changed')

    def coordinates(self, arrays):
        self._same_basis(arrays)
        change = arrays['points']-self.origin
        z = self._inner(self.direction, change)/self.norm2
        residual = change-z*self.direction
        diagnostics = dict(coordinate=float(z))
        for label, weight in self.weight.items():
            diagnostics[label+'_residual_rms_m'] = float(np.sqrt(np.sum(weight[:, None]*residual**2)))
        if not np.isfinite(z) or not np.isfinite(list(diagnostics.values())).all():
            raise InvalidReuse('nonfinite geometry projection')
        return z, diagnostics

    def predict(self, arrays, time):
        age = float(time-self.time)
        if abs(age) > LIMITS['maximum_age_seconds']+1e-14:
            raise InvalidReuse('operator age limit')
        z, diagnostics = self.coordinates(arrays)
        if abs(z) > LIMITS['maximum_normalized_coordinate']:
            raise InvalidReuse('operator extrapolation limit')
        for label in self.weight:
            if diagnostics[label+'_residual_rms_m'] > LIMITS['maximum_'+label+'_residual_rms_m']:
                raise InvalidReuse(label+' geometry leaves the local direction')
        diagnostics.update(age_seconds=age, anchor_rhs_index=self.rhs_index)
        return self.energy+z*self.coefficients[0]+.5*z*z*self.coefficients[1], diagnostics

    def check_direction(self, plus, minus):
        self._same_basis(plus)
        self._same_basis(minus)
        change = (plus['points']-minus['points'])/2
        z = self._inner(self.direction, change)/self.norm2
        residual = change-z*self.direction
        diagnostics = {}
        for label, weight in self.weight.items():
            norm2 = np.sum(weight[:, None]*change**2)
            error2 = np.sum(weight[:, None]*residual**2)
            fraction = float(np.sqrt(error2/max(norm2, 1e-30)))
            diagnostics[label+'_direction_residual_fraction'] = fraction
            if not np.isfinite(fraction) or fraction > LIMITS['maximum_direction_residual_fraction']:
                raise InvalidReuse(label+' pressure direction leaves the local direction')
        return diagnostics

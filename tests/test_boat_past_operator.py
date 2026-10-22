"""Two-direction similarity on the actual closed boat fluid boundary."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

from topwater.rankine import TrianglePanels
from .experiments.free_surface.boat_world_boundary import construct
from .experiments.free_surface.operator import assemble_energy
from .experiments.free_surface.boat_past_operator import fit, predict, PastLocalOperator
from .experiments.free_surface.boat_operator_reuse import InvalidReuse, area_weights


def test_past_fit_separates_dilation_from_rigid_translation(tmp_path):
    _, base, _ = construct(tmp_path/'boat.json', export_template=True)
    shift = np.array([.0005, -.0003, .0002])
    epsilon = 1e-4

    def geometry(scale, translation):
        return dict(base, points=scale*base['points']+translation)

    def energy(arrays):
        return assemble_energy(arrays['faces'],
            TrianglePanels(arrays['points'][arrays['faces']], 1), order=6)[0]

    previous = [geometry(.999+delta, shift) for delta in (0., epsilon, -epsilon)]
    current = [geometry(1+delta, np.zeros(3)) for delta in (0., epsilon, -epsilon)]
    matrices = [energy(arrays) for arrays in [*previous, *current]]
    data, diagnostics = fit([*previous, *current], matrices)
    target = geometry(1.0004, .6*shift)
    actual = energy(target)
    coordinate_roundoff = 64*np.finfo(float).eps*np.max(abs(base['points']))
    # Translation changes neither the Laplace operator nor fixed-material
    # energy. Dilation multiplies it by length. These independent physical
    # identities distinguish two geometry directions without an invented matrix.
    assert_allclose(actual, 1.0004*matrices[3], rtol=0, atol=2e-11)
    for ordinary in (False, True):
        predicted, projection = predict(data, target, ordinary=ordinary)
        assert_allclose(predicted, actual, rtol=0, atol=2e-9)
        assert projection['body_residual_rms_m'] < coordinate_roundoff
        assert projection['free_residual_rms_m'] < coordinate_roundoff
    assert diagnostics['constrained_current_odd_chord_relative_error'] < 2e-14
    curved, _ = predict(data, target, curved=True)
    assert_allclose(curved, actual, rtol=0, atol=2e-9)
    for arrays, exact in zip(current, matrices[3:]):
        assert_allclose(predict(data, arrays, curved=True)[0], exact, rtol=0, atol=2e-13)
    runtime = PastLocalOperator((previous, matrices[:3], .05, 1), (current, matrices[3:], .1, 2), .1, 2)
    predicted, _ = runtime.predict(target, .13)
    assert_allclose(predicted, actual, rtol=0, atol=2e-9)
    runtime.check_direction(current[1], current[2])
    assert 'past_ordinary_coefficients' not in runtime.data
    curved_runtime = PastLocalOperator((previous, matrices[:3], .05, 1),
        (current, matrices[3:], .1, 2), .1, 2, curved=True)
    assert_allclose(curved_runtime.predict(target, .13)[0], actual, rtol=0, atol=2e-9)
    assert 'past_constrained_coefficients' not in curved_runtime.data
    ownership = base['body_source_group'].copy()
    index = np.flatnonzero(ownership != ownership[0])[0]
    ownership[0], ownership[index] = ownership[index], ownership[0]
    with pytest.raises(InvalidReuse, match='ownership'):
        runtime.predict(dict(base, body_source_group=ownership), .1)
    with pytest.raises(ValueError, match='independent second direction'):
        fit([*current, *current], [*matrices[3:], *matrices[3:]])

    # A parabolic similarity path has two independent geometry directions but
    # cannot distinguish transverse response from longitudinal curvature. Its
    # exact energy is still known from dilation/translation, so no invented
    # matrix is needed to test this identifiability failure.
    weights = area_weights(base)
    metric = .5*(weights['body']+weights['free'])
    center = metric@base['points']
    transverse = np.cross(center, [1., 0., 0.])
    transverse /= np.linalg.norm(transverse)
    parameters = (-4., -3.8, -4.2, 0., 1., -1.)
    parabola = [geometry(1+epsilon*t, epsilon*t*t*transverse) for t in parameters]
    exact = [(1+epsilon*t)*matrices[3] for t in parameters]
    with pytest.raises(ValueError, match='unobservable after current constraints'):
        fit(parabola, exact, ordinary=False)

    # Current rigid translation has zero operator increment in exact arithmetic.
    # Its independently assembled roundoff must not be judged by a purely
    # relative increment tolerance. The previous dilation still identifies the
    # transverse response.
    translated = [geometry(1., delta*transverse) for delta in (0., epsilon, -epsilon)]
    translated_energy = [matrices[3], energy(translated[1]), energy(translated[2])]
    invariant, _ = fit([*previous, *translated], [*matrices[:3], *translated_energy], ordinary=False)
    for arrays, exact in zip(translated, translated_energy):
        assert_allclose(predict(invariant, arrays, curved=True)[0], exact, rtol=0, atol=2e-13)

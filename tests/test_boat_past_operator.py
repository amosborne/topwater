"""Two-direction similarity on the actual closed boat fluid boundary."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

from topwater.rankine import TrianglePanels
from .experiments.free_surface.boat_world_boundary import construct
from .experiments.free_surface.operator import assemble_energy
from .experiments.free_surface.boat_past_operator import fit, predict, PastLocalOperator
from .experiments.free_surface.boat_operator_reuse import InvalidReuse


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
    runtime = PastLocalOperator((previous, matrices[:3], .05, 1), (current, matrices[3:], .1, 2), .1, 2)
    predicted, _ = runtime.predict(target, .13)
    assert_allclose(predicted, actual, rtol=0, atol=2e-9)
    runtime.check_direction(current[1], current[2])
    assert 'past_ordinary_coefficients' not in runtime.data
    ownership = base['body_source_group'].copy()
    index = np.flatnonzero(ownership != ownership[0])[0]
    ownership[0], ownership[index] = ownership[index], ownership[0]
    with pytest.raises(InvalidReuse, match='ownership'):
        runtime.predict(dict(base, body_source_group=ownership), .1)
    with pytest.raises(ValueError, match='independent second direction'):
        fit([*current, *current], [*matrices[3:], *matrices[3:]])

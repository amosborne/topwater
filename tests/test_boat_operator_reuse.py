"""Physical similarity and ownership checks on the actual closed boat domain."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

from topwater.rankine import TrianglePanels
from .experiments.free_surface.boat_operator_reuse import LocalOperator, InvalidReuse
from .experiments.free_surface.boat_world_boundary import construct
from .experiments.free_surface.operator import assemble_energy


def test_local_model_preserves_laplace_similarity_and_rejects_unlearned_motion(tmp_path):
    _, base, _ = construct(tmp_path/'boat.json', export_template=True)

    def scaled(scale):
        return dict(base, points=scale*base['points'])

    def energy(arrays):
        return assemble_energy(arrays['faces'],
            TrianglePanels(arrays['points'][arrays['faces']], 1), order=6)[0]

    # In 3D, fixed material Dirichlet data has energy proportional to length.
    # Train with independent actual operators; the similarity law supplies the
    # expected held-out energy and its slope, not the predictor's own formula.
    epsilon = 1e-4
    train = [base, scaled(1+epsilon), scaled(1-epsilon)]
    matrices = [energy(arrays) for arrays in train]
    model = LocalOperator(train, matrices, .1, 1)
    scale = 1+4*epsilon
    predicted, diagnostics = model.predict(scaled(scale), .104)
    actual = energy(scaled(scale))
    assert_allclose(actual, scale*matrices[0], rtol=0, atol=2e-11)
    assert_allclose(predicted, actual, rtol=0, atol=2e-9)
    potential = base['points']@np.array([.7, -.2, .4])
    assert_allclose(potential@predicted@potential,
                    scale*(potential@matrices[0]@potential), rtol=2e-11)
    assert diagnostics['body_residual_rms_m'] < 2e-15
    assert diagnostics['free_residual_rms_m'] < 2e-15
    model.check_direction(train[1], train[2])

    # Identical geometry with different physical panel ownership is not the
    # same operator basis. Nor can dilation training certify rigid translation.
    ownership = base['body_source_group'].copy()
    different = np.flatnonzero(ownership != ownership[0])[0]
    ownership[0], ownership[different] = ownership[different], ownership[0]
    with pytest.raises(InvalidReuse, match='ownership'):
        model.predict(dict(base, body_source_group=ownership), .1)
    with pytest.raises(InvalidReuse, match='pressure direction'):
        model.check_direction(dict(base, points=base['points']+[.001, 0, 0]),
                              dict(base, points=base['points']-[.001, 0, 0]))
    with pytest.raises(InvalidReuse, match='age'):
        model.predict(base, .12)

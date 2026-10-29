"""Native coefficient exchange only; no Capytaine/MATLAB model or force solve."""
import numpy as np
import pytest
from scipy.io import loadmat, savemat

from topwater.solve.wecsim.exchange import _native_excitation_payload


def coefficient_fixture():
    omega = np.array([.1, .2])
    headings = np.array([0., np.pi/4, np.pi/2, 2*np.pi])
    marker = np.arange(48, dtype=float).reshape(2, 4, 6)
    fk = 2+marker+1j*(1+.5*marker)
    sc = -1+.2*marker-1j*(.3+.1*marker)
    fk[:, -1] = fk[:, 0]
    sc[:, -1] = sc[:, 0]
    total = fk+sc
    stiffness = np.diag(np.arange(1., 7.))
    stiffness[2, 3] = stiffness[3, 2] = .1
    stiffness[2, 4] = stiffness[4, 2] = -.3
    return omega, headings, fk, sc, total, stiffness


def test_native_excitation_payload_preserves_axes_phase_and_static_limit():
    omega, headings, fk, sc, total, stiffness = coefficient_fixture()
    payload = _native_excitation_payload(omega, headings, fk, sc, total, stiffness)
    assert set(payload) == {
        'excitation_omega_rad_s', 'excitation_headings_rad',
        'excitation_froude_krylov_native', 'excitation_diffraction_native',
        'excitation_native', 'excitation_zero_limit_native',
    }
    np.testing.assert_array_equal(payload['excitation_omega_rad_s'], omega)
    np.testing.assert_array_equal(payload['excitation_headings_rad'], headings)
    np.testing.assert_array_equal(payload['excitation_froude_krylov_native'], fk)
    np.testing.assert_array_equal(payload['excitation_diffraction_native'], sc)
    np.testing.assert_array_equal(payload['excitation_native'], total)
    np.testing.assert_array_equal(payload['excitation_zero_limit_native'], stiffness[:, 2])
    assert payload['excitation_native'].shape == (2, 4, 6)
    assert np.any(payload['excitation_native'].imag != 0), 'Native phase must not be discarded'


def test_mat_exchange_keeps_omega_heading_mode_axes_and_native_complex_values(tmp_path):
    values = coefficient_fixture()
    payload = _native_excitation_payload(*values)
    path = tmp_path/'native.mat'
    savemat(path, payload)
    restored = loadmat(path)
    for name in ('excitation_froude_krylov_native', 'excitation_diffraction_native', 'excitation_native'):
        np.testing.assert_array_equal(restored[name], payload[name])
    for name in ('excitation_omega_rad_s', 'excitation_headings_rad', 'excitation_zero_limit_native'):
        np.testing.assert_array_equal(restored[name].ravel(), payload[name])


@pytest.mark.parametrize('index,name', [
    (2, 'excitation_froude_krylov_native'),
    (3, 'excitation_diffraction_native'),
    (4, 'excitation_native'),
])
def test_nonfinite_dense_coefficient_cannot_cross_native_boundary(index, name):
    values = list(coefficient_fixture())
    values[index][1, 1, 2] = complex(np.nan, 0.)
    with pytest.raises(ValueError, match=name+'.*finite'):
        _native_excitation_payload(*values)


def test_periodic_seam_and_genuine_component_sum_are_required():
    values = list(coefficient_fixture())
    values[2][:, -1, 0] += 1.
    with pytest.raises(ValueError, match='FK excitation periodic seam'):
        _native_excitation_payload(*values)
    values = list(coefficient_fixture())
    values[4][1, 1, 2] += 1.
    with pytest.raises(ValueError, match='FK plus diffraction'):
        _native_excitation_payload(*values)


def test_exchange_rejects_transposed_axes_and_changed_low_frequency_policy():
    values = list(coefficient_fixture())
    values[2] = values[2].transpose(1, 0, 2)
    with pytest.raises(ValueError, match='excitation_froude_krylov_native.*shape'):
        _native_excitation_payload(*values)
    values = list(coefficient_fixture())
    values[0] = np.array([.2, .3])
    with pytest.raises(ValueError, match='0.1 rad/s first coefficient'):
        _native_excitation_payload(*values)
    values = list(coefficient_fixture())
    values[1][-1] = np.pi
    with pytest.raises(ValueError, match='0/2pi periodic seam'):
        _native_excitation_payload(*values)

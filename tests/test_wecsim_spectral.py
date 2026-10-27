from types import SimpleNamespace

import numpy as np
import pytest

from topwater._wecsim_spectral import prepare_spectral_inputs


def case():
    omega = np.array([0., .025, .1, .15, .35])
    C = np.zeros((2, 2, len(omega)), complex)
    C[0, 0, 1:] = [-.2, -.1+.2j, .3-.4j, .08+.02j]
    C[1, 0, 1:] = [.1, .03+.04j, -.1+.07j, .02-.01j]
    C[0, 1, 1:] = [.03-.02j, .04+.01j, -.05+.02j, .01+.01j]
    C[1, 1, 1:] = [.02+.01j, -.01+.03j, .01-.01j, .02-.03j]
    train = SimpleNamespace(omega=omega, theta=np.array([np.pi/4, np.pi/4]), C=C,
                            carriers=SimpleNamespace(amplitudes=np.array([.2, .3])))
    source_w = np.array([.1, .2, .4, 1.])
    theta = np.linspace(0, 2*np.pi, 5)
    H = np.empty((len(source_w), len(theta), 6), complex)
    for mode in range(6):
        H[:, :, mode] = ((mode+1)*(source_w[:, None]**2+.1*np.cos(theta)[None, :])
                        + 1j*(.5*mode+source_w[:, None]*(1+np.sin(theta)[None, :])))
    H[:, -1] = H[:, 0]
    K = np.diag([2., 3., 40., 5., 6., 7.])
    K[:, 2] = [2., -3., 40., .5, -.8, 1.2]
    K[2] = K[:, 2]
    hydro = SimpleNamespace(excitation_omega=source_w, excitation_headings=theta,
                            excitation=H, K=K, density=1025., gravity=9.81)
    return SimpleNamespace(train=train), hydro


def dimensional_coefficients(inputs, hydro, frequency_index, heading_index):
    table = inputs.hydro['hydro_coeffs']['excitation']
    return ((table['re'][:, heading_index, frequency_index]
             - 1j*table['im'][:, heading_index, frequency_index])*hydro.density*hydro.gravity)


def test_coherent_branches_and_native_spectrum_preserve_physical_components():
    setup, hydro = case()
    original = setup.train.C.copy()
    raw_H = hydro.excitation.copy()
    data = prepare_spectral_inputs(setup, hydro)
    # Sorted physical directions are opposite branch -135deg, then +45deg.
    expected = np.c_[original[:, 1, 1:].sum(axis=0), original[:, 0, 1:].sum(axis=0)]
    np.testing.assert_array_equal(data.coefficients, expected)
    np.testing.assert_array_equal(data.wave['directions'].reshape(-1), [-135., 45.])
    native_amplitude = np.sqrt(2*(data.wave['spread']*data.wave['spectrum']/(2*np.pi))
                               *data.domega_rad_s[:, None]*data.dtheta_rad[None, :])
    rebuilt = native_amplitude*np.exp(-1j*data.wave['phase'])
    np.testing.assert_allclose(rebuilt, expected, rtol=0, atol=2e-15)
    times = np.array([0., .137, 2.])
    physical = (np.exp(-1j*times[:, None]*setup.train.omega[1:])
                @original[:, :, 1:].sum(axis=(0, 1))).real
    prepared = (np.exp(-1j*times[:, None]*data.native_omega_rad_s)@rebuilt.sum(axis=1)).real
    np.testing.assert_allclose(prepared, physical, rtol=0, atol=2e-15)
    np.testing.assert_array_equal(setup.train.C, original)
    np.testing.assert_array_equal(hydro.excitation, raw_H)
    assert set(data.wave) == {'frequencies', 'directions', 'spectrum', 'spread', 'phase'}
    assert set(data.hydro) == {'simulation_parameters', 'hydro_coeffs'}
    assert set(data.hydro['simulation_parameters']) == {'w', 'direction'}
    assert set(data.hydro['hydro_coeffs']) == {'excitation', 'mean_drift'}
    assert set(data.hydro['hydro_coeffs']['excitation']) == {'re', 'im'}


def test_standard_coefficients_match_independent_bilinear_cell_and_quadratic_bridge():
    setup, hydro = case()
    data = prepare_spectral_inputs(setup, hydro)
    w = data.hydro['simulation_parameters']['w'].reshape(-1)
    theta = data.hydro['simulation_parameters']['direction'].reshape(-1)
    column = int(np.flatnonzero(theta == 45.)[0])
    # Heading pi/4 is the midpoint of two genuine source columns.
    anchor = (hydro.excitation[0, 0]+hydro.excitation[0, 1])/2
    upper = (hydro.excitation[1, 0]+hydro.excitation[1, 1])/2
    query = data.native_omega_rad_s[2]
    row = int(np.flatnonzero(w == query)[0])
    fraction = (query-.1)/(.2-.1)
    expected = anchor+fraction*(upper-anchor)
    np.testing.assert_allclose(dimensional_coefficients(data, hydro, row, column), expected, rtol=0, atol=2e-14)
    query = data.native_omega_rad_s[0]
    row = int(np.flatnonzero(w == query)[0])
    expected = hydro.K[:, 2]+(query/.1)**2*(anchor-hydro.K[:, 2])
    np.testing.assert_allclose(dimensional_coefficients(data, hydro, row, column), expected, rtol=0, atol=2e-14)
    np.testing.assert_array_equal(data.hydro['hydro_coeffs']['excitation']['re'][:, :, 0],
        np.repeat((hydro.K[:, 2]/(hydro.density*hydro.gravity))[:, None], len(theta), axis=1))
    assert np.all(data.hydro['hydro_coeffs']['excitation']['im'][:, :, 0] == 0)
    assert np.all(data.hydro['hydro_coeffs']['mean_drift'] == 0)
    np.testing.assert_array_equal(data.hydro['hydro_coeffs']['excitation']['re'][:, 0], data.hydro['hydro_coeffs']['excitation']['re'][:, -1])
    np.testing.assert_array_equal(data.hydro['hydro_coeffs']['excitation']['im'][:, 0], data.hydro['hydro_coeffs']['excitation']['im'][:, -1])


def test_hz_native_order_is_the_standard_query_identity_without_new_wave_field():
    setup, hydro = case()
    data = prepare_spectral_inputs(setup, hydro)
    native = (data.wave['frequencies'].reshape(-1)*2)*np.pi
    np.testing.assert_array_equal(data.native_omega_rad_s, native)
    np.testing.assert_array_equal(data.omega_roundtrip_error_rad_s, native-setup.train.omega[1:])
    assert np.isin(native, data.hydro['simulation_parameters']['w']).all()
    assert data.hydro['simulation_parameters']['w'][-1, 0] == hydro.excitation_omega[-1]
    assert len(data.wave['frequencies']) == len(setup.train.omega)-1
    assert 'omega_rad_s' not in data.wave


def test_nearby_distinct_physical_directions_are_not_quantized_or_merged():
    setup, hydro = case()
    setup.train.theta[:] = [.2, .2+1e-9]
    data = prepare_spectral_inputs(setup, hydro)
    assert data.coefficients.shape[1] == 4
    assert len(np.unique(data.wave['directions'])) == 4


def test_calm_and_zero_energy_rows_need_no_random_phases_or_denominator_floor():
    setup, hydro = case()
    setup.train.C[:] = 0
    setup.train.carriers.amplitudes[:] = 0
    data = prepare_spectral_inputs(setup, hydro)
    for field in ('spectrum', 'spread', 'phase'):
        assert np.all(data.wave[field] == 0)
    assert data.coefficients.shape == (4, 2)
    assert np.all(data.coefficients == 0)
    assert np.isfinite(data.hydro['hydro_coeffs']['excitation']['re']).all()


@pytest.mark.parametrize('failure', ['dc', 'nonfinite', 'outside', 'seam', 'inactive', 'underflow'])
def test_incomplete_or_unrepresentable_physical_handoffs_fail_clearly(failure):
    setup, hydro = case()
    if failure == 'dc':
        setup.train.C[0, 0, 0] = 1e-30
    elif failure == 'nonfinite':
        setup.train.C[0, 1, 2] = np.nan
    elif failure == 'outside':
        setup.train.omega[-1] = 2.
    elif failure == 'seam':
        hydro.excitation[0, -1, 0] += .01
    elif failure == 'inactive':
        setup.train.carriers.amplitudes[0] = 0
    else:
        setup.train.C[:, :, 1:] = 1e-200
    with pytest.raises(ValueError):
        prepare_spectral_inputs(setup, hydro)

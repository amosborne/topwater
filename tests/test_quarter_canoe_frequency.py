"""Exact scalar-input retry integrity; no PDE assembly or factorization."""
import numpy as np
import pytest
from scipy.sparse import csc_matrix

from tests.experiments.free_surface.quarter_canoe_frequency import save_input, reuse_input


def artifact(tmp_path):
    matrix = csc_matrix(np.array([[2.+1j, .3-.2j], [.3-.2j, 4.+2j]]))
    load = np.arange(12, dtype=float).reshape(2, 6)/7
    rhs = load[:, 2].astype(complex)
    report = dict(source_sha256={'operator': 'source', 'quarter_canoe_frequency.py': 'writer'}, mesh_report_sha256='mesh-report',
        mesh_npz_sha256='mesh-values', capytaine_report_sha256='reference', omega=2.,
        density=1000., gravity=9.81, parity='EE', computed_modes=['heave'],
        pml_width=[7., 7.], pml_strength=8., pml_onset=[2., 1.],
        volume_quadrature_order=12, surface_quadrature_order=12)
    path = tmp_path/'fixed.input.npz'
    save_input(path, matrix, load, rhs, report)
    return path, report, matrix, load, rhs


def test_persisted_complex_csc_and_rhs_reproduce_exact_input(tmp_path):
    path, report, matrix, load, rhs = artifact(tmp_path)
    actual, recovered_rhs, _ = reuse_input(path, report, load)
    np.testing.assert_array_equal(actual.data, matrix.data)
    np.testing.assert_array_equal(actual.indices, matrix.indices)
    np.testing.assert_array_equal(actual.indptr, matrix.indptr)
    np.testing.assert_array_equal(actual@np.array([1.-2j, 3.+4j]), matrix@np.array([1.-2j, 3.+4j]))
    np.testing.assert_array_equal(recovered_rhs, rhs)


@pytest.mark.parametrize('key,replacement', [('omega', 4.), ('source_sha256', {'operator': 'changed', 'quarter_canoe_frequency.py': 'writer'}),
                                         ('mesh_npz_sha256', 'different'), ('capytaine_report_sha256', 'other-reference')])
def test_reuse_rejects_changed_equations_geometry_or_reference(tmp_path, key, replacement):
    path, report, _, load, _ = artifact(tmp_path)
    report[key] = replacement
    with pytest.raises(ValueError, match='source, mesh, reference or fixed physical case differs'):
        reuse_input(path, report, load)


def test_reuse_rejects_corrupted_file_or_changed_physical_load(tmp_path):
    path, report, _, load, _ = artifact(tmp_path)
    changed = load.copy(); changed[0, 2] += .01
    with pytest.raises(ValueError, match='body load differs'):
        reuse_input(path, report, changed)
    with path.open('ab') as stream:
        stream.write(b'corruption')
    with pytest.raises(ValueError, match='file hash differs'):
        reuse_input(path, report, load)


def test_writer_ordering_change_reuses_same_input_without_overwriting_provenance(tmp_path):
    path, report, matrix, load, _ = artifact(tmp_path)
    report['source_sha256']['quarter_canoe_frequency.py'] = 'new-writer'
    report['lu_ordering'] = 'MMD_AT_PLUS_A'
    actual, _, manifest = reuse_input(path, report, load)
    np.testing.assert_array_equal(actual.data, matrix.data)
    assert manifest['writer_source_sha256'] == 'writer'
    with pytest.raises(FileExistsError, match='persisted input already exists'):
        save_input(path, matrix, load, load[:, 2].astype(complex), report)

"""Saved field values alone cannot establish their physical node identities."""
from types import SimpleNamespace

import numpy as np
import pytest

from .experiments.free_surface.p4_tetra import make_mesh
from .experiments.free_surface.p4_canoe_stage import verify_master_identity


def test_constant_surface_data_cannot_hide_reordered_mesh_masters(tmp_path):
    mesh = make_mesh(np.array([[0., 0., 0.], [1., 0., 0.],
                               [0., 1., 0.], [0., 0., -1.]]),
                     np.array([[0, 1, 2, 3]]))
    free = np.flatnonzero(mesh.points[:, 2] == 0)
    owner = SimpleNamespace(mesh=mesh, free_nodes=free)
    path = tmp_path/'stage.npz'
    # A flat, quiescent field has identical entries under every permutation.
    # Its values would pass a physical-value comparison with the wrong IDs.
    fields = dict(eta=np.zeros(len(free)), psi=np.zeros(len(free)))
    identity = dict(reference_points=mesh.points, mesh_elements=mesh.elements,
                    free_nodes=free)
    np.savez(path, **fields, **identity)
    assert 'exactly match' in verify_master_identity(owner, path)
    reordered = free[::-1]
    np.testing.assert_array_equal(fields['eta'], fields['eta'][::-1])
    np.savez(path, **fields, **dict(identity, free_nodes=reordered))
    with pytest.raises(ValueError, match='mismatch: free_nodes'):
        verify_master_identity(owner, path)
    np.savez(path, **fields, free_nodes=free)
    with pytest.raises(ValueError, match='all-or-none'):
        verify_master_identity(owner, path)
    np.savez(path, **fields)
    assert 'Legacy' in verify_master_identity(owner, path)

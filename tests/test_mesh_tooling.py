"""Bounded optional-package probes, executed only under coordinator allocation.

No dependency is installed or adopted here. Missing optional tooling is an
explicit skip. libigl execution additionally requires distribution-review
approval via TOPWATER_LIBIGL_PROBE=1; API availability is not inferred from C++.
"""

import importlib.metadata
import os
import time
import tracemalloc

import numpy as np
import pytest


def _gmsh_facet_probe(rotation, scale):
    gmsh = pytest.importorskip('gmsh')
    if gmsh.isInitialized():
        pytest.fail('probe requires its own Gmsh instance; will not modify an existing instance')
    outer = np.array([[-2., -2., 0.], [0., -2., 0.], [2., -2., 0.],
                      [2., 2., 0.], [-2., 2., 0.]])
    hole = np.array([[-.5, -.5, 0.], [-.5, .5, 0.], [.5, .5, 0.], [.5, -.5, 0.]])
    original = np.vstack((outer, hole))
    points = (original @ rotation.T + [3., -1., .7])*scale
    started = time.perf_counter()
    tracemalloc.start()
    gmsh.initialize(['topwater-bounded-facet-probe', '-nopopup'])
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.option.setNumber('Mesh.MaxNumThreads1D', 1)
        gmsh.option.setNumber('Mesh.MaxNumThreads2D', 1)
        gmsh.option.setNumber('Mesh.ElementOrder', 1)
        gmsh.option.setNumber('Mesh.MeshSizeMin', .25*scale)
        gmsh.option.setNumber('Mesh.MeshSizeMax', .5*scale)
        gmsh.model.add('one-facet-with-hole')
        point_tags = [gmsh.model.geo.addPoint(*p) for p in points]
        loops, curves = [], []
        for indices in (np.arange(5), np.arange(5, 9)):
            lines = [gmsh.model.geo.addLine(point_tags[i], point_tags[j])
                     for i, j in zip(indices, np.roll(indices, -1))]
            curves.extend(zip(lines, indices, np.roll(indices, -1)))
            loops.append(gmsh.model.geo.addCurveLoop(lines))
        surface = gmsh.model.geo.addPlaneSurface(loops)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        coordinates = np.asarray(coordinates).reshape(-1, 3)
        lookup = {int(tag): p for tag, p in zip(node_tags, coordinates)}
        assert len(lookup) == len(node_tags)
        for tag, p in zip(point_tags, points):
            _, xyz, _ = gmsh.model.mesh.getNodes(0, tag)
            np.testing.assert_allclose(np.asarray(xyz).reshape(-1, 3), p[None], atol=scale*1e-10, rtol=0)
        for curve, i, j in curves:
            types, _, connectivity = gmsh.model.mesh.getElements(1, curve)
            assert len(types) == 1 and types[0] == 1, 'first-order constrained edges required'
            segments = np.asarray(connectivity[0]).reshape(-1, 2)
            segment_points = np.array([[lookup[int(tag)] for tag in edge] for edge in segments])
            length = np.linalg.norm(points[j]-points[i])
            np.testing.assert_allclose(np.linalg.norm(segment_points[:, 1]-segment_points[:, 0], axis=1).sum(),
                                       length, atol=scale*1e-9, rtol=0)
            direction = (points[j]-points[i])/length
            projection = (segment_points-points[i]) @ direction
            residual = segment_points-points[i]-projection[:, :, None]*direction
            assert np.max(np.linalg.norm(residual, axis=2)) <= scale*1e-9
            assert projection.min() >= -scale*1e-9 and projection.max() <= length+scale*1e-9
            _, degree = np.unique(segments, return_counts=True)
            assert np.count_nonzero(degree == 1) == 2 and np.all((degree == 1) | (degree == 2))
        types, _, connectivity = gmsh.model.mesh.getElements(2, surface)
        assert len(types) == 1 and types[0] == 2, 'triangles required'
        faces = np.asarray(connectivity[0]).reshape(-1, 3)
        triangles = np.array([[lookup[int(tag)] for tag in face] for face in faces])
        cross = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
        normal = rotation @ np.array([0., 0., 1.])
        assert np.all(cross @ normal > 0)
        np.testing.assert_allclose(np.linalg.norm(cross, axis=1).sum()/2, 15.*scale**2,
                                   rtol=1e-10, atol=scale**2*1e-10)
        local = ((triangles/scale-[3., -1., .7]) @ rotation)
        centroids = local.mean(axis=1)
        assert not np.any((abs(centroids[:, 0]) < .5-1e-10) & (abs(centroids[:, 1]) < .5-1e-10))
        assert len(np.unique(np.sort(faces, axis=1), axis=0)) == len(faces)
        assert set(np.unique(faces)) == set(int(x) for x in node_tags)
        return dict(version=gmsh.__version__, triangles=len(faces), vertices=len(node_tags))
    finally:
        gmsh.finalize()
        _, python_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        print(dict(probe='gmsh-constrained-facet', wall_seconds=time.perf_counter()-started,
                   python_peak_bytes=python_peak,
                   memory_limit='Python allocations only; native Gmsh peak memory not measured'))


def test_gmsh_constrained_facet_preserves_hole_and_collinear_boundary():
    from scipy.spatial.transform import Rotation
    first = _gmsh_facet_probe(np.eye(3), 1.)
    second = _gmsh_facet_probe(Rotation.from_rotvec([.4, -.2, .1]).as_matrix(), 1e-3)
    print(dict(original=first, transformed=second))


def test_libigl_distribution_and_required_binding_availability():
    try:
        distribution = importlib.metadata.distribution('libigl')
    except importlib.metadata.PackageNotFoundError:
        pytest.skip('libigl Python distribution not installed; no installation performed')
    if os.environ.get('TOPWATER_LIBIGL_PROBE') != '1':
        pytest.skip(f'libigl {distribution.version} installed; distribution-review approval required before execution')
    started = time.perf_counter()
    import igl
    required = ('boundary_loop', 'is_edge_manifold', 'is_vertex_manifold')
    availability = {name: callable(getattr(igl, name, None)) for name in required}
    print(dict(probe='libigl-python-bindings', distribution=distribution.metadata['Name'],
               version=distribution.version, available=availability,
               geometric_self_intersections=callable(getattr(igl, 'selfintersect', None))))
    assert all(availability.values()), 'required topology bindings unavailable in this Python distribution'
    # One triangle has one known boundary; topology capability does not imply
    # robust geometric self-intersection detection or a Boolean implementation.
    faces = np.array([[0, 1, 2]], dtype=np.int64)
    assert set(np.asarray(igl.boundary_loop(faces)).ravel()) == {0, 1, 2}
    print(dict(probe='libigl-boundary-loop', wall_seconds=time.perf_counter()-started))

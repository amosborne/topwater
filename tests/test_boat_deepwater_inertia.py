"""The new boat comparison preserves physical geometry and seawater units."""
import json
import numpy as np
from numpy.testing import assert_allclose
import pytest

from .reference_boat import equilibrium, reference_body, WATER_DENSITY, GRAVITY
from .experiments.capytaine_mesh_cost import subdivide
from .capytaine_reference import instantaneous_inertia


def test_three_panel_levels_preserve_wet_surface_and_rigid_flux_moments():
    pose, hydro = equilibrium()
    triangles = hydro.wetted_triangles
    previous = None
    for level, expected in enumerate((17, 68, 272)):
        assert len(triangles) == expected
        assert triangles[:, :, 2].max() <= 1e-14
        normal_area = np.cross(triangles[:, 1]-triangles[:, 0],
                               triangles[:, 2]-triangles[:, 0])/2
        area = np.linalg.norm(normal_area, axis=1)
        centroid = triangles.mean(axis=1)
        # Independent geometric integrals include the rotational Neumann load,
        # not merely matching total triangle count/area.
        moments = np.r_[area.sum(), (area[:, None]*centroid).sum(axis=0),
            normal_area.sum(axis=0), np.cross(centroid-pose.position, normal_area).sum(axis=0)]
        if previous is not None:
            assert_allclose(moments, previous, atol=3e-14, rtol=2e-14)
        previous = moments
        triangles = subdivide(triangles)


def test_capytaine_infinite_frequency_seawater_scaling_and_gravity_independence():
    pytest.importorskip('capytaine')
    pose, hydro = equilibrium()
    body = reference_body()
    default, damping = instantaneous_inertia(hydro.wetted_triangles, pose.position, body.mass)
    seawater, seawater_damping = instantaneous_inertia(hydro.wetted_triangles, pose.position,
        body.mass, density=WATER_DENSITY, gravity=9.7)
    # At omega=infinity gravity drops out of the Dirichlet-image problem;
    # changing rho scales ALL force and moment couplings, not just diagonals.
    assert_allclose(seawater, default*(WATER_DENSITY/1000.), rtol=2e-14, atol=2e-11)
    assert np.array_equal(damping, np.zeros((6, 6)))
    assert np.array_equal(seawater_damping, damping)
    for density, gravity in ((0., GRAVITY), (WATER_DENSITY, np.nan)):
        with pytest.raises(ValueError, match='positive finite'):
            instantaneous_inertia(hydro.wetted_triangles, pose.position, body.mass,
                                  density=density, gravity=gravity)


def test_component_failure_retains_geometry_and_incomplete_checkpoint(tmp_path, monkeypatch):
    pytest.importorskip('capytaine')
    from .experiments import boat_deepwater_inertia as experiment
    failure = ArithmeticError('synthetic assembly failure')

    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(experiment, 'instantaneous_added_mass', fail)
    output = tmp_path/'failed.json'
    with pytest.raises(ArithmeticError) as caught:
        experiment.run(output)
    assert caught.value is failure
    report = json.loads(output.read_text())
    assert report['complete'] is False and report['qualified'] is False
    assert report['levels'][0]['complete'] is False
    assert report['error']['message'] == str(failure)
    with np.load(output.with_suffix('.npz'), allow_pickle=False) as archive:
        assert archive['triangles_0'].shape == (17, 3, 3)
        assert 'topwater_0' not in archive

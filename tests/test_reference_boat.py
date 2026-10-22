"""Physical checks on the authored boat, without a fluid operator solve."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation

from topwater import HullMesh, InitialState
from topwater.hydrostatics import FloodingError, LossOfContactError, hydrostatics
from .reference_boat import (
    BEAM, CHINE_BEAM_FRACTION, DEADRISE_DEGREES, DEPTH, GRAVITY, LENGTH,
    MOTOR_MASS, MOTOR_POSITION, PAYLOAD_MASS, PAYLOAD_POSITION, SHELL_MASS,
    WATER_DENSITY, equilibrium, mesh_arrays, reference_body,
)


def test_full_envelope_volume_and_cap_have_physical_identity():
    vertices, faces, rim, cap = mesh_arrays()
    mesh = HullMesh(vertices, faces, rim, cap)
    triangles = vertices[faces]
    volume = np.einsum('fi,fi->f', triangles[:, 0],
                      np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    b = BEAM/2
    c = CHINE_BEAM_FRACTION*b
    chine_z = c*np.tan(np.deg2rad(DEADRISE_DEGREES))
    section_area = DEPTH*(b+c)-b*chine_z
    # A constant section over 3/4 L plus a geometrically similar tapered cone.
    assert_allclose(volume, section_area*LENGTH*(3/4+1/12), rtol=2e-15)
    assert_allclose(vertices[rim, 2], DEPTH, atol=0)
    assert_allclose(vertices[faces[cap], 2], DEPTH, atol=0)
    stern = np.all(vertices[faces, 0] == -LENGTH/2, axis=1)
    assert stern.sum() == 3
    assert not np.intersect1d(np.flatnonzero(stern), cap).size
    assert len(mesh.physical_faces) == 15


def test_mass_tensor_matches_independent_triangle_quadrature():
    body = reference_body()
    triangles = body.mesh.vertices[body.mesh.faces[body.mesh.physical_faces]]
    area = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                   triangles[:, 2]-triangles[:, 0]), axis=1)/2
    # Independent degree-two quadrature, not the shell raw-moment formula.
    barycentric = np.full((3, 3), 1/6)
    np.fill_diagonal(barycentric, 2/3)
    samples = np.einsum('qa,fai->fqi', barycentric, triangles).reshape(-1, 3)
    dm = np.repeat(SHELL_MASS*area/area.sum()/3, 3)
    samples = np.vstack((samples, PAYLOAD_POSITION, MOTOR_POSITION))
    dm = np.r_[dm, PAYLOAD_MASS, MOTOR_MASS]
    cg = np.sum(dm[:, None]*samples, axis=0)/dm.sum()
    r = samples-cg
    inertia = np.sum(dm[:, None, None]*(np.sum(r*r, axis=1)[:, None, None]*np.eye(3)
                                        -np.einsum('ni,nj->nij', r, r)), axis=0)
    assert body.mass == 1250
    assert_allclose(body.cg, cg, atol=3e-16)
    assert_allclose(body.inertia_cg, inertia, atol=3e-12, rtol=2e-15)
    assert abs(body.inertia_cg[0, 2]) > 1  # Preserve longitudinal/vertical coupling.
    assert_allclose(body.inertia_cg[[0, 1], [1, 2]], 0, atol=3e-14)


def test_planar_subdivision_preserves_mass_geometry_and_rest():
    original = reference_body()
    rest, result = equilibrium()
    for refinement in (1, 2):
        body = reference_body(refinement)
        assert len(body.mesh.faces) == 18*4**refinement
        assert_allclose(body.cg, original.cg, atol=4e-16)
        assert_allclose(body.inertia_cg, original.inertia_cg, atol=3e-12, rtol=2e-15)
        refined_rest, refined_result = equilibrium(refinement)
        assert_allclose(refined_rest.position, rest.position, atol=2e-12)
        assert_allclose(refined_rest.orientation, rest.orientation, atol=2e-12)
        assert_allclose(refined_result.displaced_volume, result.displaced_volume, atol=2e-12)
        assert_allclose(refined_result.stiffness, result.stiffness, atol=1e-7, rtol=1e-11)


def test_actual_mesh_floats_with_pitch_balance_freeboard_and_positive_stability():
    body = reference_body()
    rest, result = equilibrium()
    assert_allclose(result.displaced_volume, body.mass/WATER_DENSITY, atol=2e-12)
    assert_allclose(result.force, [0, 0, body.mass*GRAVITY], atol=3e-8)
    assert_allclose(result.moment, 0, atol=3e-8)
    assert result.rim_clearance > .5
    # Heave/roll/pitch stability at equilibrium, retaining cross couplings.
    block = result.stiffness[np.ix_([2, 3, 4], [2, 3, 4])]
    assert np.linalg.eigvalsh(block).min() > 0
    assert not np.isin(result.source_faces, body.mesh.artificial_faces).any()
    # Independent pressure-moment identity about the actual world CG.
    assert_allclose(result.moment, np.cross(result.center_of_buoyancy-rest.position,
                                          [0, 0, body.mass*GRAVITY]), atol=3e-8)
    angle = .001
    rotation = Rotation.from_quat(np.roll(rest.orientation, -1))
    for axis in (0, 1):
        tilted = InitialState(rest.position, np.roll(
            (Rotation.from_rotvec(np.eye(3)[axis]*angle)*rotation).as_quat(), 1))
        assert hydrostatics(body, tilted, density=WATER_DENSITY, gravity=GRAVITY).moment[axis] < 0


def test_open_boat_detects_rim_immersion_and_loss_of_water_contact():
    body = reference_body()
    rest, result = equilibrium()
    submerged = InitialState(rest.position-[0, 0, result.rim_clearance+.01], rest.orientation)
    with pytest.raises(FloodingError):
        hydrostatics(body, submerged, density=WATER_DENSITY, gravity=GRAVITY)
    airborne = InitialState(rest.position+[0, 0, 2.], rest.orientation)
    with pytest.raises(LossOfContactError):
        hydrostatics(body, airborne, density=WATER_DENSITY, gravity=GRAVITY)

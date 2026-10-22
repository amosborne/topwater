"""Roll reporting follows the same active body-to-world quaternion convention."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation

from .experiments.free_surface.boat_motion_metrics import attitude_degrees


@pytest.mark.parametrize('roll', [-30., -20., -.1, 0., .1, 20., 30.])
def test_roll_is_independent_of_heading_pitch_quaternion_sign_and_norm(roll):
    rotation = (Rotation.from_rotvec([0., 0., 1.7])
        *Rotation.from_rotvec([0., -.2, 0.])
        *Rotation.from_rotvec([np.deg2rad(roll), 0., 0.]))
    q = rotation.as_quat(scalar_first=True)
    assert_allclose(attitude_degrees(q), [roll, np.rad2deg(-.2), np.rad2deg(1.7)], atol=2e-14)
    assert_allclose(attitude_degrees(-2*q), attitude_degrees(q), atol=2e-14)


def test_heading_wrap_does_not_change_small_attitude_error():
    a = attitude_degrees(Rotation.from_rotvec([0., 0., np.deg2rad(179.98)]).as_quat(scalar_first=True))
    b = attitude_degrees(Rotation.from_rotvec([0., 0., np.deg2rad(-179.98)]).as_quat(scalar_first=True))
    assert_allclose((b-a+180)%360-180, [0., 0., .04], atol=2e-14)

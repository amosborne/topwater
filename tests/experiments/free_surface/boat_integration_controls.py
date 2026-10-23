"""Shared exploratory local RK controls; these are NOT motion-accuracy targets.

Body order: world position[m], raw wxyz quaternion[1], world velocity[m/s],
world angular velocity[rad/s]. Fluid blocks retain their own physical units.
SciPy RK23 uses an RMS error norm over the entire supplied state vector.
"""
import numpy as np

RELATIVE_TOLERANCE = .001
BODY_ABSOLUTE_TOLERANCE = np.r_[np.full(3, .001), np.full(4, .0001),
                              np.full(3, .001), np.full(3, .0001)]
BODY_ABSOLUTE_TOLERANCE.setflags(write=False)
SURFACE_HEIGHT_ABSOLUTE_TOLERANCE_M = .001
SURFACE_POTENTIAL_ABSOLUTE_TOLERANCE_M2_PER_S = .01
ERROR_NORM = 'SciPy RK23 RMS over the complete supplied state; raw quaternion and world angular velocity'

"""Analytic physical fields; independent of panel integration and solves."""

import numpy as np


def harmonic_field(points):
    x, y, _ = np.asarray(points).T
    potential = np.column_stack((x/3.5, (x*x-y*y)/3.5**2))
    gradient = np.zeros((len(points), 2, 3))
    gradient[:, 0, 0] = 1/3.5
    gradient[:, 1, 0], gradient[:, 1, 1] = 2*x/3.5**2, -2*y/3.5**2
    return potential, gradient


def point_source(points, source):
    delta = np.asarray(points)-source
    distance = np.linalg.norm(delta, axis=1)
    return (1/distance)[:, None], (-delta/distance[:, None]**3)[:, None, :]

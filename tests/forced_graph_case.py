"""Declared fixed/forced boundary data for the moving-body coupling experiment.

Only the initial observed-sea tangent sets geometry. This is not a sea-state
realization evolved by a free-surface dynamic boundary condition.
"""
import numpy as np
from topwater.constants import GRAVITY
from tests.experiments.free_surface.boat_hos_continuation import tilted_geometry, tangent_graph


def boundary_data(inputs, grid=129, length=64.):
    cg = inputs.initial_state.position
    _, level, slope = tilted_geometry(inputs)
    normal = np.r_[-slope, 1.]
    normal /= np.linalg.norm(normal)
    tangent = np.array([1., 0., slope[0]])
    tangent /= np.linalg.norm(tangent)
    frame = np.vstack((tangent, np.cross(normal, tangent), normal))
    origin = cg[:2]-length/2
    x, y = np.meshgrid(origin[0]+np.arange(grid)*length/grid,
                       origin[1]+np.arange(grid)*length/grid)
    eta, inner, outer = tangent_graph(inputs, x, y, level, slope)
    if np.any(outer >= length/2):
        raise ValueError('forced graph support must fit strictly within its periodic box')
    wavevector = 2*np.pi*np.array([3., 1.])/length
    return dict(eta=eta, origin=origin, length=length, grid_xy=np.column_stack((x.ravel(), y.ravel())),
        plane_origin=np.r_[cg[:2], level], frame=frame, inner=inner,
        amplitude=.1, wavevector=wavevector, frequency=float(np.sqrt(GRAVITY*np.linalg.norm(wavevector))))

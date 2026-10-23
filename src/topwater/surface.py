"""Exact planar-triangle load and volume integrals, independent of hull fixtures."""
import numpy as np


def body_load(mesh, panels, body, center):
    """Exact P1 weak load for six rigid motions, including waterline nodes."""
    points = mesh.vertices[mesh.faces[body]]
    normal = panels.normals[body, None]
    modes = np.concatenate((np.broadcast_to(normal, points.shape),
                            np.cross(points-center, normal)), axis=2)
    load = np.zeros((len(mesh.vertices), 6))
    total = modes.sum(axis=1)
    for corner in range(3):
        values = panels.areas[body, None]*(total+modes[:, corner])/12
        np.add.at(load, mesh.faces[body, corner], values)
    return load


def surface_integrals(vertices, faces):
    """Independent divergence-theorem reference integrals of planar triangles.

    Uses z*n_z for volume (the input validator uses tetrahedral determinants).
    First moments use exact triangle polynomial averages.
    """
    triangles = vertices[faces]
    p, q, r = triangles.transpose(1, 0, 2)
    normal_area = np.cross(q - p, r - p) / 2
    volume = np.sum(normal_area[:, 2] * triangles[:, :, 2].mean(axis=1))
    average_squared = (p*p + q*q + r*r + p*q + p*r + q*r) / 6
    first_moments = np.sum(normal_area * average_squared, axis=0) / 2
    return volume, first_moments / volume

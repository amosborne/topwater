"""Closed rectangular fluid domain for independent wave verification."""
import numpy as np
from topwater import HullMesh


def tank(nx, ny, nz, depth=4.):
    x = np.linspace(-2*np.pi, 2*np.pi, nx+1)
    y = np.linspace(-np.pi, np.pi, ny+1)
    z = np.linspace(-depth, 0, nz+1)
    vertices, faces, free = [], [], []
    for first, second, axis, level, reverse, surface in (
        (x, y, 2, 0, False, True), (x, y, 2, -depth, True, False),
        (x, z, 1, y[0], False, False), (x, z, 1, y[-1], True, False),
        (y, z, 0, x[0], True, False), (y, z, 0, x[-1], False, False),
    ):
        u, v = np.meshgrid(first, second, indexing='ij')
        points = np.empty((u.size, 3))
        axes = [j for j in range(3) if j != axis]
        points[:, axes[0]], points[:, axes[1]] = u.ravel(), v.ravel()
        points[:, axis] = level
        begin = sum(len(a) for a in vertices)
        vertices.append(points)
        for i in range(len(first)-1):
            for j in range(len(second)-1):
                a = begin+i*len(second)+j
                b, c, d = a+len(second), a+len(second)+1, a+1
                children = [[a, b, c], [a, c, d]]
                faces.extend([t[::-1] for t in children] if reverse else children)
                free.extend([surface, surface])
    points = np.vstack(vertices)
    unique, inverse = np.unique(points, axis=0, return_inverse=True)
    mesh = HullMesh(unique, inverse[np.array(faces)])
    nodes = inverse[:len(x)*len(y)].reshape(len(x), len(y))
    return mesh, np.array(free), nodes


def trace(points):
    x, y, _ = np.asarray(points).T
    return np.column_stack((np.ones(len(x)), np.cos(x+2*np.pi),
                            np.cos(.5*(x+2*np.pi))*np.cos(y+np.pi)))

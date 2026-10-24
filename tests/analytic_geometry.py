"""Mathematical meshes for kernel tests, without boat loading or sea inputs."""
import numpy as np


def sphere_mesh(refinement=2):
    """Outward unit octahedron, radially refined toward the unit sphere."""
    if isinstance(refinement, bool) or not isinstance(refinement, int) or refinement < 0:
        raise ValueError('nonnegative integer refinement required')
    vertices = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0],
                         [0, -1., 0], [0, 0, 1.], [0, 0, -1.]])
    faces = np.array([[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
                      [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]])
    for _ in range(refinement):
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        unique, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        midpoints = vertices[unique].mean(axis=1)
        midpoints /= np.linalg.norm(midpoints, axis=1)[:, None]
        ab, bc, ca = (inverse+len(vertices)).reshape(3, -1)
        a, b, c = faces.T
        faces = np.vstack([np.column_stack(v) for v in
                           ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
        vertices = np.vstack((vertices, midpoints))
    return np.ascontiguousarray(vertices), np.ascontiguousarray(faces)

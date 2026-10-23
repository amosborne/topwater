"""Aperiodic Poisson lift of one continuous compact structured P1 trace.

Static CPU math reference only. The source plane is z=0; arrays use (y,x).
No evolved hidden water, FFT coupling, pressure solve or time integrator.
"""
import numpy as np
from topwater.rankine import TrianglePanels, single_layer_trace
from .panel_gradient_reference import checked_influence


class CompactP1Trace:
    """Conforming tensor patch, bottom-left/top-right diagonals, zero perimeter.

    For each triangle t, b_t=grad_xy(f), K_t=int dG/dn_source and
    g_t=grad_target int G with upward source normals and G=1/(4pi r):
      Hf = -2 sum[f_t(xy) K_t + z b_t dot g_t,xy],
      u_xy = -2 sum[K_t b_t],  u_z = -2 sum[b_t dot g_t,xy].
    The velocity identities integrate by parts: shared-edge terms cancel only
    for a continuous trace with ZERO perimeter. No new layer formula is used.
    """
    def __init__(self, nodesx, nodesy, values):
        raw = [np.asarray(a) for a in (nodesx, nodesy, values)]
        if any(np.iscomplexobj(a) or not np.isfinite(a).all() for a in raw):
            raise ValueError('finite real axes and trace required')
        x, y, values = [np.array(a, dtype=np.float64, copy=True) for a in raw]
        if (x.ndim != 1 or y.ndim != 1 or min(len(x), len(y)) < 2
                or np.any(np.diff(x) <= 0) or np.any(np.diff(y) <= 0)
                or values.shape != (len(y), len(x))):
            raise ValueError('increasing axes and matching (y,x) trace required')
        xx, yy = np.meshgrid(x, y)
        vertices = np.column_stack((xx.ravel(), yy.ravel(), np.zeros(xx.size)))
        row, col = np.indices((len(y)-1, len(x)-1))
        a = (row*len(x)+col).ravel()
        faces = np.concatenate((np.column_stack((a, a+1, a+len(x)+1)),
                                np.column_stack((a, a+len(x)+1, a+len(x)))))
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        unique, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
        boundary = unique[counts == 1]
        if np.any(values.ravel()[boundary] != 0):
            raise ValueError('compact P1 trace must be exactly zero on every perimeter edge')
        panels = TrianglePanels(vertices[faces], 1)
        nodal = values.ravel()[faces]
        gradients = np.linalg.solve(panels.triangles[:, 1:, :2]-panels.triangles[:, :1, :2],
                                    (nodal[:, 1:]-nodal[:, :1])[..., None])[..., 0]
        self.x, self.y, self.values = x, y, values
        self.panels, self.faces, self.boundary_edges = panels, faces, boundary
        self._nodal, self._gradients = nodal, gradients
        for array in (x, y, values, faces, boundary, nodal, gradients):
            array.setflags(write=False)

    @staticmethod
    def _points(points, dimensions):
        raw = np.asarray(points)
        if raw.ndim != 2 or raw.shape[1] != dimensions or np.iscomplexobj(raw) or not np.isfinite(raw).all():
            raise ValueError(f'finite real (N,{dimensions}) targets required')
        return np.asarray(raw, dtype=np.float64)

    def bulk(self, points):
        """Return potential (N,) and velocity (N,3), strictly below z=0.

        Existing near-panel accuracy guards may reject unresolved targets; no
        fixed quadrature acceptance or displacement away from a boundary.
        """
        points = self._points(points, 3)
        if np.any(points[:, 2] >= 0):
            raise ValueError('bulk targets must lie strictly below the source plane')
        value, velocity = np.empty(len(points)), np.empty((len(points), 3))
        for first in range(0, len(points), 128):
            last = min(first+128, len(points))
            targets = points[first:last]
            influence = checked_influence(targets, self.panels)
            gs = influence.gradient_single_layer[..., :2]
            contraction = np.einsum('tnj,nj->tn', gs, self._gradients)
            affine = self._nodal[None, :, 0]+np.sum((targets[:, None, :2]
                     -self.panels.triangles[None, :, 0, :2])*self._gradients, axis=-1)
            value[first:last] = -2*np.sum(affine*influence.double_layer+targets[:, 2, None]*contraction, axis=1)
            velocity[first:last, :2] = -2*influence.double_layer@self._gradients
            velocity[first:last, 2] = -2*np.sum(contraction, axis=1)
        return value, velocity

    def surface_velocity(self, indices):
        """Velocity at explicitly owned triangle centroids, never nodes/edges.

        P1 gradient jumps give logarithmic DNO edge singularities. These point
        values are not a conservative nodal evolution or weak flux integral.
        """
        indices = np.asarray(indices)
        if (indices.ndim != 1 or indices.dtype.kind not in 'iu'
                or np.any(indices < 0) or np.any(indices >= len(self.panels))):
            raise ValueError('valid panel indices required')
        velocity = np.empty((len(indices), 3))
        velocity[:, :2] = self._gradients[indices]
        for first in range(0, len(indices), 128):
            last = min(first+128, len(indices))
            trace = single_layer_trace(self.panels, indices[first:last])
            velocity[first:last, 2] = -2*np.einsum('tnj,nj->t',
                trace.gradient_single_layer[..., :2], self._gradients)
        return velocity

    def outside_velocity(self, points_xy):
        """Global return velocity on z=0 strictly outside the closed patch.

        The tangential trace derivative is zero here, but the normal return
        generally is nonzero. It is never truncated at the patch boundary.
        """
        xy = self._points(points_xy, 2)
        outside = ((xy[:, 0] < self.x[0]) | (xy[:, 0] > self.x[-1])
                   | (xy[:, 1] < self.y[0]) | (xy[:, 1] > self.y[-1]))
        if not np.all(outside):
            raise ValueError('surface targets must be strictly outside the closed patch')
        velocity = np.zeros((len(xy), 3))
        for first in range(0, len(xy), 128):
            last = min(first+128, len(xy))
            points = np.column_stack((xy[first:last], np.zeros(last-first)))
            trace = checked_influence(points, self.panels)
            velocity[first:last, 2] = -2*np.einsum('tnj,nj->t',
                trace.gradient_single_layer[..., :2], self._gradients)
        return velocity

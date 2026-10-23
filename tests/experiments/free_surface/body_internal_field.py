"""Independent zero-charge fields inside a convex displaced solid.

This verification fixture owns no body BVP solver. Freeze its source positions
before comparing numerical partitions of the same physical wetted boundary.
"""
from dataclasses import dataclass
import numpy as np
from scipy.spatial import ConvexHull
from topwater.rankine import TrianglePanels, panel_influence


def single_first_moment(targets, panels):
    """Return integral G dA and integral s G dA per target/panel, G=1/(4pi r)."""
    targets = np.asarray(targets, dtype=float)
    single = panel_influence(targets, panels).single_layer
    offset = targets[:, None, :]-panels.triangles[None, :, 0]
    height = np.einsum('tfi,fi->tf', offset, panels.geometric_normals)
    projection = targets[:, None, :]-height[..., None]*panels.geometric_normals
    first = projection*single[..., None]
    for edge in range(3):
        start = panels.triangles[:, edge][None]-targets[:, None]
        tangent = panels.tangents[:, edge]
        a = np.einsum('tfi,fi->tf', start, tangent)
        perpendicular = start-a[..., None]*tangent
        b = np.linalg.norm(perpendicular, axis=-1)
        if np.any(b <= 0):
            raise ValueError('verification targets must avoid supporting edge lines')
        end = a+panels.edge_lengths[:, edge]
        def primitive(value):
            return .5*(value*np.hypot(value, b)+b*b*np.arcsinh(value/b))
        integral = primitive(end)-primitive(a)
        first += integral[..., None]*panels.edge_outward[:, edge]/(4*np.pi)
    return single, first


@dataclass(frozen=True, eq=False)
class InternalField:
    points: np.ndarray
    strengths: np.ndarray
    center: np.ndarray
    inscribed_radius: float

    def field(self, points):
        """Potential (...,2), Cartesian gradient (...,2,3), independent of BEM."""
        points = np.asarray(points, dtype=float)
        delta = points[..., None, :]-self.points
        radius = np.linalg.norm(delta, axis=-1)
        if np.any(radius == 0) or not np.isfinite(radius).all():
            raise ValueError('finite evaluation points distinct from sources required')
        potential = (1/(4*np.pi*radius))@self.strengths
        gradient = np.einsum('...si,sc->...ci', -delta/(4*np.pi*radius[..., None]**3), self.strengths)
        return potential, gradient

    def pressure_load(self, panels, cg, density):
        """All six body loads for phi_t=field, pressure=-density*phi_t.

        panels.normals MUST point out of physical fluid, into the solid. The
        returned columns correspond to the two independent source-sink pairs.
        """
        cg = np.asarray(cg, dtype=float)
        if cg.shape != (3,) or not np.isfinite(cg).all() or not np.isfinite(density) or density <= 0:
            raise ValueError('finite CG and positive density required')
        single, first = single_first_moment(self.points, panels)
        normal = panels.normals
        force = np.einsum('sf,fi->si', single, normal)
        moment = np.cross(first-cg*single[..., None], normal).sum(axis=1)
        return -density*np.vstack((force.T, moment.T))@self.strengths


def internal_field(wetted_triangles):
    """Construct fixed vertical/transverse source-sink pairs strictly in solid.

    The convex hull supplies the water cap. Its exact tetrahedron volume formula
    determines the volume centroid, not the vertex or panel-centroid average.
    """
    panels = TrianglePanels(wetted_triangles)
    vertices = np.unique(panels.triangles.reshape(-1, 3), axis=0)
    hull = ConvexHull(vertices)
    tolerance = 256*np.finfo(float).eps*max(1., np.max(abs(vertices)))
    plane_values = panels.triangles@hull.equations[:, :3].T+hull.equations[:, 3]
    if not np.all(np.any(np.max(abs(plane_values), axis=1) <= tolerance, axis=1)):
        raise ValueError('wetted triangles must lie on the convex displaced-solid boundary')
    inside = vertices.mean(axis=0)
    faces = vertices[hull.simplices]
    volume6 = np.abs(np.linalg.det(faces-inside))
    center = np.sum(volume6[:, None]*(faces.sum(axis=1)+inside)/4, axis=0)/volume6.sum()
    normals, offset = hull.equations[:, :3], hull.equations[:, 3]
    radius = np.min(-(normals@center+offset)/np.linalg.norm(normals, axis=1))
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('nondegenerate convex displaced solid required')
    axes = np.array([[0., 0., 1.], [0., 0., -1.], [0., 1., 0.], [0., -1., 0.]])
    points = center+.25*radius*axes
    if np.any(points@normals.T+offset >= 0):
        raise ValueError('all sources must lie strictly inside displaced solid')
    strengths = np.array([[1., 0.], [-1., 0.], [0., 1.], [0., -1.]])
    for array in (points, strengths, center):
        array.setflags(write=False)
    return InternalField(points, strengths, center, float(radius))

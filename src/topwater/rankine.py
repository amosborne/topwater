"""Constant-triangle positive Rankine kernels; no boundary solve or wave model.

Targets are evaluated in bounded blocks. Geometry is prepared once; S and its
gradient ignore orientation, while K and its gradient use fluid-normal signs.
"""

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from ._validation import array
from .numerics import LaplaceControls


class QuadratureError(RuntimeError):
    """A fallback rule did not meet the requested error within its work budget."""


@dataclass(frozen=True, eq=False)
class TrianglePanels:
    """Cached flat numerical panels, independent of the full closed hull class.

    triangles has shape (N,3,3). normal_sign is +1 for the triangle winding,
    -1 for its reverse, or one such sign per panel. Use -1 for fluid normals on
    a body-outward hull; use +1 for the manufactured bounded interior test.
    """

    triangles: object
    normal_sign: object = 1

    def __post_init__(self):
        triangles = array(self.triangles, (None, 3, 3), "triangles")
        if len(triangles) == 0:
            raise ValueError("at least one numerical panel is required")
        edges = np.roll(triangles, -1, axis=1) - triangles
        lengths = np.linalg.norm(edges, axis=2)
        cross = np.cross(edges[:, 0], -edges[:, 2])
        twice_area = np.linalg.norm(cross, axis=1)
        if np.any(twice_area <= 64*np.finfo(float).eps*lengths[:, 0]*lengths[:, 2]):
            raise ValueError("numerical panels must be nondegenerate triangles")
        signs = np.asarray(self.normal_sign)
        if signs.dtype.kind not in "iuf" or signs.shape not in ((), (len(triangles),)):
            raise ValueError("normal_sign must be +1/-1, scalar or one value per panel")
        if not np.all((signs == 1) | (signs == -1)):
            raise ValueError("normal_sign must contain only +1/-1")
        signs = np.broadcast_to(signs, (len(triangles),)).astype(float, copy=True)
        normals = cross / twice_area[:, None]
        tangents = edges / lengths[:, :, None]
        values = dict(
            triangles=triangles, normal_sign=signs, edge_lengths=lengths,
            tangents=tangents, geometric_normals=normals,
            edge_outward=np.cross(tangents, normals[:, None, :]),
            normals=normals*signs[:, None], areas=twice_area/2,
            diameters=lengths.max(axis=1), centroids=triangles.mean(axis=1),
        )
        for name, value in values.items():
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    @classmethod
    def from_mesh(cls, mesh, normal_sign=1):
        """Use the full closed shell, including any artificial closure.

        This constructs manufactured closed-domain problems. Physical wetted
        panels instead use TrianglePanels(result.wetted_triangles, -1).
        """
        return cls(mesh.vertices[mesh.faces], normal_sign)

    def __len__(self):
        return len(self.triangles)


@dataclass(frozen=True, eq=False)
class PanelInfluence:
    single_layer: np.ndarray
    double_layer: np.ndarray
    gradient_single_layer: np.ndarray | None = None
    gradient_double_layer: np.ndarray | None = None
    fallback_pairs: int = 0
    quadrature_evaluations: int = 0


def _fallback(target, triangle, normal, controls, gradients, single_trace=False):
    """Resolve analytic cancellation with bounded additional work.

    Fixed-origin Gaussian rules can agree while missing a near-singular peak.
    Never use their order agreement as evidence for an unresolved close pair.
    """
    area2 = np.linalg.norm(np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0]))
    diameter = np.linalg.norm(np.roll(triangle, -1, axis=0)-triangle, axis=1).max()
    if np.linalg.norm(target-triangle.mean(axis=0)) < 2*diameter:
        if gradients:
            raise QuadratureError("analytic accuracy unresolved for a close panel gradient")
        # A close peak can defeat quadrature order agreement. Re-evaluate the
        # analytic integral instead; edge/vertex/trace checks run before this.
        from decimal import DecimalException
        from ._rankine_extended import triangle_values

        previous = None
        last_error = None
        for refinement in range(min(controls.max_refinement+1, 5)):
            try:
                values, rounding = triangle_values(target, triangle, normal, 32+16*refinement)
            except (DecimalException, ValueError, ZeroDivisionError) as error:
                # Cancellation can round a positive denominator to zero at a
                # lower precision. Retry within the same bounded work budget.
                last_error = error
                previous = None
                continue
            scaled = np.asarray(values)/[diameter, 1.]
            error = np.asarray(rounding)/[diameter, 1.]
            tolerance = controls.quadrature_atol+controls.quadrature_rtol*np.abs(scaled)
            if (previous is not None and np.all(np.isfinite(scaled)) and values[0] > 0
                    and np.all(np.abs(scaled-previous)+error <= tolerance)):
                return values, 0
            previous = scaled
        raise QuadratureError("close panel analytic evaluation exceeded max_refinement or float precision") from last_error
    previous = None
    evaluations = 0
    for refinement in range(controls.max_refinement + 1):
        order = 4*(refinement + 1)
        nodes, weights = leggauss(order)
        u, v = np.meshgrid((nodes+1)/2, (nodes+1)/2, indexing="ij")
        y = triangle[0] + u[:, :, None] * (
            (1-v[:, :, None])*(triangle[1]-triangle[0])
            + v[:, :, None]*(triangle[2]-triangle[0])
        )
        delta = target-y
        distance = np.linalg.norm(delta, axis=2)
        measure = np.outer(weights, weights) * u * area2 / (16*np.pi)
        height = delta @ normal
        result = [np.sum(measure/distance), np.sum(measure*height/distance**3)]
        if gradients:
            result.extend(np.sum(-measure[:, :, None]*delta/distance[:, :, None]**3, axis=(0, 1)))
            if not single_trace:
                result.extend(np.sum(measure[:, :, None]*(normal/distance[:, :, None]**3
                               - 3*height[:, :, None]*delta/distance[:, :, None]**5), axis=(0, 1)))
        result = np.asarray(result)
        evaluations += order*order
        scaled = result.copy()
        scaled[0] /= diameter
        if gradients and not single_trace:
            scaled[5:8] *= diameter
        if previous is not None and np.all(
            np.abs(scaled-previous) <= controls.quadrature_atol + controls.quadrature_rtol*np.abs(scaled)
        ):
            return result, evaluations
        previous = scaled
    raise QuadratureError("panel quadrature failed to converge within max_refinement")


def panel_influence(targets, panels, *, controls=None, gradients=False, self_indices=None):
    """Return S, principal-value K and optional bulk gradients for one block.

    Explicit self_indices identify centroid principal values. Edges, vertices,
    boundary gradients and unresolved near-panel cancellation fail explicitly.
    """
    return _evaluate_block(targets, panels, controls=controls, gradients=gradients,
                           self_indices=self_indices)


def single_layer_trace(panels, indices, *, controls=None):
    """S, principal-value K and domain-side grad(S) at indexed centroids.

    The domain lies toward -panels.normals. Only the owned self pair receives
    the +0.5*normal jump; no double-layer gradient is requested or evaluated.
    Targets must fit controls.target_block_size, as for panel_influence.
    """
    if not isinstance(panels, TrianglePanels):
        raise ValueError('expected cached TrianglePanels')
    indices = np.asarray(indices)
    if (indices.ndim != 1 or indices.dtype.kind not in 'iu'
            or np.any(indices < 0) or np.any(indices >= len(panels))):
        raise ValueError('indices must contain valid source-panel indices')
    result = _evaluate_block(panels.centroids[indices], panels, controls=controls,
                             gradients=True, self_indices=indices, single_trace=True)
    result.gradient_single_layer[np.arange(len(indices)), indices] += .5*panels.normals[indices]
    return result


def _evaluate_block(targets, panels, *, controls=None, gradients=False, self_indices=None,
                    outputs=None, single_trace=False):
    """Evaluate one target block against cached panels.

    S and K are (M,N), with no trace jump in K. Optional bulk gradients are
    (M,N,3). Supply self_indices (source index or -1 per target) for centroid
    principal values; boundary gradients and edge/vertex targets are rejected.
    An assembly caller must iterate blocks instead of requesting global NN3 data.
    """
    controls = LaplaceControls() if controls is None else controls
    if not isinstance(controls, LaplaceControls) or not isinstance(panels, TrianglePanels):
        raise ValueError("expected LaplaceControls and cached TrianglePanels")
    targets = array(targets, (None, 3), "targets")
    if len(targets) > controls.target_block_size:
        raise ValueError("target count exceeds target_block_size; evaluate bounded blocks")
    indices = np.full(len(targets), -1, dtype=np.int64)
    if self_indices is not None:
        indices = np.asarray(self_indices)
        if indices.shape != (len(targets),) or indices.dtype.kind not in "iu":
            raise ValueError("self_indices must be an integer index or -1 per target")
        if np.any((indices < -1) | (indices >= len(panels))):
            raise ValueError("self_indices contains an invalid panel index")
        rows = np.flatnonzero(indices >= 0)
        columns = indices[rows]
        error = np.linalg.norm(targets[rows]-panels.centroids[columns], axis=1)
        if np.any(error > 64*np.finfo(float).eps*panels.diameters[columns]):
            raise ValueError("self_indices must identify centroid targets")

    from ._rankine_kernel import analytic_block

    shape = (len(targets), len(panels))
    single, double = ((np.empty(shape), np.empty(shape)) if outputs is None else outputs)
    gradient_single = np.empty((*shape, 3 if gradients else 0))
    gradient_double = np.empty((*shape, 3 if gradients and not single_trace else 0))
    unsafe, bad = analytic_block(
        targets, panels.triangles, panels.geometric_normals, panels.edge_outward,
        panels.tangents, panels.edge_lengths, panels.areas, panels.diameters,
        panels.normal_sign, indices, gradients, controls.quadrature_atol,
        controls.quadrature_rtol, single, double, gradient_single, gradient_double, single_trace,
    )
    if bad & 1:
        raise ValueError("targets on panel vertices are unsupported")
    if bad & 2:
        raise ValueError("targets on panel edges are unsupported")
    if bad & 4:
        raise ValueError("boundary centroid targets require explicit self_indices")
    if bad & 8:
        raise ValueError("boundary gradients require separate trace limits")
    fallbacks = 0
    evaluations = 0
    for row, column in zip(*np.nonzero(unsafe)):
        if indices[row] == column:
            raise QuadratureError("invalid analytic self-panel integral")
        value, count = _fallback(targets[row], panels.triangles[column], panels.normals[column], controls, gradients, single_trace)
        single[row, column], double[row, column] = value[:2]
        if gradients:
            gradient_single[row, column] = value[2:5]
            if not single_trace:
                gradient_double[row, column] = value[5:8]
        fallbacks += 1
        evaluations += count
    return PanelInfluence(single, double, gradient_single if gradients else None,
                          gradient_double if gradients and not single_trace else None, fallbacks, evaluations)

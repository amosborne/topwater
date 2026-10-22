"""Physical-facet coverage for the contained canonical-canoe radiation case.

Different numerical subdivisions must cover the same clipped planar facets.
This does not waive the pinned 704-panel Capytaine hash or prove field accuracy.
Candidate triangles use BODY-outward orientation and original hull face IDs.
"""
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

from .moving_canoe import make_domain, geometry
from .p4_radiation import read_reference


def pinned_reference(path, frequencies=(2., 4., 6.)):
    """Regenerate the ORIGINAL wet2 contact mesh, without TetGen/fluid solves."""
    domain = make_domain(12, False, 2)
    points, _, _ = geometry(domain, domain.rest.position, domain.rest.orientation,
                            np.zeros(len(domain.free_nodes)))
    triangles = points[domain.local_mapping][domain.contact.faces[domain.contact.body_faces]][:, ::-1]
    digest = sha256(np.asarray(triangles, dtype='<f8').tobytes()).hexdigest()
    matrices, settings = read_reference(Path(path), digest, domain.rest.position, frequencies)
    return domain, triangles, matrices, dict(wetted_triangles_sha256=digest,
        wetted_panels=len(triangles), rotation_center=domain.rest.position.tolist(),
        water_depth=4., density=1000., gravity=9.81,
        physical_constants_provenance='Explicit rho=1000,g=9.81 in tests.capytaine_reference.frequency_response; '
            'legacy reference report omits these two metadata fields.',
        capytaine_source_sha256=sha256(Path(__file__).parents[2].joinpath(
            'capytaine_reference.py').read_bytes()).hexdigest(),
        capytaine_green_function_settings=settings)


def _cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def _area(polygon):
    if len(polygon) < 3:
        return 0.
    return float(abs(np.sum(_cross(polygon, np.roll(polygon, -1, axis=0))))/2)


def _intersection_area(first, second):
    """Convex triangle intersection, including shared edges of zero area."""
    polygon = np.asarray(first)
    for a, b in zip(second, np.roll(second, -1, axis=0)):
        if not len(polygon):
            return 0.
        output = []
        previous = polygon[-1]
        old = _cross(b-a, previous-a)
        for current in polygon:
            new = _cross(b-a, current-a)
            if (old >= 0) != (new >= 0):
                output.append(previous+old/(old-new)*(current-previous))
            if new >= 0:
                output.append(current)
            previous, old = current, new
        polygon = np.asarray(output)
    return _area(polygon)


def verify_coverage(triangles, parents, reference, reference_parents):
    """Containment + per-parent area + pairwise overlap at roundoff tolerance.

    Reference parent regions must be convex (flat-water clipping of an original
    triangle). Explicit overlap checks prevent equal-area holes/duplicates from
    passing. Tolerances bound numerical geometry checks, not hydrodynamic error.
    """
    triangles, reference = np.asarray(triangles), np.asarray(reference)
    parents, reference_parents = np.asarray(parents), np.asarray(reference_parents)
    for values, ids in ((triangles, parents), (reference, reference_parents)):
        if values.ndim != 3 or values.shape[1:] != (3, 3) or not len(values) or ids.shape != (len(values),) or not np.all(np.isfinite(values)):
            raise ValueError('finite triangles with one original parent ID each required')
    if not np.array_equal(np.unique(parents), np.unique(reference_parents)):
        raise ValueError('candidate loses or adds a physical parent facet')
    records = []
    for parent in np.unique(reference_parents):
        target = reference[reference_parents == parent]
        candidate = triangles[parents == parent]
        origin = target[0, 0]
        edge = target[0, 1]-origin
        normal = np.cross(edge, target[0, 2]-origin)
        normal /= np.linalg.norm(normal)
        axis = edge/np.linalg.norm(edge)
        frame = np.array([axis, np.cross(normal, axis)])
        projected = (target-origin)@frame.T
        hull = ConvexHull(np.unique(projected.reshape(-1, 2), axis=0))
        expected = _area(hull.points[hull.vertices])
        scale = max(1., float(np.max(abs(target))), float(np.ptp(target.reshape(-1, 3), axis=0).max()))
        length_tolerance = 256*np.finfo(float).eps*scale
        area_tolerance = 1024*np.finfo(float).eps*scale**2*len(candidate)
        planar_error = float(np.max(abs((candidate-origin)@normal)))
        local = (candidate-origin)@frame.T
        outside = float(np.max(local@hull.equations[:, :2].T+hull.equations[:, 2]))
        signed = _cross(local[:, 1]-local[:, 0], local[:, 2]-local[:, 0])/2
        if planar_error > length_tolerance or outside > length_tolerance or np.any(signed <= 0):
            raise ValueError(f'candidate is outside or reverses physical parent {parent}')
        area_error = float(abs(signed.sum()-expected))
        overlap = sum(_intersection_area(a, b) for i, a in enumerate(local) for b in local[i+1:])
        if area_error > area_tolerance or overlap > area_tolerance:
            raise ValueError(f'gap or overlap in physical parent {parent}: area error {area_error}, overlap {overlap}')
        records.append(dict(parent=int(parent), triangles=len(candidate), reference_area=expected,
            absolute_area_error=area_error, overlapping_area=float(overlap),
            maximum_plane_error=planar_error, maximum_outside_distance=outside,
            length_tolerance=length_tolerance, area_tolerance=area_tolerance))
    return dict(candidate_triangles=len(triangles), physical_parents=len(records),
                parent_coverage=records, field_accuracy_proven=False)

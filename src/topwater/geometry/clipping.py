"""Flat-water hull clipping and exact triangle pressure/volume integrals."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def clip_flat_water(vertices, faces, physical_faces, level, save_geometry):
    """Vertices use world axes relative to CG; an artificial cap is never read."""
    capacity = 2*len(physical_faces) if save_geometry else 0
    wet = np.empty((capacity, 3, 3))
    source = np.empty(capacity, dtype=np.int64)
    segments = np.empty((len(physical_faces) if save_geometry else 0, 2, 3))
    polygon = np.empty((4, 3))
    force, moment, first = np.zeros(3), np.zeros(3), np.zeros(3)
    volume = area = 0.
    count = segment_count = 0
    for face in physical_faces:
        size = 0
        previous = faces[face, 2]
        previous_inside = vertices[previous, 2] < level
        for corner in range(3):
            current = faces[face, corner]
            current_inside = vertices[current, 2] < level
            if current_inside != previous_inside:
                # Incident faces traverse a shared edge in opposite directions.
                # Always interpolate from wet to dry so the result is identical.
                wet_vertex = current if current_inside else previous
                dry_vertex = previous if current_inside else current
                if vertices[dry_vertex, 2] == level:
                    polygon[size] = vertices[dry_vertex]
                else:
                    fraction = ((level-vertices[wet_vertex, 2])
                                / (vertices[dry_vertex, 2]-vertices[wet_vertex, 2]))
                    for axis in range(3):
                        polygon[size, axis] = (vertices[wet_vertex, axis]
                            + fraction*(vertices[dry_vertex, axis]-vertices[wet_vertex, axis]))
                polygon[size, 2] = level
                size += 1
            if current_inside:
                polygon[size] = vertices[current]
                size += 1
            previous, previous_inside = current, current_inside
        if size < 3:
            continue
        if save_geometry:
            for edge in range(size):
                following = (edge+1) % size
                if polygon[edge, 2] == level and polygon[following, 2] == level:
                    if np.any(polygon[edge] != polygon[following]):
                        segments[segment_count, 0] = polygon[edge]
                        segments[segment_count, 1] = polygon[following]
                        segment_count += 1
        for corner in range(1, size-1):
            p, q, r = polygon[0], polygon[corner], polygon[corner+1]
            # Relative-plane determinant needs no triangle copies or pressure
            # geometry. Root iterations execute only this volume calculation.
            pz, qz, rz = p[2]-level, q[2]-level, r[2]-level
            dv = (p[0]*(q[1]*rz-qz*r[1])
                  -p[1]*(q[0]*rz-qz*r[0])
                  +pz*(q[0]*r[1]-q[1]*r[0]))/6
            if not save_geometry:
                volume += dv
                continue
            a, b = q-p, r-p
            normal_area = np.cross(a, b)/2
            magnitude = np.sqrt(np.dot(normal_area, normal_area))
            if magnitude == 0:
                continue
            # Close the wetted surface implicitly at z=level. Tetrahedra from
            # that plane have zero contribution from its artificial closure.
            volume += dv
            for axis in range(3):
                first[axis] += dv*(p[axis]+q[axis]+r[axis]
                                  -(3*level if axis == 2 else 0))/4
            depths = np.array([level-p[2], level-q[2], level-r[2]])
            depth_sum = depths.sum()
            force -= (depth_sum/3)*normal_area
            weighted = ((p+q+r)*depth_sum
                        + p*depths[0]+q*depths[1]+r*depths[2])/12
            moment -= np.cross(weighted, normal_area)
            area += magnitude
            wet[count, 0], wet[count, 1], wet[count, 2] = p, q, r
            source[count] = face
            count += 1
    return wet[:count], source[:count], segments[:segment_count], volume, first, force, moment, area

"""Independent quadrature reference, never a runtime dependency or fast path.

Fan triangles about the target's plane projection have analytically integrated
radial coordinates. SciPy integrates the remaining coordinate. Gradients use
separate edge quadrature and a radial normal-derivative integral.
"""

import numpy as np
from scipy.integrate import quad


def panel_reference(target, triangle, normal_sign=1):
    target, triangle = np.asarray(target), np.asarray(triangle)
    cross = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
    normal = cross/np.linalg.norm(cross)
    height = np.dot(target-triangle[0], normal)
    projection = target-height*normal
    exterior_projection = np.any(np.cross(np.roll(triangle, -1, axis=0)-triangle,
                                          projection-triangle) @ normal < 0)
    single = double = grad_normal_double = 0.0
    gradient_single = np.zeros(3)
    gradient_double = np.zeros(3)
    for i in range(3):
        start, end = triangle[i], triangle[(i+1) % 3]
        left, right = start-projection, end-projection
        signed_area2 = np.dot(np.cross(left, right), normal)
        direction = right-left
        closest = -np.dot(left, direction)/np.dot(direction, direction)
        points = [closest] if 0 < closest < 1 else None

        def radius(v):
            ray = left+v*direction
            return np.sqrt(height*height+np.dot(ray, ray))

        def integrate(function):
            return quad(function, 0, 1, epsabs=2e-13, epsrel=2e-13,
                        limit=300, points=points)[0]

        if signed_area2 != 0:
            single += signed_area2 * integrate(lambda v: 1/(radius(v)+abs(height)))
            if height != 0:
                double += signed_area2*np.sign(height) * integrate(
                    lambda v: 1/(radius(v)*(radius(v)+abs(height))))
            # At coplanar exterior targets K vanishes, but its normal derivative
            # does not. Self-panel gradients are outside this reference's API.
            if height != 0 or exterior_projection:
                grad_normal_double -= signed_area2 * integrate(lambda v: 1/radius(v)**3)
        tangent = (end-start)/np.linalg.norm(end-start)
        outward = np.cross(tangent, normal)
        edge_length = np.linalg.norm(end-start)

        def distance(v):
            delta = target-(start+v*(end-start))
            return np.linalg.norm(delta)

        gradient_single -= outward*edge_length * quad(
            lambda v: 1/distance(v), 0, 1, epsabs=2e-13, epsrel=2e-13, limit=300)[0]
        if height != 0:
            gradient_double -= height*outward*edge_length * quad(
                lambda v: 1/distance(v)**3, 0, 1, epsabs=2e-13, epsrel=2e-13, limit=300)[0]
    gradient_single -= normal*double
    gradient_double += normal*grad_normal_double
    return (single/(4*np.pi), normal_sign*double/(4*np.pi),
            gradient_single/(4*np.pi), normal_sign*gradient_double/(4*np.pi))

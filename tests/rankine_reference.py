"""Independent quadrature reference, never a runtime dependency or fast path.

Fan triangles about the target's plane projection have analytically integrated
radial coordinates. SciPy integrates the remaining coordinate. Gradients use
separate edge quadrature and a radial normal-derivative integral.
"""

import numpy as np
from scipy.integrate import quad


def edge_aligned_values(target, triangle, normal_sign=1):
    """Independent off-plane S/K reference for thin waterline triangles.

    Integrate along the longest edge analytically, then adapt across its width.
    Endpoint-transition splits resolve the narrow peak which global Gaussian
    refinement misses. This is not a trace or gradient reference.
    """
    triangle = np.asarray(triangle)
    lengths = np.linalg.norm(np.roll(triangle, -1, axis=0)-triangle, axis=1)
    a, b, c = np.roll(triangle, -int(np.argmax(lengths)), axis=0)
    length = np.linalg.norm(b-a)
    tangent = (b-a)/length
    normal = np.cross(b-a, c-a)
    normal /= np.linalg.norm(normal)
    short = np.cross(normal, tangent)
    width, skew = np.dot(c-a, short), np.dot(c-a, tangent)
    along, across = np.dot(target-a, tangent), np.dot(target-a, short)
    height = np.dot(target-a, normal)
    if height == 0:
        raise ValueError('edge-aligned reference requires an off-plane target')
    split = []
    if skew != 0:
        split.append(along/skew)
    if skew != length:
        split.append((along-length)/(skew-length))
    split = [v for v in split if 0 < v < 1]

    def integrand(v, double):
        lo = v*skew-along
        hi = lo+(1-v)*length
        perpendicular2 = (v*width-across)**2+height**2
        if double:
            radial = (hi/np.sqrt(hi*hi+perpendicular2)
                      - lo/np.sqrt(lo*lo+perpendicular2))/perpendicular2
            return normal_sign*width*height*radial/(4*np.pi)
        radial = (np.arcsinh(hi/np.sqrt(perpendicular2))
                  - np.arcsinh(lo/np.sqrt(perpendicular2)))
        return width*radial/(4*np.pi)

    values = []
    for double in (False, True):
        result = quad(integrand, 0, 1, args=(double,), points=split,
                      epsabs=1e-13 if double else 1e-17,
                      epsrel=1e-12, limit=300, full_output=1)
        if len(result) != 3:
            raise RuntimeError(result[3])
        values.append(result[0])
    return tuple(values)


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

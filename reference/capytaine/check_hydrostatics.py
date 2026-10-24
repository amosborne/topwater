"""Analytic polynomial-moment regression for the hydrostatics-only correction.

No BEM solve. Execute only in the coordinator's allocated correction/check job.
"""
import argparse
import itertools
import json
import os
from pathlib import Path

for thread_variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[thread_variable] = '1'

import numpy as np
import capytaine as cpt
from scipy.io import loadmat


def moments(triangle):
    """Exact area, integral r, integral rr' for a planar linear triangle."""
    a, b, c = triangle
    area = np.linalg.norm(np.cross(b-a, c-a))/2
    total = triangle.sum(axis=0)
    return area, area*total/3, area*(np.outer(total, total)+triangle.T@triangle)/12


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corrected', type=Path, required=True)
    args = parser.parse_args()
    unit = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    skew = np.array([[-1.2, .4, -2.], [2.3, -.7, -.3], [.2, 1.5, -.8]])
    rows = []
    for name, triangle in [('unit_right', unit), ('skew_translated', skew)]:
        for permutation in itertools.permutations(range(3)):
            triangle_order = triangle[list(permutation)]
            mesh = cpt.Mesh(vertices=triangle_order, faces=[[0, 1, 2]],
                            quadrature_method='Gauss-Legendre 2', auto_clean=False, auto_check=False)
            points, weights = mesh.quadrature_points
            expected = moments(triangle_order)
            actual = (weights.sum(), np.einsum('fqi,fq->i', points, weights),
                      np.einsum('fqi,fqj,fq->ij', points, points, weights))
            errors = [float(np.max(np.abs(a-b))) for a, b in zip(actual, expected)]
            passed = all(np.allclose(a, b, rtol=2e-13, atol=2e-13) for a, b in zip(actual, expected))
            rows.append(dict(triangle=name, vertex_permutation=list(permutation),
                             absolute_errors=errors, passed=passed))
    centroid = cpt.Mesh(vertices=unit, faces=[[0, 1, 2]], auto_clean=False, auto_check=False)
    p, w = centroid.quadrature_points
    centroid_x_squared = float(np.sum(p[..., 0]**2*w))
    exact_x_squared = 1/12
    original_rule_rejected = abs(centroid_x_squared-exact_x_squared) > 1e-3
    # Independent closed-surface polynomial moments on the actual unchanged hull.
    with np.load(args.corrected/'native.npz') as raw:
        triangles = raw['triangles']
    mat = loadmat(args.corrected/'native.mat')
    volume_integral = 0.
    second_normal_integral = np.zeros(3)
    waterplane_area = 0.
    waterplane_first = np.zeros(3)
    waterplane_second = np.zeros((3, 3))
    for triangle in triangles:
        area, first, second = moments(triangle)
        cross = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
        normal = cross/np.linalg.norm(cross)
        volume_integral += np.dot(normal, first)/3
        second_normal_integral += normal*np.diag(second)
        waterplane_area -= normal[2]*area
        waterplane_first -= normal[2]*first
        waterplane_second -= normal[2]*second
    cb = second_normal_integral/(2*volume_integral)
    rho = mat['rho'].item()
    gravity = mat['g'].item()
    mass = mat['M'][0, 0]
    cg = mat['cg'].ravel()
    xc, yc, zc = cg
    normalized = np.zeros((6, 6))
    normalized[2, 2] = waterplane_area
    normalized[2, 3] = normalized[3, 2] = waterplane_first[1]-yc*waterplane_area
    normalized[2, 4] = normalized[4, 2] = -(waterplane_first[0]-xc*waterplane_area)
    normalized[3, 3] = (waterplane_second[1, 1]-2*yc*waterplane_first[1]+yc**2*waterplane_area
                         + volume_integral*(cb[2]-zc)-mass/rho*(cg[2]-zc))
    normalized[4, 4] = (waterplane_second[0, 0]-2*xc*waterplane_first[0]+xc**2*waterplane_area
                         + volume_integral*(cb[2]-zc)-mass/rho*(cg[2]-zc))
    normalized[3, 4] = normalized[4, 3] = -(waterplane_second[0, 1]-xc*waterplane_first[1]
                                          -yc*waterplane_first[0]+xc*yc*waterplane_area)
    normalized[3, 5] = -volume_integral*(cb[0]-xc)+mass/rho*(cg[0]-xc)
    normalized[4, 5] = -volume_integral*(cb[1]-yc)+mass/rho*(cg[1]-yc)
    stiffness = rho*gravity*normalized
    checks = dict(triangle_moments=all(row['passed'] for row in rows),
                  centroid_rule_rejected=original_rule_rejected,
                  native_volume_matches_analytic=bool(np.allclose(mat['volume'].item(), volume_integral, rtol=2e-13, atol=2e-13)),
                  native_cb_matches_analytic=bool(np.allclose(mat['cb'].ravel(), cb, rtol=2e-13, atol=2e-13)),
                  native_K_matches_analytic=bool(np.allclose(mat['K'], stiffness, rtol=2e-13, atol=2e-8)),
                  lateral_buoyancy_symmetry=bool(abs(cb[1]) < 2e-12),
                  physical_equilibrium_cb_xy=bool(np.max(np.abs(cb[:2]-cg[:2])) < 2e-10),
                  reflection_odd_K_terms_zero=bool(np.max(np.abs(stiffness[[2, 3, 4], [3, 4, 5]])) < 2e-8))
    result = dict(checks=checks, all_passed=all(checks.values()), triangle_regressions=rows,
                  unit_triangle_centroid_x_squared=centroid_x_squared,
                  unit_triangle_exact_x_squared=exact_x_squared, analytic_volume_m3=float(volume_integral),
                  analytic_cb_m=cb.tolist(), analytic_K=stiffness.tolist(),
                  native_K_max_absolute_error=float(np.max(np.abs(mat['K']-stiffness))),
                  derivation='Triangle area A, s=a+b+c, Q=aa^T+bb^T+cc^T: integral r=A*s/3, integral rr^T=A*(s*s^T+Q)/12. Divergence theorem gives V, CB; -n_z-weighted surface moments give waterplane moments.')
    (args.corrected/'hydrostatics-checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('checks', 'all_passed', 'native_K_max_absolute_error')}, indent=2))
    return 0 if result['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

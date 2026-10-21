"""Independent checks of nonzero-surface/body-motion coupling on the canoe.

Uniform flows prescribe matching outer-wall flux. They are manufactured mixed
boundary problems, not uniform flows through impermeable walls or evolving waves.
"""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from threadpoolctl import threadpool_limits

from tests.laplace_reference import triangle_samples


def check_flow(mesh, panels, body, free, center, energy, load, potential, normal_load):
    """Check the four solved columns phi=x,y,z,1; no assembly/factorization."""
    exact = np.column_stack((mesh.vertices, np.ones(len(mesh.vertices))))
    fixed = np.unique(mesh.faces[free])
    # Remove ALL prescribed body/wall flux, including shared surface-node rows.
    flux = energy[fixed]@potential-normal_load[fixed]
    local = np.full(len(mesh.vertices), -1, dtype=int)
    local[fixed] = np.arange(len(fixed))
    faces = local[mesh.faces[free]]
    mass = np.zeros((len(fixed), len(fixed)))
    for face, area in zip(faces, panels.areas[free]):
        mass[np.ix_(face, face)] += area*(np.ones((3, 3))+np.eye(3))/12
    elevation_rate = solve(mass, flux, assume_a='pos')
    expected_rate = np.tile([0., 0., 1., 0.], (len(fixed), 1))
    error = elevation_rate-expected_rate
    area = panels.areas[free].sum()
    # Whole-surface continuous P1 L2 error, not an interior-node crop.
    rate_rms = np.sqrt(np.einsum('ik,ij,jk->k', error, mass, error)/area)
    # Independent integral(phi*Q); prescribed pressure p=-phi gives its negative.
    # This is not Bernoulli pressure of uniform flow. Quadrature is exact here.
    samples, weights = triangle_samples(panels.triangles[body], 3)
    values = np.concatenate((samples, np.ones((*samples.shape[:2], 1))), axis=2)
    q = np.concatenate((np.broadcast_to(panels.normals[body, None], samples.shape),
                        np.cross(samples-center, panels.normals[body, None])), axis=2)
    reference = np.einsum('f,q,fqk,fqj->jk', panels.areas[body], weights, values, q)
    actual = load.T@potential
    length = 3.5
    scale = np.array([1., 1., 1., length, length, length])[:, None]*length**3
    return dict(fields=['x', 'y', 'z', '1'],
        potential_max_length_normalized=(np.max(np.abs(potential-exact), axis=0)/length).tolist(),
        whole_surface_rate_rms=rate_rms.tolist(), integrated_flux=flux.sum(axis=0).tolist(),
        expected_integrated_flux=(area*np.array([0., 0., 1., 0.])).tolist(),
        potential_load=actual.tolist(), reference_potential_load=reference.tolist(),
        exact_trace_load_max_length_normalized=float(np.max(np.abs((load.T@exact-reference)/scale))),
        solved_load_max_length_normalized=float(np.max(np.abs((actual-reference)/scale))))


if __name__ == '__main__':
    from topwater import SimulationInputs
    from topwater.hydrostatics import equilibrate_heave
    from tests.reference_hull import input_data
    from .canoe import fluid_domain, solve_inertia

    data = SimulationInputs.from_json(json.dumps(input_data(4)))
    rows = []
    with threadpool_limits(1):
        for resolution in ((8, 4, 4), (16, 8, 4)):
            rest, wet = equilibrate_heave(data.body, data.initial_state)
            mesh, body, free = fluid_domain(wet, resolution)
            _, _, diagnostics = solve_inertia(mesh, body, free, rest.position)
            rows.append(dict(qualified=False, resolution=resolution, **diagnostics['coupling'],
                scope='manufactured prescribed wall flux; not nonlinear pressure or wave evolution'))
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-coupling.json').write_text(
        json.dumps(rows, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(rows, indent=2))

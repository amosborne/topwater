"""Current mixed BEM field and fixed-world surface rates; diagnostic tank only.

No pressure time derivative, trajectory, or open-ocean coupling is implied.
The fine surface flux is recovered before nonlinear Bernoulli projection.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from topwater import HullMesh, sea_states
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.reference_boat import GRAVITY, WATER_DENSITY
from .canoe import body_load
from .incident_baseline_benchmark import coefficients, fused_cpu
from .p4_fluid_ordering import fresh


def surface_rates(points, free_faces, free_nodes, parent_map, parent_psi, flux, *, gravity=GRAVITY):
    """Project physical kinematic/Bernoulli rates onto exposed fixed P1 parents.

    parent_map is the fine-free-node interpolation of fixed world parent state.
    flux is the weak nodal normal reaction after prescribed hull/wall flux.
    Resolve this flux on the FINE cut mesh before evaluating velocity squared;
    projecting first would hide unresolved radiation in nonlinear Bernoulli.
    """
    local = np.full(len(points), -1, dtype=int)
    local[free_nodes] = np.arange(len(free_nodes))
    faces = local[free_faces]
    if np.any(faces < 0):
        raise ValueError('free faces reference unsupported nodes')
    triangles = points[free_faces]
    edges = triangles[:, 1:, :2]-triangles[:, :1, :2]
    determinant = np.linalg.det(edges)
    if np.any(determinant <= 0):
        raise ValueError('surface projection must be upward and nondegenerate')
    area = determinant/2
    opposite = np.roll(triangles[:, :, :2], -2, axis=1)-np.roll(triangles[:, :, :2], -1, axis=1)
    gradient = np.stack((-opposite[:, :, 1], opposite[:, :, 0]), axis=-1)/determinant[:, None, None]
    slope = np.einsum('fci,fc->fi', gradient, triangles[:, :, 2])
    psi = parent_map@parent_psi
    tangent = np.einsum('fci,fc->fi', gradient, psi[faces])
    mass = np.zeros((len(free_nodes), len(free_nodes)))
    unit = (np.ones((3, 3))+np.eye(3))/12
    for i in range(3):
        for j in range(3):
            np.add.at(mass, (faces[:, i], faces[:, j]), area*unit[i, j])
    reduced_mass = parent_map.T@mass@parent_map
    fine_condition = float(np.linalg.cond(mass))
    parent_condition = float(np.linalg.cond(reduced_mass))
    # Algebraic roundoff budget, not a physical discretization-error tolerance.
    if (not np.isfinite([fine_condition, parent_condition]).all()
            or max(fine_condition, parent_condition)*np.finfo(float).eps > 1e-8):
        raise ValueError('exposed surface mass exceeds 1e-8 conditioning/roundoff budget')
    fine_factor = cho_factor(mass)
    physical_flux = cho_solve(fine_factor, flux)
    parent_factor = cho_factor(reduced_mass)
    eta_rate = cho_solve(parent_factor, parent_map.T@flux)
    unresolved = physical_flux-parent_map@eta_rate
    shape, weights = triangle_samples(np.eye(3)[None], 3)
    shape = shape[0]
    normal = np.einsum('qc,fc->fq', shape, physical_flux[faces])
    uz = (normal+np.sum(slope*tangent, axis=1)[:, None])/(1+np.sum(slope*slope, axis=1))[:, None]
    z = np.einsum('qc,fc->fq', shape, triangles[:, :, 2])
    bernoulli = (-gravity*z-.5*np.sum(tangent*tangent, axis=1)[:, None]
                 +.5*(1+np.sum(slope*slope, axis=1))[:, None]*uz**2)
    load = np.zeros(len(free_nodes))
    for corner in range(3):
        np.add.at(load, faces[:, corner], area*np.einsum('q,q,fq->f', weights, shape[:, corner], bernoulli))
    psi_rate = cho_solve(parent_factor, parent_map.T@load)
    flux_norm = float(physical_flux@mass@physical_flux)
    if not np.isfinite(eta_rate).all() or not np.isfinite(psi_rate).all():
        raise FloatingPointError('nonfinite fixed-world surface rates')
    diagnostics = dict(fine_mass_condition=fine_condition,
        parent_mass_condition=parent_condition, conditioning_roundoff_budget=1e-8,
        unresolved_kinematic_mass_relative=float(np.sqrt(max(0., unresolved@mass@unresolved)/flux_norm)) if flux_norm else 0.,
        maximum_parent_height_rate=float(np.max(abs(eta_rate))),
        maximum_parent_potential_rate=float(np.max(abs(psi_rate))))
    return eta_rate, psi_rate, diagnostics


def current_field(arrays, energy, values, twist, *, acceleration_columns=True, time=0.):
    """One current mixed factor; six responses computed only when requested."""
    points, faces, fixed = arrays['points'], arrays['faces'], arrays['free_nodes']
    mesh = HullMesh(points, faces)
    panels = TrianglePanels(points[faces], 1)
    load = body_load(mesh, panels, arrays['body_faces'], arrays['pose_position'])
    mapping = np.zeros((len(fixed), len(arrays['parent_nodes'])))
    for column in range(3):
        ids = arrays['free_parent_indices'][:, column]
        valid = ids >= 0
        mapping[np.flatnonzero(valid), ids[valid]] += arrays['free_parent_weights'][valid, column]
    if not np.allclose(mapping.sum(axis=1), 1., rtol=0, atol=1e-14):
        raise ValueError('parent weights fail partition of unity')
    right = load@twist
    outer = np.r_[arrays['wall_faces'], arrays['bottom_faces']]
    samples, weights = triangle_samples(points[faces[outer]], 3)
    shape, _ = triangle_samples(np.eye(3)[None], 3)
    velocity = fused_cpu(samples.reshape(-1, 3), values, time)[:, 6:].reshape(samples.shape)
    flux = np.einsum('fqi,fi->fq', velocity, panels.normals[outer])
    for corner in range(3):
        np.add.at(right, faces[outer, corner], panels.areas[outer]
                  *np.einsum('q,q,fq->f', weights, shape[0, :, corner], flux))
    unknown = np.setdiff1d(np.arange(len(points)), fixed)
    potential = np.zeros(len(points))
    potential[fixed] = mapping@arrays['parent_psi']
    factor = cho_factor(energy[np.ix_(unknown, unknown)])
    potential[unknown] = cho_solve(factor, (right-energy@potential)[unknown])
    reaction = energy@potential-right
    result = dict(potential=potential, parent_map=mapping, body_load=load,
                  flux=reaction[fixed], right=right)
    if acceleration_columns:
        auxiliary = np.zeros((len(points), 6))
        auxiliary[unknown] = cho_solve(factor, load[unknown])
        result['acceleration_potentials'] = auxiliary
    result['relative_field_residual'] = float(np.linalg.norm(reaction[unknown])/max(np.linalg.norm(right[unknown]), 1e-30))
    return result


def run(snapshot, output):
    from .cuda_operator import prepare_energy, cp
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    snapshot = Path(snapshot)
    geometry = json.loads(snapshot.read_text())
    if not geometry.get('complete') or not geometry.get('geometry_roundoff_checks_passed'):
        raise ValueError('requires accepted complete boundary reference')
    arrays = dict(np.load(snapshot.with_suffix('.npz'), allow_pickle=False))
    waves = sea_states.loads_realization(json.dumps(geometry['actual_sea_realization']))
    values = coefficients(waves)
    points, faces = arrays['points'], arrays['faces']
    report = dict(complete=False, qualified=False,
        scope='Current nonflat diagnostic tank field and surface rates only; no pressure derivative, ocean coupling or trajectory',
        snapshot_sha256=sha256(snapshot.read_bytes()).hexdigest(),
        snapshot_arrays_sha256=sha256(snapshot.with_suffix('.npz').read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), timings=[])
    saved = {}
    start = perf_counter()
    try:
        for repetition in range(2):
            before = perf_counter()
            current = prepare_energy(faces, TrianglePanels(points[faces], 1))
            cp.cuda.get_current_stream().synchronize()
            report['timings'].append(dict(operation='energy', repetition=repetition, seconds=perf_counter()-before))
        twist = np.array([.03, -.02, .01, .01, -.015, .02])
        for repetition in range(2):
            before = perf_counter()
            field = current_field(arrays, current.energy, values, twist)
            report['timings'].append(dict(operation='field_and_six_responses', repetition=repetition, seconds=perf_counter()-before))
        before = perf_counter()
        single = current_field(arrays, current.energy, values, twist, acceleration_columns=False)
        report['timings'].append(dict(operation='field_without_six_responses', seconds=perf_counter()-before))
        np.testing.assert_allclose(single['potential'], field['potential'], rtol=1e-13, atol=1e-13)
        before = perf_counter()
        eta_rate, psi_rate, diagnostics = surface_rates(points, faces[arrays['free_faces']],
            arrays['free_nodes'], field['parent_map'], arrays['parent_psi'], field['flux'])
        report['timings'].append(dict(operation='surface_rates', seconds=perf_counter()-before))
        added_mass = WATER_DENSITY*field['body_load'].T@field['acceleration_potentials']
        asymmetry = np.linalg.norm(added_mass-added_mass.T)/np.linalg.norm(added_mass)
        if field['relative_field_residual'] > 1e-9 or asymmetry > 1e-10:
            raise ValueError('current field/response algebraic check failed')
        report.update(complete=True, seconds=perf_counter()-start, nodes=len(points), faces=len(faces),
            surface=diagnostics, relative_field_residual=field['relative_field_residual'],
            added_mass_reciprocity_relative=float(asymmetry), added_mass=added_mass.tolist())
        saved.update({k: v for k, v in field.items() if isinstance(v, np.ndarray)})
        saved.update(eta_rate=eta_rate, psi_rate=psi_rate, twist=twist)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))

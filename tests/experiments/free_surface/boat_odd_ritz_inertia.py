"""Saved Ritz six-port instantaneous added mass against independent solvers."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from topwater.constants import WATER_DENSITY, GRAVITY
from topwater.frames import rotation_matrix
from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from tests.capytaine_reference import instantaneous_inertia
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import RegularizedSnapshot, owned
from .odd_variational_energy import OddRitzEnergy


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def ratio(a, b):
    return float(a/b) if b else (0. if a == 0 else None)


def run(args, report):
    ritz, refinement, odd = owned(args.ritz), owned(args.refinement), owned(args.odd)
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    for name in ('water', 'phase', 'mode'):
        if ritz['parents'][name] != sha256(getattr(args, name).read_bytes()).hexdigest():
            raise ValueError('saved Ritz physical parents differ')
    if (refinement['parent_phase_sha256'] != sha256(args.phase.read_bytes()).hexdigest()
            or refinement['input_sha256'] != snapshot.mode['input_sha256']
            or refinement['panel_sha256'] != snapshot.mode['panel_sha256']):
        raise ValueError('archived body refinement is not the same phase/input/panels')
    # The Ritz report owns its body16 parent; that report in turn owns the odd
    # collocation matrix. Follow the explicit CLI parent instead of guessing.
    owner = owned(args.body16)
    if (sha256(args.body16.read_bytes()).hexdigest() != ritz['parents']['body16']
            or sha256(args.odd.read_bytes()).hexdigest() != owner['parents']['odd']):
        raise ValueError('collocation density comparison parent differs')
    relative = Path('tests/experiments/free_surface/odd_variational_energy.py')
    if sha256(relative.read_bytes()).hexdigest() != ritz['producer_sha256'].get(str(relative)):
        raise ValueError('saved Ritz energy owner source differs')
    with np.load(args.ritz.with_suffix('.npz')) as z:
        values = {name: z[name].copy() for name in
            ('HH', 'HR', 'RH', 'RR', 'BH', 'BR', 'source_areas', 'body_R', 'body_J',
             'body_weights', 'body_T', 'body_points')}
    with np.load(args.odd.with_suffix('.npz')) as z:
        jd = z['jd'].copy()
    panels = TrianglePanels(snapshot.triangles, -1)
    n, b = len(snapshot.patch.active_nodes), len(panels)
    np.testing.assert_array_equal(values['source_areas'], panels.areas)
    if values['body_R'].shape != values['body_J'].shape or values['body_R'].shape != (b*256, b):
        raise ValueError('saved body16 source trace dimensions differ')
    center = snapshot.inputs.initial_state.position
    owners = np.repeat(np.arange(b), 256)
    expected_T = np.column_stack((panels.normals[owners],
        np.cross(values['body_points']-center, panels.normals[owners])))
    np.testing.assert_array_equal(values['body_T'], expected_T)
    energy = OddRitzEnergy(*(values[name] for name in ('HH', 'HR', 'RH', 'RR', 'BH', 'BR', 'source_areas')))
    sigma = energy.solve_density(np.zeros((n, 6)), np.eye(6))
    potential = values['body_R']@sigma
    normal = values['body_J']@sigma
    weight, T = values['body_weights'], values['body_T']
    mass_ritz = WATER_DENSITY*energy.A
    quadrature_mass = WATER_DENSITY*T.T@(weight[:, None]*potential)
    body_error = normal-T
    surface = float(weight.sum())
    rms = np.sqrt((weight[:, None]*body_error**2).sum(axis=0)/surface)
    reference_rms = np.sqrt((weight[:, None]*T**2).sum(axis=0)/surface)
    mode = np.column_stack((panels.normals, np.cross(panels.centroids-center, panels.normals)))
    sigma_collocation = np.linalg.solve(jd, mode)
    sigma_norm = np.sqrt((panels.areas[:, None]*sigma**2).sum(axis=0))
    sigma_collocation_norm = np.sqrt((panels.areas[:, None]*sigma_collocation**2).sum(axis=0))
    difference_norm = np.sqrt((panels.areas[:, None]*(sigma-sigma_collocation)**2).sum(axis=0))
    rotation = rotation_matrix(snapshot.inputs.initial_state.orientation)
    dry = np.zeros((6, 6))
    dry[:3, :3] = snapshot.inputs.body.mass*np.eye(3)
    dry[3:, 3:] = rotation@snapshot.inputs.body.inertia_cg@rotation.T
    factor = cholesky(dry, lower=True)

    def whiten(matrix):
        left = solve_triangular(factor, matrix, lower=True)
        return solve_triangular(factor, left.T, lower=True).T

    def describe(matrix):
        white = whiten(matrix)
        return dict(diagonal=matrix.diagonal().tolist(), dry_mass_whitened=white.tolist(),
            dry_mass_whitened_raw_asymmetry_max=float(abs(white-white.T).max()))

    def compare(matrix):
        difference = mass_ritz-matrix
        white = whiten(difference)
        reference = whiten(matrix)
        return dict(difference=difference.tolist(), diagonal_relative_difference=[
            ratio(a, b) for a, b in zip(difference.diagonal(), matrix.diagonal(), strict=True)],
            dry_mass_whitened_difference=white.tolist(),
            dry_mass_whitened_difference_max=float(abs(white).max()),
            dry_mass_whitened_relative_operator_norm=ratio(np.linalg.norm(white, 2), np.linalg.norm(reference, 2)))

    arrays = dict(ritz_added_mass=mass_ritz, quadrature_added_mass=quadrature_mass,
        sigma_ritz=sigma, sigma_collocation=sigma_collocation, dry_mass=dry,
        body_potential=potential, body_neumann=normal, body_neumann_error=body_error,
        triangles=snapshot.triangles, world_cg=center)
    report.update(complete=False, qualified=False, density_kg_m3=WATER_DENSITY, gravity=GRAVITY,
        body_panels=b, world_cg=center.tolist(), physical_inputs_sha256=snapshot.mode['input_sha256'],
        panel_sha256=snapshot.mode['panel_sha256'], ritz=describe(mass_ritz),
        body_neumann_rms=rms.tolist(), body_rigid_normal_rms=reference_rms.tolist(),
        body_neumann_relative_rms=[ratio(a, b) for a, b in zip(rms, reference_rms, strict=True)],
        body_neumann_max=abs(body_error).max(axis=0).tolist(),
        area_weighted_sigma_norm=sigma_norm.tolist(), area_weighted_collocation_sigma_norm=sigma_collocation_norm.tolist(),
        area_weighted_sigma_relative_difference=[ratio(a, b) for a, b in zip(difference_norm, sigma_collocation_norm, strict=True)],
        instantaneous_pressure_port_bookkeeping_max=float(abs(quadrature_mass-mass_ritz).max()),
        archived_body_refinement=refinement['body_only_halfspace_refinement'],
        reference_refinement_scope='Owned same-phase 68/272/1088 direct/source diagonals; not a certified continuum interval',
        stages=[])
    checkpoint(args, report, arrays)
    started = perf_counter()
    direct = instantaneous_added_mass(panels, center, density=WATER_DENSITY)
    arrays['direct_added_mass'] = direct.added_mass
    report.update(direct=describe(direct.added_mass), ritz_minus_direct=compare(direct.added_mass),
                  direct_diagnostics=direct.diagnostics)
    report['stages'].append(dict(name='production_instantaneous_added_mass', seconds=perf_counter()-started))
    checkpoint(args, report, arrays)
    started = perf_counter()
    capy, damping = instantaneous_inertia(snapshot.triangles, center, snapshot.inputs.body.mass,
                                         density=WATER_DENSITY, gravity=GRAVITY)
    if not np.isfinite(capy).all() or not np.isfinite(damping).all():
        raise FloatingPointError('nonfinite independent Capytaine reference')
    arrays.update(capytaine_added_mass=capy, capytaine_damping=damping)
    report.update(capytaine=describe(capy), ritz_minus_capytaine=compare(capy),
        capytaine_damping_max=float(abs(damping).max()), process_peak_rss=peak_rss_bytes(), complete=True,
        limits=['psi=0 gives zero acceleration potential on the WHOLE plane; no wave-patch far-port ambiguity',
            'Auxiliary unit-acceleration pressure reaction is minus added_mass times acceleration',
            'Cached physical body quadrature validates the potential port, not an additional pressure load',
            'All matrices/references retain their raw asymmetry; no eigenvalue filtering or force repair',
            'No new refinement, evolved free surface, moving-domain, nonlinear6DOF, 60s or runtime qualification'])
    report['stages'].append(dict(name='capytaine_six_radiation_columns_infinity', seconds=perf_counter()-started))
    checkpoint(args, report, arrays)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'ritz', 'body16', 'odd', 'refinement', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output and source archive required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True)
    hashes = {}
    for source in [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'), *Path(__file__).parent.glob('*.py')]:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    report = dict(complete=False, producer_sha256=hashes, parents={name: sha256(getattr(args, name).read_bytes()).hexdigest()
        for name in ('water', 'phase', 'mode', 'ritz', 'body16', 'odd', 'refinement')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))

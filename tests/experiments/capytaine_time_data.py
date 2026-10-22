"""Export raw finite-depth Capytaine data for a validation-only time reference.

Requires an explicit positive frequency grid; it supplies no zero-frequency
extrapolation, interpolation, high-frequency tail, kernel or convergence claim.
Run with the correlation dependency group (Capytaine 3.0.0). --reference is the
pinned 704-panel reference JSON; --mesh is its actual saved wet-triangle NPZ
(these need not have the same stem). --refine-once produces exact coplanar
subdivision children, not a claim of identity with independently clipped meshes.
Completed six-column results survive a later frequency/A-infinity failure.
New exports default to isolated Prony seed 0, reset per frequency/A-infinity
call; --prony-seed chooses another reproducible realization without changing
Capytaine's approximation backend or tolerance.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from numbers import Integral
import subprocess

import numpy as np

from tests.capytaine_reference import MODES, frequency_response
from tests.experiments.capytaine_mesh_cost import subdivide
from tests.experiments.free_surface.canoe_wet_reference import pinned_reference


def _file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _array_hash(values):
    return sha256(np.asarray(values, dtype='<f8').tobytes()).hexdigest()


def _frequencies(values):
    values = np.asarray(values, dtype=float)
    if (values.ndim != 1 or not len(values) or not np.all(np.isfinite(values))
            or np.any(values <= 0) or np.any(np.diff(values) <= 0)):
        raise ValueError('supply a nonempty strictly increasing positive omega grid in rad/s')
    return values


def _matrix(values):
    values = np.asarray(values, dtype=float)
    if values.shape != (6, 6) or not np.all(np.isfinite(values)):
        raise ValueError('Capytaine must return all 36 finite real coefficients')
    return values


def run(reference, mesh, frequencies, output, *, refine_once=False, grid_file=None, prony_seed=0):
    """One pinned physical case, checkpointed after each complete helper call."""
    reference, mesh, output = Path(reference), Path(mesh), Path(output).with_suffix('.json')
    archive = output.with_suffix('.npz')
    if output.exists() or archive.exists():
        raise FileExistsError('choose a fresh output stem; existing results are preserved')
    if isinstance(prony_seed, bool) or not isinstance(prony_seed, Integral) or prony_seed < 0:
        raise ValueError('new datasets require an explicit nonnegative integer Prony seed')
    prony_seed = int(prony_seed)
    frequencies = _frequencies(frequencies)
    output.parent.mkdir(parents=True, exist_ok=True)
    domain, pinned, _, provenance = pinned_reference(reference)
    with np.load(mesh, allow_pickle=False) as saved:
        coarse = np.asarray(saved['wetted_triangles'], dtype=float)
    if not np.array_equal(coarse, pinned):
        raise ValueError('saved wet triangles differ from the exact canonical 704-panel pin')
    triangles = subdivide(coarse) if refine_once else coarse.copy()
    center, mass = domain.rest.position.copy(), float(domain.body.mass)
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), root/'capytaine_reference.py',
        root/'experiments/capytaine_mesh_cost.py',
        root/'experiments/free_surface/canoe_wet_reference.py',
        root/'experiments/free_surface/p4_radiation.py',
        root/'experiments/free_surface/moving_canoe.py',
        root/'experiments/free_surface/canoe_contact.py', root/'reference_hull.py']
    source_hashes = {str(path.relative_to(root.parent)): _file_hash(path) for path in sources}
    input_hashes = dict(reference_json=_file_hash(reference), saved_mesh_npz=_file_hash(mesh))
    if grid_file is not None:
        input_hashes['frequency_grid_file'] = _file_hash(grid_file)
    import capytaine as cpt
    if cpt.__version__ != '3.0.0':
        raise ValueError('this reference producer requires Capytaine 3.0.0')
    report = dict(qualified=False, complete=False,
        scope='Raw Capytaine data for a time-domain validation reference; not a Topwater solver.',
        source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'],
            cwd=root.parent, text=True).strip(), source_sha256=source_hashes,
        input_sha256=input_hashes, input_paths=dict(reference=str(reference), mesh=str(mesh),
            frequency_grid_file=None if grid_file is None else str(grid_file)),
        physical_provenance=provenance, capytaine_version=cpt.__version__,
        mode_order=list(MODES), matrix_axes='row=response force/moment; column=imposed rigid mode',
        rotation_center=center.tolist(), body_mass=mass, water_depth=4., density=1000., gravity=9.81,
        forward_speed=0., coarse_wetted_panels=len(coarse), coarse_triangles_sha256=_array_hash(coarse),
        wetted_panels=len(triangles), wetted_triangles_sha256=_array_hash(triangles),
        refine_once=bool(refine_once),
        mesh_provenance='Exact saved canonical pin; optional one coplanar four-child subdivision. '
            'No equality with a separately clipped 2816-panel hash is assumed.',
        requested_omega=frequencies.tolist(), completed_omega=[], prony_seed=prony_seed,
        prony_seed_policy='Default seed 0; reset independently for each completed frequency/A-infinity call. '
            'No backend/tolerance change; historical difference attribution remains unproven.',
        prony_realizations=[], infinity_prony_realization=None,
        frequency_grid_sha256=_array_hash(frequencies),
        convention='Re(amplitude*exp(-i*omega*t)); radiation force=(omega^2*A+i*omega*B)*displacement; '
            'time force=-A_infinity*acceleration-integral(K(t-s)*velocity(s),s=0..t).',
        units=dict(omega='rad/s', translation='m', rotation='rad', force='N', moment='N m',
            added_mass='row force unit / column acceleration unit',
            radiation_damping='row force unit / column velocity unit'),
        added_mass_infinite_depth=4.,
        infinity_provenance='Separate omega=inf call to frequency_response at water_depth=4; '
            'the infinite-depth instantaneous_inertia helper is never used. Capytaine 3.0 '
            'automatically selects high_freq Green singularities at infinity; its exported '
            'settings retain the configured singularity value.',
        infinity_effective_green_singularities='high_freq',
        limitations='No zero frequency, extrapolation, interpolation, tail closure, fitting to Topwater, '
            'IRF or force-history convergence. Frequency, time and spatial refinement remain required.')
    added, damping, infinite, infinite_damping, settings = [], [], None, None, None

    def checkpoint():
        # Replace only this run's own checkpoint. JSON binds the completed NPZ
        # bytes; an interrupted two-file replacement is detectable by its hash.
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = archive.with_suffix('.npz.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream, omega=np.asarray(report['completed_omega']),
                requested_omega=frequencies, added_mass=np.asarray(added).reshape(-1, 6, 6),
                radiation_damping=np.asarray(damping).reshape(-1, 6, 6),
                added_mass_infinite=np.empty((0, 6, 6)) if infinite is None else infinite[None],
                radiation_damping_infinite=np.empty((0, 6, 6)) if infinite_damping is None else infinite_damping[None],
                wetted_triangles=triangles, coarse_wetted_triangles=coarse, rotation_center=center,
                mode_order=np.asarray(MODES), water_depth=4., density=1000., gravity=9.81,
                source_sha256_json=json.dumps(source_hashes, sort_keys=True),
                input_sha256_json=json.dumps(input_hashes, sort_keys=True),
                prony_seed_decimal=str(prony_seed), prony_realizations_json=json.dumps(report['prony_realizations'], sort_keys=True),
                infinity_prony_realization_json=json.dumps(report['infinity_prony_realization'], sort_keys=True),
                capytaine_version=cpt.__version__, green_function_settings_json=json.dumps(settings, sort_keys=True))
        temporary.replace(archive)
        report.update(added_mass_infinite_complete=infinite is not None,
            capytaine_green_function_settings=settings, arrays_sha256=_file_hash(archive))
        temporary = output.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        temporary.replace(output)

    def solve(omega):
        nonlocal settings
        result = frequency_response(triangles, center, mass, [omega], [], water_depth=4., prony_seed=prony_seed)
        current = result['green_function_settings']
        if current != provenance['capytaine_green_function_settings']:
            raise ValueError('runtime Green settings differ from the pinned Capytaine reference')
        if settings is not None and settings != current:
            raise ValueError('Green settings changed between completed calls')
        a = _matrix(result['added_mass'][0])
        b = _matrix(result['radiation_damping'][0])
        settings = current
        return a, b, result['prony_realization']

    try:
        for omega in frequencies:
            a, b, realization = solve(float(omega))
            added.append(a)
            damping.append(b)
            report['completed_omega'].append(float(omega))
            report['prony_realizations'].append(dict(omega=float(omega), realization=realization))
            checkpoint()
            print(json.dumps(dict(stage='completed frequency checkpoint', omega=float(omega))), flush=True)
        infinite, infinite_damping, realization = solve(np.inf)
        report['infinity_prony_realization'] = realization
        report['complete'] = True
        checkpoint()
        print(json.dumps(dict(stage='completed finite-depth A-infinity checkpoint')), flush=True)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        try:
            checkpoint()
        except BaseException as save_error:
            error.add_note('Could not update failure checkpoint: '+repr(save_error))
        raise
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    grid = parser.add_mutually_exclusive_group(required=True)
    grid.add_argument('--omega', type=float, nargs='+', help='explicit increasing positive rad/s values')
    grid.add_argument('--omega-file', type=Path, help='plain text one-column positive rad/s grid')
    parser.add_argument('--refine-once', action='store_true')
    parser.add_argument('--prony-seed', type=int, default=0,
                        help='nonnegative isolated Capytaine Prony RNG seed (default: 0)')
    parser.add_argument('--output', type=Path, required=True, help='fresh JSON/NPZ output stem')
    args = parser.parse_args()
    omega = np.loadtxt(args.omega_file, ndmin=1) if args.omega_file else args.omega
    run(args.reference, args.mesh, omega, args.output,
        refine_once=args.refine_once, grid_file=args.omega_file, prony_seed=args.prony_seed)

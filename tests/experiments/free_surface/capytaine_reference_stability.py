"""Offline704-to2816-panel Capytaine radiation-reference stability check.

Only existing reports and saved wet triangles are read. No Capytaine or fluid
solve is performed. Two reference resolutions cannot establish
reference convergence, Topwater accuracy, or nonlinear qualification.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from tests.experiments.capytaine_mesh_cost import subdivide
from .validate_radiation import coefficient_gates, MODES


FREQUENCIES = [2., 4., 6.]
LENGTH = 3.5
DENSITY = 1000.
GRAVITY = 9.81


def _triangles(values, report, count):
    values = np.asarray(values, dtype=float)
    if values.shape != (count, 3, 3) or not np.all(np.isfinite(values)):
        raise ValueError(f'{count} finite saved wet triangles required')
    digest = sha256(np.asarray(values, dtype='<f8').tobytes()).hexdigest()
    if report['wetted_panels'] != count or report['wetted_triangles_sha256'] != digest:
        raise ValueError('saved wet triangle count/hash differs from origin report')
    return values


def compare(coarse, fine, origin, coarse_triangles, fine_triangles):
    """Keep all raw matrix entries; fine family defines the declared floor."""
    for report in (coarse, origin):
        if (report['mode_order'] != list(MODES) or report['water_depth'] != 4.
                or [case['omega'] for case in report['cases']] != FREQUENCIES):
            raise ValueError('origin mode/depth/frequency metadata differs from fixed family')
    if (fine['frequencies'] != FREQUENCIES or fine['water_depth'] != 4.
            or coarse['rotation_center'] != origin['rotation_center']
            or np.shape(coarse['rotation_center']) != (3,)
            or not np.all(np.isfinite(coarse['rotation_center']))):
        raise ValueError('fine reference or origin CG/depth/frequency metadata differs')
    if coarse['wet_subdivisions'] != 2 or origin['wet_subdivisions'] != 3:
        raise ValueError('origin numerical subdivisions differ from the704-to2816 reference pair')
    for key in ('rotation_center', 'mode_order'):
        if key in fine and fine[key] != coarse[key]:
            raise ValueError(f'explicit fine {key} differs from origin/coarse metadata')
    for report in (coarse, fine, origin):
        for name, expected in (('density', DENSITY), ('rho', DENSITY), ('gravity', GRAVITY), ('g', GRAVITY)):
            if name in report and report[name] != expected:
                raise ValueError('explicit physical constants differ from the reference normalization')
    if fine['green_function_settings'] != coarse['capytaine_green_function_settings']:
        raise ValueError('coarse/fine Green-function settings differ')
    if ('capytaine_green_function_settings' in origin
            and origin['capytaine_green_function_settings'] != fine['green_function_settings']):
        raise ValueError('fine producer-origin Green-function settings differ')
    a = _triangles(coarse_triangles, coarse, 704)
    b = _triangles(fine_triangles, origin, 2816)
    if fine['panels'] != len(b) or fine['wetted_triangles_sha256'] != origin['wetted_triangles_sha256']:
        raise ValueError('fine Capytaine triangle hash/count differs from its producer origin')
    expected = subdivide(a)
    difference = b-expected
    error = float(np.max(abs(difference)))
    scale = max(1., LENGTH, float(np.max(abs(a))), float(np.max(abs(b))))
    bound = 256*np.finfo(float).eps*scale
    normal = np.cross(b[:, 1]-b[:, 0], b[:, 2]-b[:, 0])
    expected_normal = np.cross(expected[:, 1]-expected[:, 0], expected[:, 2]-expected[:, 0])
    if error > bound or np.any(np.einsum('fi,fi->f', normal, expected_normal) <= 0):
        raise ValueError('fine wet triangles are not ordered coplanar subdivision children within roundoff')
    fine_a, fine_b = np.asarray(fine['added_mass']), np.asarray(fine['radiation_damping'])
    if fine_a.shape != (3, 6, 6) or fine_b.shape != (3, 6, 6) or not np.all(np.isfinite(fine_a)) or not np.all(np.isfinite(fine_b)):
        raise ValueError('finite complete fine6x6 reference family required')
    d = np.array([1., 1., 1., LENGTH, LENGTH, LENGTH])
    scale = DENSITY*LENGTH**3*d[:, None]*d[None, :]
    cases = []
    for i, omega in enumerate(FREQUENCIES):
        ca = np.asarray(coarse['cases'][i]['capytaine']['added_mass'])
        cb = np.asarray(coarse['cases'][i]['capytaine']['radiation_damping'])
        if ca.shape != (6, 6) or cb.shape != (6, 6) or not np.all(np.isfinite(ca)) or not np.all(np.isfinite(cb)):
            raise ValueError('finite complete coarse6x6 reference family required')
        cases.append(dict(omega=omega, coarse_added_mass=ca.tolist(), fine_added_mass=fine_a[i].tolist(),
            coarse_radiation_damping=cb.tolist(), fine_radiation_damping=fine_b[i].tolist(),
            added_mass_difference_coarse_minus_fine=(ca-fine_a[i]).tolist(),
            radiation_damping_difference_coarse_minus_fine=(cb-fine_b[i]).tolist(),
            added_mass_nondimensional=(ca/scale).tolist(),
            reference_added_mass_nondimensional=(fine_a[i]/scale).tolist(),
            damping_nondimensional=(cb/(omega*scale)).tolist(),
            reference_damping_nondimensional=(fine_b[i]/(omega*scale)).tolist(),
            added_mass_difference_nondimensional=((ca-fine_a[i])/scale).tolist(),
            damping_difference_nondimensional=((cb-fine_b[i])/(omega*scale)).tolist()))
    gates = coefficient_gates(cases)
    return dict(qualified=False, reference_stability_gate_passed=all(value['passed'] for value in gates.values()),
        reference_resolution_levels=2, length=LENGTH, density=DENSITY, gravity=GRAVITY,
        mode_order=list(MODES), rotation_center=coarse['rotation_center'], water_depth=4.,
        green_function_settings=fine['green_function_settings'],
        geometry=dict(coarse_panels=704, fine_panels=2816,
            coarse_sha256=coarse['wetted_triangles_sha256'], fine_sha256=fine['wetted_triangles_sha256'],
            ordered_subdivision_exact=bool(np.array_equal(b, expected)),
            maximum_coordinate_roundoff_difference=error,
            coordinate_difference_rms=float(np.sqrt(np.mean(difference**2))),
            geometry_roundoff_bound=bound, winding_preserved=True,
            correspondence='Same ordered four-child subdivision of every coarse triangle'),
        cases=cases, gates=gates,
        legacy_provenance='Fine Capytaine JSON omits CG, mode order, rho and g. CG/mode order '
            'are linked through its triangle-hash-matched origin report and reviewed producer/helper. '
            'rho=1000 and g=9.81 follow the helper source; these legacy files do not certify '
            'the producer runtime source revision or independently serialize all runtime constants.',
        scope='Reference-solver stability only, using all36 entries at all3frequencies and the '
            'existing family gate with fine as reference. Passing two resolutions does not establish '
            'asymptotic reference convergence, Topwater coefficient agreement or nonlinear qualification.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse', type=Path, required=True)
    parser.add_argument('--fine', type=Path, required=True)
    parser.add_argument('--fine-origin', type=Path, required=True)
    parser.add_argument('--coarse-mesh', type=Path)
    parser.add_argument('--fine-mesh', type=Path)
    parser.add_argument('--producer', type=Path, default=Path('artifacts/capytaine_fine_frequency_family.py'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('choose a new stability report; existing reports are preserved')
    coarse_mesh = args.coarse_mesh or args.coarse.with_suffix('.npz')
    fine_mesh = args.fine_mesh or args.fine_origin.with_suffix('.npz')
    with np.load(coarse_mesh) as saved:
        coarse_triangles = saved['wetted_triangles']
    with np.load(fine_mesh) as saved:
        fine_triangles = saved['wetted_triangles']
    report = compare(json.loads(args.coarse.read_text()), json.loads(args.fine.read_text()),
                     json.loads(args.fine_origin.read_text()), coarse_triangles, fine_triangles)
    sources = dict(coarse=args.coarse, fine=args.fine, fine_origin=args.fine_origin,
        coarse_mesh=coarse_mesh, fine_mesh=fine_mesh, producer=args.producer,
        reference_helper=Path(__file__).parents[2]/'capytaine_reference.py',
        verifier=Path(__file__), subdivision=Path(__file__).parents[1]/'capytaine_mesh_cost.py',
        family_gate=Path(__file__).with_name('validate_radiation.py'))
    report['input_provenance'] = {name: dict(path=str(path.resolve()), sha256=sha256(path.read_bytes()).hexdigest())
                                  for name, path in sources.items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(qualified=False, reference_stability_gate_passed=report['reference_stability_gate_passed'],
        geometry=report['geometry'], maximum_tolerance_ratios={name: value['maximum_tolerance_ratio']
        for name, value in report['gates'].items()})))


if __name__ == '__main__':
    main()

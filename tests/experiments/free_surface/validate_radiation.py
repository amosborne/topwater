"""Offline full-matrix coefficient gates for six-mode radiation reports.

Passing these coefficient gates never establishes convergence or qualification.
No fluid or Capytaine solves are performed here.
"""
import argparse
import json
from pathlib import Path

import numpy as np

MODES = ('surge', 'sway', 'heave', 'roll', 'pitch', 'yaw')
SETTINGS = ('n', 'gpu', 'wet_subdivisions', 'contact_mesh', 'exterior_scale', 'periods',
            'ramp_periods', 'fit_periods', 'amplitude', 'rtol', 'derivative_step',
            'max_step', 'water_depth', 'sponge_strength', 'sponge_buffer',
            'domain_half_widths', 'sponge_widths', 'maximum_surface_edge',
            'formulation', 'integrated_fluid_variables')
REFERENCE_SETTINGS = ('omega', 'n', 'wet_subdivisions', 'contact_mesh', 'water_depth')


def _array(value, shape, label):
    array = np.asarray(value, dtype=float)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f'{label} must be finite with shape {shape}')
    return array


def _same_settings(report, expected, keys):
    for key in keys:
        if key not in report or key not in expected or report[key] != expected[key]:
            raise ValueError(f'inconsistent or missing setting: {key}')


def coefficient_gates(cases):
    """Apply the declared family tolerance to complete nondimensional matrices.

    This gate is shared by harmonic transients and exact rest linearization;
    it makes no convergence or nonlinear-qualification assertion.
    """
    if not cases:
        raise ValueError('nonempty complete coefficient family required')
    gates={}
    for family,actual_key,reference_key in (
            ('added_mass','added_mass_nondimensional','reference_added_mass_nondimensional'),
            ('damping','damping_nondimensional','reference_damping_nondimensional')):
        actual=np.array([_array(case[actual_key],(6,6),actual_key) for case in cases])
        reference=np.array([_array(case[reference_key],(6,6),reference_key) for case in cases])
        # One floor over every supplied frequency and all 36 matrix entries.
        floor=max(1e-10,.005*float(np.max(np.abs(reference))))
        tolerance=.05*np.abs(reference)+floor
        ratio=np.abs(actual-reference)/tolerance
        worst=np.unravel_index(np.argmax(ratio),ratio.shape)
        gates[family]=dict(passed=bool(np.all(ratio<=1)),family_floor=floor,
            maximum_tolerance_ratio=float(ratio[worst]),
            worst_omega=cases[worst[0]]['omega'],worst_force_mode=MODES[worst[1]],
            worst_radiating_mode=MODES[worst[2]],tolerance_ratios=ratio.tolist())
    return gates


def aggregate(reports, *, references=(), length=3.5, density=1000.):
    """Compare complete matrices over exactly the supplied frequency set.

    Optional reference-only reports may be short Capytaine smoke runs. Every
    frequency must have one reference, either embedded or explicitly supplied.
    Reference reports need matching hull resolution/depth/frequency, whereas
    all six evolved reports must also have matching run controls and domain.
    """
    if not np.isfinite(length) or length <= 0 or not np.isfinite(density) or density <= 0:
        raise ValueError('positive finite normalization length and density required')
    groups = {}
    for report in reports:
        omega = report.get('omega')
        if not isinstance(omega, (int, float)) or not np.isfinite(omega) or omega <= 0:
            raise ValueError('positive finite omega required')
        mode = report.get('mode')
        if mode not in MODES or report.get('mode_order') != list(MODES):
            raise ValueError('invalid mode or mode ordering')
        expected_unit = 'm' if MODES.index(mode) < 3 else 'rad'
        if report.get('amplitude_unit') != expected_unit:
            raise ValueError('invalid forcing amplitude unit')
        group = groups.setdefault(omega, {})
        if mode in group:
            raise ValueError(f'duplicate mode {mode} at omega={omega}')
        group[mode] = report
    if not groups:
        raise ValueError('no radiation reports supplied')
    reference_reports = list(references)
    if any(report.get('omega') not in groups for report in reference_reports):
        raise ValueError('reference frequency is outside the supplied radiation set')
    cases = []
    for omega, group in sorted(groups.items()):
        missing = set(MODES)-set(group)
        if missing:
            raise ValueError(f'missing modes at omega={omega}: {sorted(missing)}')
        base = group['surge']
        actual_a, actual_b = np.empty((6, 6)), np.empty((6, 6))
        ref_a = ref_b = None
        for column, mode in enumerate(MODES):
            report = group[mode]
            _same_settings(report, base, SETTINGS)
            if not np.isfinite(report['amplitude']) or report['amplitude'] <= 0:
                raise ValueError('positive finite forcing amplitude required')
            if report['periods'] < report['ramp_periods']+report['fit_periods']:
                raise ValueError('insufficient post-ramp periods for the requested fit')
            fits = report.get('column_harmonic_fits')
            if not isinstance(fits, list) or len(fits) != 6:
                raise ValueError(f'missing complete fitted column for {mode}')
            actual_a[:, column] = _array([row.get('added_mass') for row in fits], (6,), 'fitted added mass')
            actual_b[:, column] = _array([row.get('radiation_damping') for row in fits], (6,), 'fitted damping')
        candidates = list(group.values())+[report for report in reference_reports if report['omega'] == omega]
        for report in candidates:
            reference = report.get('capytaine')
            if reference is None:
                continue
            _same_settings(report, base, REFERENCE_SETTINGS)
            a = _array(reference.get('added_mass'), (6, 6), 'reference added mass')
            b = _array(reference.get('radiation_damping'), (6, 6), 'reference damping')
            if ref_a is not None and (not np.array_equal(a, ref_a) or not np.array_equal(b, ref_b)):
                raise ValueError(f'inconsistent reference matrices at omega={omega}')
            ref_a, ref_b = a, b
        if ref_a is None:
            raise ValueError(f'missing Capytaine reference at omega={omega}')
        d = np.array([1., 1., 1., length, length, length])
        scale = density*length**3*d[:, None]*d[None, :]
        cases.append(dict(omega=omega, settings={key: base[key] for key in SETTINGS},
            added_mass=actual_a.tolist(), radiation_damping=actual_b.tolist(),
            reference_added_mass=ref_a.tolist(), reference_radiation_damping=ref_b.tolist(),
            added_mass_nondimensional=(actual_a/scale).tolist(),
            reference_added_mass_nondimensional=(ref_a/scale).tolist(),
            damping_nondimensional=(actual_b/(omega*scale)).tolist(),
            reference_damping_nondimensional=(ref_b/(omega*scale)).tolist()))
    gates=coefficient_gates(cases)
    return dict(qualified=False, coefficient_gates_passed=all(gate['passed'] for gate in gates.values()),
        length=length, density=density, mode_order=list(MODES), cases=cases, gates=gates,
        unresolved=['amplitude convergence', 'time-step convergence', 'spatial convergence',
                    'domain and sponge convergence', 'startup/transient convergence',
                    'reference-solver convergence'],
        scope='Only the supplied complete frequency set is tested. Coefficient agreement alone '
              'does not establish nonlinear wave/body qualification. Input reports lack an exact '
              'geometry fingerprint; geometry provenance still requires independent review.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reports', nargs='+', type=Path)
    parser.add_argument('--reference', action='append', type=Path, default=[])
    parser.add_argument('--length', type=float, default=3.5)
    parser.add_argument('--density', type=float, default=1000.)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        report = aggregate([json.loads(path.read_text()) for path in args.reports],
            references=[json.loads(path.read_text()) for path in args.reference],
            length=args.length, density=args.density)
    except (ValueError, KeyError, TypeError) as error:
        parser.error(str(error))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(qualified=False, coefficient_gates_passed=report['coefficient_gates_passed'],
        frequencies=[case['omega'] for case in report['cases']],
        maximum_tolerance_ratios={name: gate['maximum_tolerance_ratio']
                                 for name, gate in report['gates'].items()},
        unresolved=report['unresolved'])))
    return 0 if report['coefficient_gates_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

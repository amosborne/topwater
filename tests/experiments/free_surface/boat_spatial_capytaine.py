"""Independent same-input/same-panel radiation reference for the spatial probe."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil

import numpy as np

from topwater import SimulationInputs
from topwater.constants import WATER_DENSITY, GRAVITY
from topwater.hydrostatics import hydrostatics
from tests.capytaine_reference import radiation_columns
from .boat_spatial_work import subdivide


def encode(value):
    if isinstance(value, np.ndarray):
        return dict(real=value.real.tolist(), imag=value.imag.tolist()) if np.iscomplexobj(value) else value.tolist()
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    return value


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--work-report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    inputs = SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8'))
    work = json.loads(args.work_report.read_text(encoding='utf-8'))
    if sha256(inputs.to_json().encode()).hexdigest() != work['input_sha256']:
        parser.error('work report and physical inputs differ')
    if args.output.exists():
        parser.error('output already exists')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    hashes = {}
    for source in [*Path('src/topwater').glob('*.py'), Path('tests/capytaine_reference.py'),
                   Path('tests/laplace_reference.py'), Path(__file__),
                   Path(__file__).with_name('boat_spatial_work.py'),
                   Path(__file__).with_name('panel_gradient_reference.py')]:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    wet = hydrostatics(inputs.body, inputs.initial_state, density=WATER_DENSITY, gravity=GRAVITY)
    rows = []
    for refinement in sorted({row['refinement'] for row in work['rows']}):
        triangles = wet.wetted_triangles
        for _ in range(refinement):
            triangles = subdivide(triangles)
        panel_hash = sha256(triangles.tobytes()).hexdigest()
        matching = [row for row in work['rows'] if row['refinement'] == refinement]
        if any(row['panel_sha256'] != panel_hash for row in matching):
            raise ValueError('actual physical panels differ from work report')
        for omega in [np.inf, 2*np.pi*work['frequency_hz']]:
            result = radiation_columns(triangles, inputs.initial_state.position, inputs.body.mass,
                omega, ('Heave', 'Roll'), density=WATER_DENSITY, gravity=GRAVITY)
            comparisons = {}
            if np.isinf(omega):
                reference = result['added_mass'][[2, 3]]
                for row in matching:
                    val = np.array(row['added_mass_heave_roll'])
                    comparisons[row['representation']] = dict(
                        matrix=val.tolist(), absolute_error=(val-reference).tolist(),
                        diagonal_relative_error=(np.diag(val-reference)/np.diag(reference)).tolist())
            rows.append(dict(refinement=refinement, panel_count=len(triangles), panel_sha256=panel_hash,
                reference=encode(result), comparisons=comparisons))
            print(json.dumps(rows[-1]), flush=True)
    args.output.write_text(json.dumps(dict(complete=True, qualified=False,
        scope='Independent same-panel radiation coefficients; no time-domain correlation yet',
        work_report_sha256=sha256(args.work_report.read_bytes()).hexdigest(),
        input_sha256=work['input_sha256'], producer_sha256=hashes, rows=rows), indent=2, allow_nan=False), encoding='utf-8')

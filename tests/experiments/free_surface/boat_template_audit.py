"""Two cold exports and frozen-snapshot guarded-template falsifiers; no PDE."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
import numpy as np

from topwater import InitialState
from tests.reference_boat import reference_body
from .boat_world_boundary import construct
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .p4_fluid_ordering import fresh


def run(input_stem, output):
    stem, output = Path(input_stem), fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(complete=False, qualified=False, scope='Guarded current-geometry proof only; no field or trajectory',
                  cases=[], comparisons=[], input_sha256={}, source_sha256={
                      name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                      for name in ('boat_template_audit.py', 'boat_world_boundary.py',
                                   'boat_parent_clip_proof.py', 'boat_geometry_template.py')})
    start = perf_counter()
    saved = {}

    def inputs(name):
        path = Path(str(stem)+'-'+name+'.npz')
        report['input_sha256'][path.name] = sha256(path.read_bytes()).hexdigest()
        return dict(np.load(path, allow_pickle=False))

    body = reference_body()

    def cold(name, arrays):
        before = perf_counter()
        _, result, metadata = construct(output.with_name(output.stem+'-'+name+'.json'), body=body,
            pose=InitialState(arrays['pose_position'], arrays['pose_orientation']),
            parent_xy=arrays['parent_nodes'][:, :2], parents=arrays['parents'],
            parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'], export_template=True)
        report['cases'].append(dict(name=name, cold_seconds=perf_counter()-before,
                                   template_export=metadata['template_export']))
        return result

    def check(name, template, target, expected_valid):
        before = perf_counter()
        try:
            current, diagnostics = evaluate_template(template,
                pose=InitialState(target['pose_position'], target['pose_orientation']),
                parent_eta=target['parent_nodes'][:, 2], parent_psi=target['parent_psi'])
        except InvalidTemplate as error:
            report['comparisons'].append(dict(name=name, valid=False, expected_valid=expected_valid,
                                             seconds=perf_counter()-before, reason=str(error)))
            if expected_valid:
                raise
            return
        if not expected_valid:
            raise AssertionError(f'{name}: changed arrangement was accepted')
        if not np.array_equal(current['faces'], target['faces']):
            raise AssertionError('cold reference mesh topology differs from frozen target')
        errors = dict(points_m=float(np.max(abs(current['points']-target['points']))),
            material_m=float(np.max(abs(current['body_material_vertices']-target['body_material_vertices']))))
        # Compare full fixed-parent maps independently of sparse row ordering.
        mappings = []
        for arrays in (current, target):
            mapping = np.zeros((len(arrays['points']), len(arrays['parent_nodes'])))
            for k in range(3):
                ids = arrays['free_parent_indices'][:, k]; valid = ids >= 0
                mapping[arrays['free_nodes'][valid], ids[valid]] = arrays['free_parent_weights'][valid, k]
            mappings.append(mapping)
        errors['parent_map'] = float(np.max(abs(mappings[0]-mappings[1])))
        if max(errors.values()) > 1e-9:
            raise AssertionError(f'geometry proof 1e-9 absolute coordinate/map budget failed: {errors}')
        report['comparisons'].append(dict(name=name, valid=True, expected_valid=True,
            seconds=perf_counter()-before, errors=errors, diagnostics=diagnostics))
        saved[name+'_points'] = current['points']

    try:
        plus, half, minus, base = (inputs(name) for name in
            ('level0-sign1', 'level1-sign1', 'level0-sign-1', 'base'))
        template = cold('plus-template', plus)
        check('self', template, plus, True)
        check('same-side-half', template, half, True)
        check('cross-side', template, minus, False)
        base_template = cold('base-template', base)
        check('base-to-plus', base_template, plus, False)
        report.update(complete=True, seconds=perf_counter()-start)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-stem', required=True)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))

"""Matched exact-predicate cache timings on actual failed-stage boat states."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import InitialState
from .boat_world_boundary import construct
from .p4_fluid_ordering import fresh


def run(output, repeats=2):
    if repeats not in (1, 2, 3):
        raise ValueError('bounded benchmark supports one to three pairs')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fixture = Path(__file__).resolve().parents[2]/'reference_boat.py'
    fixture_hash = sha256(fixture.read_bytes()).hexdigest()
    report = dict(complete=False, qualified=False,
        scope='Exact geometry work only on prescribed failed-stage canonical-boat states; no fluid or simulation throughput claim.',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_geometry_cost.py', 'boat_world_boundary.py',
                         'boat_parent_clip_proof.py', 'boat_geometry_template.py')},
        fixture_sha256=fixture_hash, states=[], rows=[])
    try:
        for name in ('boat-state-275ms.json', 'boat-state-337ms.json'):
            path = fixture.with_name('data')/name
            case = json.loads(path.read_text())
            if case['fixture_sha256'] != fixture_hash:
                raise ValueError('state mass/geometry fixture changed')
            report['states'].append(dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest()))
            state = np.array(case['state'])
            pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
            reference = None
            for repeat in range(repeats):
                for enabled in ((False, True) if repeat % 2 == 0 else (True, False)):
                    target = output.with_name(output.stem+f'-{name[:-5]}-r{repeat}-cache{int(enabled)}.json')
                    before = perf_counter()
                    _, arrays, geometry = construct(target, pose=pose,
                        parent_eta=state[13:28], parent_psi=state[28:],
                        export_template=True, cache_signs=enabled)
                    seconds = perf_counter()-before
                    if reference is None:
                        reference = arrays
                    if set(reference) != set(arrays) or any(
                            not np.array_equal(reference[k], arrays[k]) for k in reference):
                        raise ValueError('exact predicate cache changed geometry/ownership arrays')
                    report['rows'].append(dict(state=name, repeat=repeat, cache_signs=enabled,
                        seconds=seconds, geometry_seconds=geometry['seconds'],
                        predicates=geometry['exact_sign_predicates'],
                        all_geometry_arrays_bitwise_identical=True,
                        geometry_path=str(target),
                        geometry_sha256=sha256(target.read_bytes()).hexdigest(),
                        arrays_sha256=sha256(target.with_suffix('.npz').read_bytes()).hexdigest()))
            rows = [row for row in report['rows'] if row['state'] == name]
            off = np.median([row['seconds'] for row in rows if not row['cache_signs']])
            on = np.median([row['seconds'] for row in rows if row['cache_signs']])
            report['states'][-1].update(uncached_median_seconds=float(off), cached_median_seconds=float(on),
                                        geometry_speedup=float(off/on))
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--repeats', type=int, default=2)
    run(**vars(parser.parse_args()))

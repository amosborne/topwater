"""A real moving chart retains original cold provenance and strict replay code."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest

from topwater import InitialState
from .experiments.free_surface import boat_contact_reaction as snapshot
from .experiments.free_surface.boat_geometry_template import evaluate_template
from .reference_boat import boundary_snapshot


def test_moving_chart_replay_owns_cold_sources_and_rejects_changed_warm_code(tmp_path, monkeypatch):
    case = json.loads(Path(__file__).with_name('data').joinpath('boat-state-337ms.json').read_text())
    state = np.array(case['state'])
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    geometry = tmp_path/'moving-geometry.json'
    _, chart, _ = boundary_snapshot(geometry, pose=pose, parent_eta=state[13:28],
                            parent_psi=state[28:], export_template=True)
    arrays, _ = evaluate_template(chart, pose=pose, parent_eta=state[13:28], parent_psi=state[28:])
    pilot = tmp_path/'pilot.json'
    pilot.write_text(json.dumps(dict(boundary_refinement=dict(levels=0))))
    audit = tmp_path/'snapshot.json'
    saved = dict(state=state, **{name: arrays[name] for name in ('points', 'faces', 'free_nodes')})
    np.savez(audit.with_suffix('.npz'), **saved)
    archive, current = tmp_path/'snapshot-source', tmp_path/'current'
    archive.mkdir(); current.mkdir()
    hashes = {}
    for name in snapshot.DEPENDENCIES:
        contents = Path(snapshot.__file__).with_name(name).read_bytes()
        hashes[name] = sha256(contents).hexdigest()
        (archive/name).write_bytes(contents)
        (current/name).write_bytes(contents)
    audit.write_text(json.dumps(dict(complete=True, source_sha256=hashes,
        input_sha256={str(path): sha256(path.read_bytes()).hexdigest()
            for path in (pilot, geometry, geometry.with_suffix('.npz'))},
        chart_identity=dict(path=str(geometry.with_suffix('.npz')), node_chart_key=str(chart['node_chart_key'])))))
    # Use isolated source identities; imported runtime replay remains unchanged.
    monkeypatch.setattr(snapshot, '__file__', str(current/'boat_contact_reaction.py'))
    snapshot.load_snapshot(audit)
    cold = current/'boat_world_boundary.py'
    cold.write_text(cold.read_text()+'\n# Cold producer changed after the owned snapshot.\n')
    template = current/'boat_geometry_template.py'
    original_bytes = template.read_bytes()
    original = template.read_text()
    changed_export = original.replace('Cold exact witness export', 'Changed cold witness export', 1)
    assert changed_export != original
    template.write_text(changed_export)
    report, _, replay, _, _ = snapshot.load_snapshot(audit)
    assert set(report['geometry_replay_verification']) == {'boat_world_boundary.py', 'boat_geometry_template.py'}
    np.testing.assert_array_equal(replay['points'], saved['points'])
    changed_definition = changed_export.replace('def export_template(registry, meta, arrays, critical_vertices):',
        'def export_template(registry, meta, arrays, critical_vertices, *, changed_default=True):', 1)
    assert changed_definition != changed_export
    template.write_text(changed_definition)
    with pytest.raises(ValueError, match='outside cold exporter'):
        snapshot.load_snapshot(audit)
    changed_warm = changed_export.replace('    kind, ids = label',
        '    raise RuntimeError("changed warm geometry")\n    kind, ids = label', 1)
    assert changed_warm != changed_export
    template.write_text(changed_warm)
    with pytest.raises(ValueError, match='outside cold exporter'):
        snapshot.load_snapshot(audit)
    template.write_text(changed_export)
    (archive/'boat_geometry_template.py').write_text('corrupted original source')
    with pytest.raises(ValueError, match='archive is absent or changed'):
        snapshot.load_snapshot(audit)
    (archive/'boat_geometry_template.py').write_bytes(original_bytes)
    saved['faces'] = saved['faces'][:, ::-1]
    np.savez(audit.with_suffix('.npz'), **saved)
    with pytest.raises(ValueError, match='saved geometry differs: faces'):
        snapshot.load_snapshot(audit)

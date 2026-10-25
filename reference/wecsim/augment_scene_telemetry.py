"""Augment accepted scene with producer-defined telemetry; never alter geometry."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[2]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def array_sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def augment(config):
    from topwater.visualization import SceneSequence
    for p,digest in config['input_sha256'].items():assert sha(ROOT/p)==digest,p
    target=ROOT/config['augmented_scene'];assert not target.exists(),'Preserve prior augmented input'
    source=SceneSequence.load(ROOT/config['scene'])
    accepted=read(ROOT/config['motion_acceptance'])
    assert accepted['decision']=='accepted_for_visualization' and accepted['reviewed_by']
    assert accepted['motion_sha256']==sha(ROOT/config['motion_mat'])
    native=loadmat(ROOT/config['motion_mat'],simplify_cells=True)
    times=np.asarray(native['time']).ravel();q=np.asarray(native['q'])
    rows=np.searchsorted(times,source.time_s)
    assert np.all(rows<len(times)) and np.array_equal(times[rows],source.time_s),'Scene times must exactly match native samples'
    assert q.ndim==2 and q.shape[1]==6 and np.isfinite(q[rows]).all()
    telemetry=dict(roll_deg=np.rad2deg(q[rows,3]),pitch_deg=np.rad2deg(q[rows,4]),heave_m=q[rows,2].copy())
    names=['time_s','hull_points','hull_faces','water_points','water_faces']
    arrays={name:getattr(source,name).copy() for name in names}
    metadata=dict(source.metadata,telemetry_units=dict(roll_deg='deg',pitch_deg='deg',heave_m='m'),telemetry_basis='Producer supplied native RzRyRx roll/pitch; known-neutral-CG vertical displacement, not firstclip-frame rezero',telemetry_source_scene_sha256=sha(ROOT/config['scene']),telemetry_source_motion_sha256=sha(ROOT/config['motion_mat']))
    augmented=SceneSequence(**arrays,metadata=metadata,telemetry=telemetry)
    target.parent.mkdir(parents=True,exist_ok=True);augmented.save(target)
    restored=SceneSequence.load(target)
    identity={name:dict(original=array_sha(getattr(source,name)),augmented=array_sha(getattr(restored,name)),dtype=str(getattr(source,name).dtype),shape=list(getattr(source,name).shape)) for name in names}
    assert all(getattr(source,name).dtype==getattr(restored,name).dtype and np.array_equal(getattr(source,name),getattr(restored,name)) for name in names),'Geometry/time/connectivity changed'
    assert all(np.array_equal(restored.telemetry[k],v) for k,v in telemetry.items())
    assert all(sha(ROOT/p)==digest for p,digest in config['input_sha256'].items())
    return dict(complete=True,source_scene_sha256=sha(ROOT/config['scene']),augmented_scene_sha256=sha(target),geometry_time_bit_exact=True,array_identity=identity,native_rows_zero_based=rows.tolist(),telemetry_names=list(telemetry),source_time_s=[float(source.time_s[0]),float(source.time_s[-1])])

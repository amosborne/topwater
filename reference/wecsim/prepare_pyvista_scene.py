"""WEC/cache adapter only: world meshes for the generic visualization library.
Uses the accepted incident coefficient bank; never runs hydrodynamics/dynamics.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from calm_buoy_train import CalmBuoyTrain
ROOT=Path(__file__).resolve().parents[2]

def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')

def prepare_scene(config,output):
    from topwater.visualization import SceneSequence
    started=time.monotonic();output=Path(output)
    assert not output.exists(),'Preserve prior prepared cache; use load/reuse in runner'
    for name,digest in config['input_sha256'].items():assert sha(ROOT/name)==digest,name
    output.mkdir(parents=True)
    with np.load(ROOT/config['accepted_cache']) as saved:cache={k:saved[k].copy() for k in saved.files}
    physical=read(ROOT/config['physical_case']);heading=read(ROOT/config['heading_file'])
    accepted=read(ROOT/config['motion_acceptance']);review=read(ROOT/config['cache_acceptance'])
    assert accepted['decision']=='accepted_for_visualization' and accepted['reviewed_by']
    assert accepted['packet_sha256']==sha(ROOT/config['packet_npz'])
    assert accepted['heading_sha256']==sha(ROOT/config['heading_file'])
    assert review['decision']=='preview_pass' and review['surface_cache_sha256']==sha(ROOT/config['accepted_cache'])
    selected=np.flatnonzero((cache['time']>=config['clip_start_s']-1e-9)&(cache['time']<=config['clip_end_s']+1e-9))
    times=cache['time'][selected];assert len(times)==config['expected_frames']
    assert np.max(abs(np.diff(times)-1/config['fps']))<1e-8
    hull=cache['hull'][selected].copy();centers=cache['center'][selected]
    geometry=physical['minimum_hull'];faces=np.asarray([f for i,f in enumerate(geometry['triangular_faces']) if i not in geometry['nonphysical_closure_faces']],dtype=np.int64)
    train=CalmBuoyTrain(ROOT/config['packet_npz'])
    R=np.asarray(train.data['global_from_boat_rotation']);translation=np.asarray(train.data['global_from_boat_translation_m'])
    assert np.array_equal(R,np.asarray(heading['global_from_boat_rotation']))
    half=config['display_halfwidth_m'];step=config['outer_spacing_m']
    outer=np.arange(-half,half+step/2,step)
    innerx=cache['dx'][0];innery=cache['dy'][:,0]
    x=np.r_[outer[outer<innerx.min()-1e-9],innerx,outer[outer>innerx.max()+1e-9]]
    y=np.r_[outer[outer<innery.min()-1e-9],innery,outer[outer>innery.max()+1e-9]]
    xx,yy=np.meshgrid(x,y);offset=np.c_[xx.ravel(),yy.ravel()]
    native_centers=(np.c_[centers,np.zeros(len(centers))]-translation)@R
    native_offset=offset@R[:2,:2]
    directions=np.c_[np.cos(train.theta),np.sin(train.theta)]
    radius=float(np.linalg.norm(offset,axis=1).max())
    lo=float((native_centers[:,:2]@directions.T).min()-radius)
    hi=float((native_centers[:,:2]@directions.T).max()+radius)
    sgrid=np.linspace(lo,hi,int(np.ceil((hi-lo)/config['profile_spacing_m']))+1)
    eta=np.zeros((len(times),len(offset)))
    active=np.flatnonzero(train.data['carrier_amplitudes_m']);batch=config['component_batch']
    for first in range(0,len(active),batch):
        assert time.monotonic()-started<config['maximum_prepare_seconds'],'Displayfield preparation budget exhausted'
        ids=active[first:first+batch];profiles=train.profiles(ids,sgrid,times)
        for j,component in enumerate(ids):
            projected_offset=native_offset@directions[component]
            projected_center=native_centers[:,:2]@directions[component]
            for frame in range(len(times)):eta[frame]+=np.interp(projected_offset+projected_center[frame],sgrid,profiles[j,frame])
    eta=eta.reshape(len(times),len(y),len(x))
    ix=np.searchsorted(x,innerx);iy=np.searchsorted(y,innery)
    generated_inner=eta[:,iy[:,None],ix[None,:]].copy()
    eta[:,iy[:,None],ix[None,:]]=cache['eta'][selected]
    assert np.array_equal(eta[:,iy[:,None],ix[None,:]],cache['eta'][selected])
    points=np.empty((len(times),len(offset),3));points[:,:,:2]=centers[:,None,:]+offset[None,:,:];points[:,:,2]=eta.reshape(len(times),-1)
    indices=np.arange(len(offset)).reshape(len(y),len(x));a=indices[:-1,:-1];b=indices[:-1,1:];c=indices[1:,:-1];d=indices[1:,1:]
    waterfaces=np.concatenate([np.stack([a,b,c],axis=-1).reshape(-1,3),np.stack([c,b,d],axis=-1).reshape(-1,3)])
    metadata=dict(units='m',coordinate_system='right-handed-z-up',world_axes='global east,north,up',source='accepted WEC world meshes, same frozen incident train',mean_heading_deg=float(np.rad2deg(heading['mean_box_yaw_global_rad'])),preparation_controls={k:config[k] for k in ['clip_start_s','clip_end_s','expected_frames','display_halfwidth_m','outer_spacing_m','profile_spacing_m']},input_sha256=config['input_sha256'],display_only=True)
    scene=SceneSequence(time_s=times,hull_points=hull,hull_faces=faces,water_points=points,water_faces=waterfaces,metadata=metadata)
    scene.save(output/'scene.npz');write(output/'scene.json',metadata)
    report=dict(complete=True,elapsed_s=time.monotonic()-started,frames=len(times),time_s=[float(times[0]),float(times[-1])],water_nodes=[len(y),len(x)],source_cache_indices=selected.tolist(),central_accepted_heights_bit_identical=True,wide_profile_vs_original_near_max_error_m=float(abs(generated_inner-cache['eta'][selected]).max()),profile_interpolation_bound_m=float(np.diff(sgrid).max()**2*np.sum(abs(train.C)*train.k[None,None,:]**2)/8),field='Wider DISPLAYONLY same accepted coefficient field; no regenerated sea/force/motion',scene_sha256=sha(output/'scene.npz'))
    assert all(sha(ROOT/name)==digest for name,digest in config['input_sha256'].items())
    write(output/'preparation.json',report);print(json.dumps(report),flush=True)
    return scene

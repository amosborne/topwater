"""Full accepted reference cache + ONLY missing border nodes + exact native HUD."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.io import loadmat
from calm_buoy_train import CalmBuoyTrain
ROOT=Path(__file__).resolve().parents[2]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')

def prepare_full_scene(c):
    from topwater.visualization import SceneSequence,RenderSettings
    started=time.perf_counter();target=ROOT/c['scene_output'];assert not target.exists()
    assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    accepted=read(ROOT/c['motion_acceptance']);review=read(ROOT/c['cache_acceptance'])
    assert accepted['decision']=='accepted_for_visualization' and accepted['motion_sha256']==sha(ROOT/c['motion_mat'])
    assert accepted['packet_sha256']==sha(ROOT/c['packet_npz'])
    assert review['decision']=='preview_pass' and review['surface_cache_sha256']==sha(ROOT/c['surface_cache'])
    with np.load(ROOT/c['surface_cache']) as z:cache={k:z[k].copy() for k in z.files}
    times=cache['time'];hull=cache['hull'];centers=cache['center'];ix=cache['dx'][0];iy=cache['dy'][:,0]
    assert len(times)==c['frames'] and times[0]==0 and times[-1]==60
    assert np.max(abs(np.diff(times)-1/c['fps']))<1e-8
    physical=read(ROOT/c['physical_case']);geometry=physical['minimum_hull']
    faces=np.asarray([f for i,f in enumerate(geometry['triangular_faces']) if i not in geometry['nonphysical_closure_faces']],dtype=np.int64)
    # IDENTICAL fixed topdown view formula in the frozen shared renderer.
    viewcenter=((hull.min(axis=1)+hull.max(axis=1))/2).mean(axis=0)[:2]
    half=float(abs(hull[:,:,:2]-viewcenter).max())*RenderSettings().framing
    needed=np.max(abs(viewcenter-centers),axis=0)+half
    steps=np.array([np.diff(ix).max(),np.diff(iy).max()]);oldhalf=np.array([abs(ix).max(),abs(iy).max()])
    pad=np.where(needed<=oldhalf,0,np.ceil((needed-oldhalf)/steps).astype(int)+1)
    x=ix[0]+np.arange(-pad[0],len(ix)+pad[0])*steps[0];y=iy[0]+np.arange(-pad[1],len(iy)+pad[1])*steps[1]
    assert np.array_equal(x[pad[0]:pad[0]+len(ix)],ix) and np.array_equal(y[pad[1]:pad[1]+len(iy)],iy)
    xx,yy=np.meshgrid(x,y);offset=np.c_[xx.ravel(),yy.ravel()]
    border=np.ones(xx.shape,bool);border[pad[1]:pad[1]+len(iy),pad[0]:pad[0]+len(ix)]=False
    eta=np.empty((len(times),len(y),len(x)));eta[:,pad[1]:pad[1]+len(iy),pad[0]:pad[0]+len(ix)]=cache['eta']
    if border.any():
        train=CalmBuoyTrain(ROOT/c['packet_npz']);R=train.data['global_from_boat_rotation'];translation=train.data['global_from_boat_translation_m']
        native_centers=(np.c_[centers,np.zeros(len(centers))]-translation)@R
        native_offsets=offset[border.ravel()]@R[:2,:2];directions=np.c_[np.cos(train.theta),np.sin(train.theta)]
        radius=np.linalg.norm(native_offsets,axis=1).max();lo=float((native_centers[:,:2]@directions.T).min()-radius);hi=float((native_centers[:,:2]@directions.T).max()+radius)
        sgrid=np.linspace(lo,hi,int(np.ceil((hi-lo)/c['profile_spacing_m']))+1);added=np.zeros((len(times),int(border.sum())))
        active=np.flatnonzero(train.data['carrier_amplitudes_m']);batch=c['component_batch']
        for first in range(0,len(active),batch):
            assert time.perf_counter()-started<c['maximum_prepare_seconds'],'Border preparation budget exhausted'
            ids=active[first:first+batch];profiles=train.profiles(ids,sgrid,times)
            for j,component in enumerate(ids):
                projection=native_offsets@directions[component];origin=native_centers[:,:2]@directions[component]
                for frame in range(len(times)):added[frame]+=np.interp(projection+origin[frame],sgrid,profiles[j,frame])
        eta[:,border]=added
    assert np.array_equal(eta[:,pad[1]:pad[1]+len(iy),pad[0]:pad[0]+len(ix)],cache['eta'])
    lower=viewcenter-half;upper=viewcenter+half
    assert np.all(centers[:,0]+x[0]<lower[0]) and np.all(centers[:,0]+x[-1]>upper[0])
    assert np.all(centers[:,1]+y[0]<lower[1]) and np.all(centers[:,1]+y[-1]>upper[1])
    points=np.empty((len(times),len(offset),3));points[:,:,:2]=centers[:,None,:]+offset[None,:,:];points[:,:,2]=eta.reshape(len(times),-1)
    idx=np.arange(len(offset)).reshape(len(y),len(x));a=idx[:-1,:-1];b=idx[:-1,1:];d=idx[1:,:-1];e=idx[1:,1:]
    waterfaces=np.concatenate([np.stack([a,b,d],axis=-1).reshape(-1,3),np.stack([d,b,e],axis=-1).reshape(-1,3)])
    native=loadmat(ROOT/c['motion_mat'],simplify_cells=True);native_time=np.asarray(native['time']).ravel();q=np.asarray(native['q']);rows=np.searchsorted(native_time,times)
    assert np.array_equal(rows,np.arange(0,len(native_time),4)) and np.array_equal(native_time[rows],times)
    telemetry=dict(roll_deg=np.rad2deg(q[rows,3]),pitch_deg=np.rad2deg(q[rows,4]),heave_m=q[rows,2].copy())
    scene=SceneSequence(time_s=times,hull_points=hull,hull_faces=faces,water_points=points,water_faces=waterfaces,telemetry=telemetry,metadata=dict(units='m',coordinate_system='right-handed-z-up',source='Acceptedfull751cache; onlymissingphysicalfieldDISPLAYborder evaluated',input_sha256=c['input_sha256'],telemetry_units=dict(roll_deg='deg',pitch_deg='deg',heave_m='m')))
    target.parent.mkdir(parents=True,exist_ok=True);scene.save(target)
    restored=SceneSequence.load(target)
    assert np.array_equal(restored.time_s,times) and np.array_equal(restored.hull_points,hull)
    old_water=np.stack([cache['dx'][None]+centers[:,None,None,0],cache['dy'][None]+centers[:,None,None,1],cache['eta']],axis=-1)
    new_inner=restored.water_points.reshape(len(times),len(y),len(x),3)[:,pad[1]:pad[1]+len(iy),pad[0]:pad[0]+len(ix)]
    assert np.array_equal(old_water,new_inner) and all(np.array_equal(restored.telemetry[k],v) for k,v in telemetry.items())
    assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    return dict(complete=True,preparation_seconds=time.perf_counter()-started,frames=len(times),source_time_s=[float(times[0]),float(times[-1])],grid_nodes_xy=[len(x),len(y)],new_border_nodes_per_frame=int(border.sum()),central_accepted_water_xyz_bit_exact=True,hull_time_bit_exact=True,telemetry_exact_native_rows=True,fixed_view_center_xy=viewcenter.tolist(),fixed_view_half_m=half,water_half_xy_m=[float(x[-1]),float(y[-1])],all_fixed_view_corners_covered=True,scene_sha256=sha(target),input_sha256=c['input_sha256'])

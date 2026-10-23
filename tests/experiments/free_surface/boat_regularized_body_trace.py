"""Fluid-side body-normal trace audit of a recorded homogeneous mode.

Offsets h and h/2 are explicit one-sided convergence probes, not solver-point
jitter or a change to the body geometry/BC. The prescribed normal velocity is0.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull
from topwater import SimulationInputs
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .physical_p1_patch import PhysicalP1Patch
from .compact_p1_maps import body_maps
from .panel_gradient_reference import checked_influence


def owned(path):
    report=json.loads(path.read_text())
    if not report['complete'] or sha256(path.with_suffix('.npz').read_bytes()).hexdigest()!=report['arrays_sha256']:
        raise ValueError('complete owned arrays required: '+str(path))
    return report


def run(args,report):
    mode=owned(args.mode);water=owned(args.water);phase=owned(args.phase)
    if (mode['parents']['water']!=sha256(args.water.read_bytes()).hexdigest()
            or mode['parents']['phase']!=sha256(args.phase.read_bytes()).hexdigest()
            or water['parent_phase_sha256']!=sha256(args.phase.read_bytes()).hexdigest()):
        raise ValueError('mode/water/phase ownership differs')
    for name in ('tests/experiments/free_surface/physical_p1_patch.py',
                 'tests/experiments/free_surface/compact_p1_trace.py',
                 'tests/experiments/free_surface/compact_p1_maps.py'):
        path=Path(name)
        if sha256(path.read_bytes()).hexdigest()!=mode['producer_sha256'].get(str(path)):
            raise ValueError('mode source/extension owner changed: '+name)
    contact_path=args.phase.with_suffix('.source')/'contact.json'
    if sha256(contact_path.read_bytes()).hexdigest()!=phase['parent_contact_sha256']:
        raise ValueError('phase/contact ownership differs')
    contact=json.loads(contact_path.read_text());audit_path=Path(contact['parents']['audit']['path']);audit=owned(audit_path)
    if sha256(audit_path.read_bytes()).hexdigest()!=contact['parents']['audit']['sha256']:
        raise ValueError('audit ownership changed')
    inputs=SimulationInputs.from_json((audit_path.with_suffix('.source')/'inputs.json').read_text())
    if sha256(inputs.to_json().encode()).hexdigest()!=mode['input_sha256']:
        raise ValueError('physical inputs changed')
    with np.load(audit_path.with_suffix('.npz')) as z: triangles=z['triangles'].copy()
    if sha256(triangles.tobytes()).hexdigest()!=mode['panel_sha256']:
        raise ValueError('body panel ownership changed')
    with np.load(args.water.with_suffix('.npz')) as z: x=z['x'].copy();y=z['y'].copy();active=z['active_nodes'].copy()
    with np.load(args.phase.with_suffix('.npz')) as z: polygon=z['polygon'].copy()
    patch=PhysicalP1Patch(x,y,polygon);np.testing.assert_array_equal(patch.active_nodes,active)
    with np.load(args.mode.with_suffix('.npz')) as z:
        eta=z['mode_0_eta'].copy();coefficient=z['mode_0_chi'].copy();sigma=z['mode_0_sigma'].copy();growth=complex(z['eigenvalues'][0])
    scale=args.eta_rms/np.sqrt(np.vdot(eta,patch.mass@eta).real/float(patch.mass.sum()))
    full=patch.expand(coefficient*scale);sigma*=scale
    panels=TrianglePanels(triangles,-1);solid=ConvexHull(triangles.reshape(-1,3))

    def evaluate(points,normals):
        p,n=body_maps(patch,points,normals)
        value=p@full;normal=n@full
        for first in range(0,len(points),128):
            last=min(first+128,len(points));target=points[first:last]
            real=checked_influence(target,panels);image=checked_influence(target*[1.,1.,-1.],panels)
            value[first:last]+=(real.single_layer+image.single_layer)@sigma
            gradient=real.gradient_single_layer+image.gradient_single_layer*[1.,1.,-1.]
            normal[first:last]+=np.einsum('mni,mi->mn',gradient,normals[first:last])@sigma
        return value,normal

    report.update(input_sha256=mode['input_sha256'],panel_sha256=mode['panel_sha256'],
        eta_mode_rms_m=args.eta_rms,mode_growth_real_per_s=growth.real,mode_growth_imag_rad_per_s=growth.imag,
        collocation_normal_residual_at_mode_scale=mode['modes'][0]['body_normal_residual_max']*scale,sets=[])
    saved={}
    for name,order in [('centroids',0),('independent Gauss2',2),('independent Gauss4',4)]:
        if order:
            points,weights=triangle_samples(triangles,order);parents=np.repeat(np.arange(len(panels)),len(weights))
            weight=(panels.areas[:,None]*weights).ravel();base=points.reshape(-1,3)
        else:
            base=panels.centroids.copy();parents=np.arange(len(panels));weight=panels.areas.copy()
        normals=panels.normals[parents];diameters=panels.diameters[parents]
        offset=args.offset_fraction*diameters
        rises=normals[:,2]<0
        offset[rises]=np.minimum(offset[rises],.25*abs(base[rises,2])/(-normals[rises,2]))
        record=dict(name=name,points=len(base),offsets=[]);report['sets'].append(record);normal_history=[]
        for factor in (1.,.5):
            step=offset*factor;target=base-step[:,None]*normals
            if np.any(target[:,2]>=0) or np.any(np.all(target==base,axis=1)):
                raise ValueError('one-sided offset is unresolved or crosses the source plane')
            separation=np.max(target@solid.equations[:,:3].T+solid.equations[:,3],axis=1)
            if np.any(separation<=0): raise ValueError('normal trace probe did not leave the solid')
            print(json.dumps(dict(set=name,offset_factor=factor,start=True)),flush=True);started=perf_counter()
            value,normal=evaluate(target,normals);normal_history.append(normal)
            entry=dict(offset_factor=factor,seconds=perf_counter()-started,
                offset_min_m=float(step.min()),offset_max_m=float(step.max()),
                offset_over_diameter_min=float((step/diameters).min()),
                offset_over_diameter_max=float((step/diameters).max()),
                normal_real_rms_m_per_s=float(np.sqrt(weight@(normal.real**2)/weight.sum())),
                normal_imag_rms_m_per_s=float(np.sqrt(weight@(normal.imag**2)/weight.sum())),
                normal_complex_rms_m_per_s=float(np.sqrt(weight@abs(normal)**2/weight.sum())),
                normal_max_m_per_s=float(abs(normal).max()),
                net_normal_volume_flux_real_m3_per_s=float((weight@normal).real),
                net_normal_volume_flux_imag_m3_per_s=float((weight@normal).imag),
                potential_complex_rms_m2_per_s=float(np.sqrt(weight@abs(value)**2/weight.sum())),
                green_body_boundary_work_J=float(WATER_DENSITY*np.real(weight@(np.conj(value)*normal))),
                pressure_work_into_fluid_W=float(WATER_DENSITY*np.real(weight@(np.conj(growth*value)*normal))))
            record['offsets'].append(entry)
            saved[f'{order}_{factor}_normal']=normal;saved[f'{order}_{factor}_potential']=value
            print(json.dumps(dict(set=name,**entry)),flush=True)
        record['offset_halving_normal_rms_change_m_per_s']=float(np.sqrt(weight@abs(normal_history[0]-normal_history[1])**2/weight.sum()))
    np.savez(args.output.with_suffix('.npz'),**saved)
    report.update(complete=True,qualified=False,scope=__doc__,arrays_sha256=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest(),
        limits=['One-sided convergence is empirical, not a certified trace error interval; all offsets and body geometry remain explicit',
            'Body normal velocity should be0 throughout the physical surface; centroid residuals alone do not enforce that condition',
            'Green work has energy units; pressure-rate diagnostic is a declared sesquilinear envelope with phi_t=lambda*phi, no time-average factor',
            'A large off-centroid leak implicates boundary consistency but does not exclude far-closure/projection work defects',
            'No stabilization, tolerance relaxation, nonlinear6DOF trajectory or throughput qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('mode','water','phase','output'): parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--eta-rms',type=float,default=.001)
    parser.add_argument('--offset-fraction',type=float,default=1e-4)
    args=parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists() or not 0<args.eta_rms<=.01
            or not 0<args.offset_fraction<.01): parser.error('fresh archive and positive bounded mode/offset scales required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,
        parents={name:sha256(getattr(args,name).read_bytes()).hexdigest() for name in ('mode','water','phase')})
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))

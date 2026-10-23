"""One bounded source-trace refinement on the SAME actual waterline test hat.

S-only exterior integration avoids a fine full-patch weak kernel. The fixed
body solution, physical support, far trace and external reference are retained.
This is a spatial diagnostic, not a new motion-error or throughput criterion.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.constants import GRAVITY
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from tests.capytaine_reference import body_from_triangles
from .boat_hos_continuation import window
from .boat_local_radiation_trace import phases,amplitude,far_trace_spectrum
from .boat_local_trace_enrichment import tensor_values
from .boat_weak_contact_rows import owned,triangle_rule,fourier_values
from .cut_hat_geometry import cut_hat
from .plane_p1_field import PlaneP1Field
from .panel_gradient_reference import checked_influence


def run(args,report):
    import capytaine as cpt
    contact=owned(args.contact)
    paths={name:Path(record['path']) for name,record in contact['parents'].items()}
    parents={}
    for name,path in paths.items():
        if sha256(path.read_bytes()).hexdigest()!=contact['parents'][name]['sha256']:
            raise ValueError('contact input ownership changed: '+name)
        parents[name]=owned(path)
    snapshot,field,audit=(parents[name] for name in ('snapshot','field','audit'))
    inputs=SimulationInputs.from_json((paths['audit'].with_suffix('.source')/'inputs.json').read_text())
    identity=sha256(inputs.to_json().encode()).hexdigest()
    if any(p['input_sha256']!=identity for p in (contact,snapshot,field,audit)):
        raise ValueError('physical inputs differ')
    with np.load(paths['snapshot'].with_suffix('.npz')) as z:
        x=z['patch_x'].copy();y=z['patch_y'].copy();delta=z['residual'].copy();density=z['enriched_density'].copy()
    with np.load(paths['audit'].with_suffix('.npz')) as z:
        triangles=z['triangles'].copy();external_sources=z['sources_0'].copy()
    panel_hash=sha256(triangles.tobytes()).hexdigest()
    if any(p['panel_sha256']!=panel_hash for p in (contact,snapshot,field,audit)):
        raise ValueError('body panels differ')
    with np.load(args.contact.with_suffix('.npz')) as z:
        polygon=z['polygon'].copy();previous_geometry=z['geometry_0_exterior'].copy()
    case=contact['cases'][0]
    if case['name']!='straight waterline': raise ValueError('owned straight-waterline diagnostic required')
    hx,hy=float(x[1]-x[0]),float(y[1]-y[0])
    geometry=cut_hat(case['centre'],hx,hy,polygon)
    np.testing.assert_array_equal(geometry['exterior']['triangles'],previous_geometry)
    np.testing.assert_allclose(geometry['physical_hat_mass'],case['physical_hat_mass'],atol=0.,rtol=0.)
    mass=geometry['physical_hat_mass'];cg=inputs.initial_state.position
    grid=snapshot['grid'];length=snapshot['length_m'];origin=cg[:2]-length/2
    with np.load(paths['field'].with_suffix('.npz')) as z: trace=z['surface_0'].reshape(grid,grid,2)
    beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    inner=np.max(abs(triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam;outer=inner+beam
    far=far_trace_spectrum(trace,origin,length,cg[:2],inner,beam)
    kx,ky=np.meshgrid(2*np.pi*np.fft.fftfreq(grid,d=length/grid),2*np.pi*np.fft.fftfreq(grid,d=length/grid))
    far_vertical=far*np.hypot(kx,ky)[:,:,None]
    omega=2*np.pi*field['frequency_hz'];panels=TrianglePanels(triangles,-1)
    external_density=phases(-external_sources/(-1j*omega))
    settings=dict(audit['rows'][0]['green_settings']);settings.pop('green_function')
    cache=Path('artifacts/capytaine-surface-audit-cache');cache.mkdir(exist_ok=True)
    body=body_from_triangles(triangles,cg,inputs.body.mass)
    solver=cpt.BEMSolver(green_function=cpt.Delhommeau(**settings,tabulation_cache_dir=str(cache)))
    external=[SimpleNamespace(body=body,sources=external_sources[:,c],free_surface=0.,water_depth=np.inf,
        encounter_wavenumber=omega**2/GRAVITY,forward_speed=0.) for c in range(2)]

    def reference(points,velocity=False):
        columns=[]
        function=solver.compute_velocity if velocity else solver.compute_potential
        with threadpool_limits(1):
            for result in external:
                values=np.concatenate([function(np.column_stack((points[first:first+256],np.zeros(len(points[first:first+256])))),result)
                    for first in range(0,len(points),256)])
                columns.append((values[:,2] if velocity else values)/(-1j*omega))
        return np.column_stack(columns)

    quadrature={}
    for order in args.orders:
        ep,ew,ev,eg,_=triangle_rule(geometry['exterior'],order)
        line=geometry['waterline'];nodes,weights=leggauss(order);t=(nodes+1)/2;weights=weights/2
        lp=(line['segments'][:,:1]+t[None,:,None]*(line['segments'][:,1:]-line['segments'][:,:1])).reshape(-1,2)
        lv=((1-t)*line['values'][:,:1]+t*line['values'][:,1:]).ravel()
        lw=(np.linalg.norm(line['segments'][:,1]-line['segments'][:,0],axis=1)[:,None]*weights).ravel()*lv
        ln=np.repeat(line['hole_outward_normals'],order,axis=0)
        source=np.zeros(4);transplanted=np.zeros(4)
        for first in range(0,len(ep),128):
            last=min(first+128,len(ep));p=np.column_stack((ep[first:last],np.zeros(last-first)))
            mapping=2*checked_influence(p,panels).gradient_single_layer[:,:,2]
            source+=(ew[first:last]*ev[first:last])@(mapping@density)
            transplanted+=(ew[first:last]*ev[first:last])@(mapping@external_density)
        r=(ew*ev)@reference(ep,velocity=True)
        quadrature[order]=dict(ep=ep,ew=ew,ev=ev,eg=eg,lp=lp,lw=lw,ln=ln,
            far=(ew*ev)@fourier_values(far_vertical,ep,origin,length),source=source,
            reference=r,reference_lift=r-amplitude(transplanted[None])[0],
            reference_trace=reference(ep),far_trace=fourier_values(far,ep,origin,length))
    report.update(input_sha256=identity,panel_sha256=panel_hash,physical_hat_mass=mass,
        physical_test_centre=case['centre'],physical_test_spacing=[hx,hy],rows=[])
    saved={}
    for refinement in args.refinements:
        scale=2**refinement
        fx=np.linspace(x[0],x[-1],(len(x)-1)*scale+1)
        fy=np.linspace(y[0],y[-1],(len(y)-1)*scale+1)
        np.testing.assert_array_equal(fx[::scale],x)
        np.testing.assert_array_equal(fy[::scale],y)
        if refinement==0: values=delta
        else:
            xx,yy=np.meshgrid(fx,fy)
            exact=phases(reference(np.column_stack((xx.ravel(),yy.ravel())))).reshape(len(fy),len(fx),4)
            taper=window(xx-cg[0],inner[0],outer[0])*window(yy-cg[1],inner[1],outer[1])
            values=(exact-tensor_values(far,fx,fy,origin,length))*taper[:,:,None]
            values[[0,-1],:]=0.;values[:,[0,-1]]=0.
            np.testing.assert_allclose(values[::scale,::scale],delta,atol=2e-12,rtol=2e-12)
        owner=PlaneP1Field(fx,fy,values,args.maximum_pairs)
        row=dict(refinement=refinement,source_triangles=len(owner.panels),
            source_spacing=[float(fx[1]-fx[0]),float(fy[1]-fy[0])],quadrature=[])
        report['rows'].append(row)
        for order in args.orders:
            print(json.dumps(dict(refinement=refinement,order=order,start=True)),flush=True)
            started=perf_counter();q=quadrature[order]
            eo=owner.point_owners(q['ep']);lo=owner.point_owners(q['lp'])
            a=owner.vector_potential(q['ep'],eo)
            area=2*np.einsum('m,mi,mic->c',q['ew'],q['eg'],a)
            a=owner.vector_potential(q['lp'],lo)
            waterline=2*np.einsum('m,mi,mic->c',q['lw'],q['ln'],a)
            local=area+waterline;lift=amplitude((q['far']+local)[None])[0]
            total=amplitude((q['far']+local+q['source'])[None])[0]
            nodal=values.reshape(-1,4)[owner.faces[eo,0]]
            point_trace=nodal+np.einsum('mi,mic->mc',q['ep']-owner.panels.triangles[eo,0,:2],owner.gradient[eo])
            trace_error=amplitude(point_trace+q['far_trace'])-q['reference_trace']
            previous=next((v for v in case['quadrature'] if v['order']==order),None)
            entry=dict(order=order,seconds=perf_counter()-started,area=area.tolist(),waterline=waterline.tolist(),
                local_moment=local.tolist(),lift_moment_real=lift.real.tolist(),lift_moment_imag=lift.imag.tolist(),
                inferred_reference_lift_real=q['reference_lift'].real.tolist(),inferred_reference_lift_imag=q['reference_lift'].imag.tolist(),
                lift_defect_absolute=abs(lift-q['reference_lift']).tolist(),
                total_defect_absolute=abs(total-q['reference']).tolist(),
                trace_weighted_rms_error=np.sqrt((q['ew']*q['ev'])@abs(trace_error)**2/mass).tolist(),
                old_full_minus_hidden_difference=abs(amplitude((local-np.array(previous['corrected_local']))[None])[0]).tolist()
                    if previous is not None and refinement==0 else None,source_pairs_so_far=owner.pairs)
            row['quadrature'].append(entry)
            saved[f'refinement_{refinement}_q{order}_local']=local
            saved[f'refinement_{refinement}_x']=fx;saved[f'refinement_{refinement}_y']=fy
            saved[f'refinement_{refinement}_trace']=values
            np.savez(args.output.with_suffix('.npz'),**saved)
            report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
            print(json.dumps(dict(refinement=refinement,**entry)),flush=True)
    report.update(complete=True,qualified=False,scope=__doc__,peak_process_rss_bytes=peak_rss_bytes(),
        limits=['Fixed body density is intentionally not re-solved: isolate the trace lift, not a new coupled BVP',
            'Target quadrature is fixed on the same physical cut cells; order sensitivity is diagnostic, not a certified bound',
            'Inferred reference lift includes external finite-panel and native kernel approximation uncertainty',
            'No time integration, body-motion accuracy, throughput or baseline promotion claim'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('contact','output'): parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--refinements',type=int,nargs='+',default=[0,1])
    parser.add_argument('--orders',type=int,nargs='+',default=[8,16])
    parser.add_argument('--maximum-pairs',type=int,default=900_000_000)
    args=parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists() or not args.refinements
            or min(args.refinements)<0 or max(args.refinements)>1 or not args.orders
            or min(args.orders)<2 or max(args.orders)>16 or args.maximum_pairs<1):
        parser.error('fresh archive, relative refinements0..1, orders2..16 and positive source-pair budget required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.contact,archive/'contact.json')
    report=dict(complete=False,producer_sha256=hashes,parent_contact_sha256=sha256(args.contact.read_bytes()).hexdigest())
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))

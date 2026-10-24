"""Data-only sampled projection-error allocation; not a mesh marking rule."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import owned


def checkpoint(args,report,arrays):
    np.savez(args.output.with_suffix('.npz'),**arrays)
    report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False))


def run(args,report):
    projection,gal=owned(args.projection),owned(args.galerkin)
    if projection['parents']['galerkin']!=sha256(args.galerkin.read_bytes()).hexdigest():
        raise ValueError('projection and Galerkin parent differ')
    for name in ('tests/laplace_reference.py','src/topwater/rankine.py','src/topwater/constants.py'):
        relative=Path(name)
        if sha256(relative.read_bytes()).hexdigest()!=projection['producer_sha256'].get(str(relative)):
            raise ValueError('sample/physical-geometry reconstruction owner differs: '+name)
    with np.load(args.galerkin.with_suffix('.npz')) as z:
        triangles,vertices,faces,lift,center=(z[name].copy() for name in
            ('triangles','vertices','faces','odd_lift','world_cg'))
    with np.load(args.projection.with_suffix('.npz')) as z:
        np.testing.assert_array_equal(triangles,z['physical_triangles'])
        np.testing.assert_array_equal(center,z['world_cg'])
        values={f'q{q}_{name}':z[f'q{q}_{name}'].copy() for q in (4,8) for name in
            ('projection_coefficients','reference_phi','projection_minus_reference_port')}
    nf=len(triangles)
    np.testing.assert_array_equal(vertices[faces[:nf]],triangles)
    panels=TrianglePanels(triangles,-1)  # Geometry only; no influence kernel.
    seam=np.any(triangles[:,:,2]==0,axis=1)
    area=float(panels.areas.sum())
    arrays=dict(triangles=triangles,world_cg=center,facet_areas=panels.areas,seam_touching=seam)
    report.update(complete=False,qualified=False,physical_facets=nf,
        seam_touching_facets=int(seam.sum()),seam_touching_area_share=float(panels.areas[seam].sum()/area),
        seam_definition='Any original physical triangle vertex exactly at z=0; diagnostic only',
        area_concentration_scope='Area of minimal-count descending-score facet set; NOT a minimum-area optimization',
        orders=[])
    sets,distribution={},{}
    for order in (4,8):
        started=perf_counter()
        bary,weight=triangle_samples(np.eye(3)[None],order);bary=bary[0]
        coefficients=values[f'q{order}_projection_coefficients']
        nodal=(lift@coefficients)[faces[:nf]]
        projected=np.einsum('qj,fjc->fqc',bary,nodal)
        reference=values[f'q{order}_reference_phi'].reshape(nf,len(weight),6)
        error=projected-reference
        w=panels.areas[:,None]*weight
        reference_norm2=np.sum(w[:,:,None]*reference**2,axis=(0,1))
        if np.any(reference_norm2<=0):
            raise ValueError('equal-six normalization unresolved for a zero reference trace column; no mode exclusion')
        energy=np.sum(w[:,:,None]*error**2,axis=1)
        normalized=energy/reference_norm2
        score=normalized.mean(axis=1)
        total=float(score.sum())
        relative=np.sqrt(normalized.sum(axis=0))
        expected=next(row for row in projection['orders'] if row['order']==order)['relative_trace_L2']['projection_minus_reference']
        np.testing.assert_allclose(relative,expected,rtol=1e-12,atol=1e-14)
        ranking=np.argsort(-score,kind='stable')
        cumulative=np.cumsum(score[ranking])
        concentration={}
        for fraction in (.5,.8,.9):
            count=int(np.searchsorted(cumulative,fraction*total,side='left')+1) if total else 0
            selected=ranking[:count]
            tag=str(int(fraction*100))
            sets[order,tag]=set(int(v) for v in selected)
            arrays[f'q{order}_diagnostic_top_{tag}_facets']=selected
            concentration[tag]=dict(facets=count,facet_share=count/nf,
                area_m2=float(panels.areas[selected].sum()),area_share=float(panels.areas[selected].sum()/area),
                captured_score_share=float(score[selected].sum()/total) if total else 0.)
        points=np.einsum('qj,fji->fqi',bary,triangles)
        normals=panels.normals[:,None,:]
        T=np.concatenate((np.broadcast_to(normals,points.shape),np.cross(points-center,normals)),axis=2)
        load=WATER_DENSITY*np.einsum('fq,fqi,fqj->fij',w,T,error)
        saved_port=values[f'q{order}_projection_minus_reference_port']
        parity=float(abs(load.sum(axis=0)-saved_port).max())
        scale=float((WATER_DENSITY*np.einsum('fq,fqi,fqj->ij',w,abs(T),abs(error))).max())
        allowance=128*np.finfo(float).eps*(nf*len(weight)+6)*scale
        if parity>allowance:
            raise ValueError('per-facet signed load contributions differ from owned aggregate port')
        diag=np.diagonal(load,axis1=1,axis2=2)
        groups={}
        for name,mask in (('all',np.ones(nf,dtype=bool)),('seam_touching',seam),('not_seam_touching',~seam)):
            signed=diag[mask].sum(axis=0)
            absolute=abs(diag[mask]).sum(axis=0)
            groups[name]=dict(facets=int(mask.sum()),area_share=float(panels.areas[mask].sum()/area),
                sampled_score_share=float(score[mask].sum()/total) if total else 0.,
                signed_diagonal_load=signed.tolist(),absolute_diagonal_load_sum=absolute.tolist(),
                cancellation_fraction=[float(1-abs(a)/b) if b else None for a,b in zip(signed,absolute,strict=True)],
                signed_full_load=load[mask].sum(axis=0).tolist())
        distribution[order]=score/total if total else np.zeros(nf)
        arrays.update({f'q{order}_reference_trace_norm_squared':reference_norm2,
            f'q{order}_facet_projection_error_squared':energy,
            f'q{order}_facet_normalized_six_errors':normalized,
            f'q{order}_facet_equal_six_score':score,f'q{order}_facet_signed_load':load})
        row=dict(order=order,seconds=perf_counter()-started,relative_projection_L2=relative.tolist(),
            total_equal_six_score=total,concentration=concentration,groups=groups,
            aggregate_load_parity_max=parity,aggregate_load_arithmetic_allowance=allowance)
        report['orders'].append(row)
        checkpoint(args,report,arrays)
        print(json.dumps(dict(order=order,concentration=concentration,seam_score_share=groups['seam_touching']['sampled_score_share'])),flush=True)
    overlap={}
    for tag in ('50','80','90'):
        a,b=sets[4,tag],sets[8,tag]
        overlap[tag]=dict(intersection=len(a&b),union=len(a|b),jaccard=len(a&b)/len(a|b) if a|b else 1.)
    report.update(complete=True,topset_q4_q8_overlap=overlap,
        normalized_score_distribution_overlap=float(np.minimum(distribution[4],distribution[8]).sum()),
        process_peak_rss=peak_rss_bytes(),limits=[
            'Sampled coarse P0 source control is not converged physical truth; its original triangulation edges may create artificial roughness',
            'Do not adapt future hull meshes to coarse-source artifacts; no automatic marking/refinement is implemented',
            'Equal-six normalized trace error is not a certified motion/load estimator or unavoidable load floor',
            'Stable input-order tie breaks are only for reporting minimal-count concentration, not a geometric region rule',
            'Signed load cancellation is retained; trace-error concentration alone does not locate load accuracy',
            'q4/q8 sensitivity is reported before any mesh recommendation; no new source kernel, solve, mesh or quadrature targets'])
    checkpoint(args,report,arrays)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('projection','galerkin','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output/source required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest()
        for name in ('projection','galerkin')})
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))

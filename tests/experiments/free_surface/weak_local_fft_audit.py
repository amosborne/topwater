"""Weak compact-P1 DtN: independent assembly, indexing, scale and actual-size cost.

No hull-hole deletion, body-source work or curved/nonlinear evolution is inferred
from this full-plane compact weak operator. Independent Poisson work is separate.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np

from tests.benchmark_rankine import peak_rss_bytes
from tests.laplace_reference import triangle_samples
from .compact_p1_trace import CompactP1Trace
from .compact_p1_weak import WeakPatch,lag_kernel,mass_matrix
from .cuda_operator import cp,_pairs
from .boat_benchmark_contract import measurement_environment


def dense_reference(shape,hx,hy,order):
    """Assemble every physical triangle pair, then scatter interior gradients.

    Shares analytic scalar panel primitive, not hat/lag aggregation or FFT.
    This checks geometry/assembly; independent Poisson work checks the formula.
    """
    ny,nx=shape
    patch=CompactP1Trace(np.arange(nx+2)*hx,np.arange(ny+2)*hy,np.zeros((ny+2,nx+2)))
    faces=patch.faces;panels=patch.panels;nf=len(panels);nv=(nx+2)*(ny+2)
    points,weights=triangle_samples(panels.triangles,order)
    opposite=np.roll(panels.triangles,-2,axis=1)-np.roll(panels.triangles,-1,axis=1)
    gradient=np.cross(panels.geometric_normals[:,None],opposite)/(2*panels.areas[:,None,None])
    geometry=[cp.asarray(np.ascontiguousarray(a)) for a in
        (panels.triangles,panels.geometric_normals,panels.edge_outward,panels.tangents,
         panels.edge_lengths,panels.areas,panels.centroids,panels.normal_sign,
         faces.astype(np.int64),gradient,points,weights)]
    v=cp.empty((nf,nf),dtype=cp.float64);unused=cp.zeros((nf,nv),dtype=cp.float64)
    _pairs(((nf*nf+127)//128,),(128,),(*geometry,np.int32(nf),np.int32(nv),np.int32(len(weights)),v,unused))
    raw=cp.asnumpy(v);v=(raw+raw.T)/2
    lookup=np.full(nv,-1,dtype=int)
    lookup.reshape(ny+2,nx+2)[1:-1,1:-1]=np.arange(nx*ny).reshape(shape)
    bx=np.zeros((nf,nx*ny));by=np.zeros_like(bx);mass=np.zeros((nx*ny,nx*ny))
    for t in range(nf):
        for a in range(3):
            i=lookup[faces[t,a]]
            if i<0: continue
            bx[t,i]=gradient[t,a,0];by[t,i]=gradient[t,a,1]
            for b in range(3):
                j=lookup[faces[t,b]]
                if j>=0: mass[i,j]+=panels.areas[t]*(2 if a==b else 1)/12
    return 2*(bx.T@v@bx+by.T@v@by),mass,float(np.max(abs(raw-raw.T)))


def dense_lags(kernel,shape):
    ny,nx=shape;yy,xx=np.indices(shape);y=yy.ravel();x=xx.ravel()
    return kernel[y[:,None]-y[None,:]+ny-1,x[:,None]-x[None,:]+nx-1]


def benchmark(patch,values):
    for _ in range(3): patch.apply(cp.asarray(values))
    cp.cuda.get_current_stream().synchronize()
    times=[];answer=None
    for _ in range(10):
        before=perf_counter();answer=patch.apply(cp.asarray(values))
        cp.cuda.get_current_stream().synchronize();times.append(perf_counter()-before)
    return dict(median_host_synchronized_seconds=float(np.median(times)),
        minimum_seconds=min(times),maximum_seconds=max(times),sample_count=len(times),
        weak_energy=float(cp.sum(cp.asarray(values)*answer).item()))


def run(args,report):
    with np.load(args.residual) as z:
        residual=z['residual'].copy();x=z['patch_x'].copy();y=z['patch_y'].copy()
    shape=(7,8);hx=.31;hy=.23
    patch=WeakPatch(shape,hx,hy)
    kernel=patch.kernel;lag_dense=dense_lags(kernel,shape)
    references=[]
    for order in [8,16,32]:
        matrix,mass,asymmetry=dense_reference(shape,hx,hy,order)
        references.append(matrix)
        report.setdefault('dense_quadrature',[]).append(dict(order=order,
            maximum_kernel_matrix_difference=float(np.max(abs(matrix-lag_dense))),
            raw_triangle_integral_asymmetry=asymmetry))
    reference=references[-1]
    np.testing.assert_allclose(patch.mass.toarray(),mass,atol=2e-16,rtol=3e-14)
    assert np.linalg.eigvalsh(reference).min()>0
    yy,xx=np.indices(shape)
    actual=residual[len(residual)//2:len(residual)//2+shape[0],
                    residual.shape[1]//2:residual.shape[1]//2+shape[1],0]
    boundary=np.zeros(shape);boundary[0,0]=1.
    cases=[('actual residual samples',actual),('near boundary hat',boundary),
           ('smooth compact trace',np.sin(np.pi*(xx+1)/(shape[1]+1))*np.sin(np.pi*(yy+1)/(shape[0]+1))),
           ('interior plateau with zero exterior',np.ones(shape))]
    report['actions']=[]
    for name,values in cases:
        expected=lag_dense@values.ravel()
        actual=cp.asnumpy(patch.apply(cp.asarray(values))).ravel()
        np.testing.assert_allclose(actual,expected,atol=2e-14,rtol=2e-12)
        error=np.linalg.norm(actual-reference@values.ravel())
        norm=np.linalg.norm(reference@values.ravel())
        assert error<2e-4*norm
        energy=float(values.ravel()@actual)
        assert energy>0
        report['actions'].append(dict(case=name,fft_vs_toeplitz_max_difference=float(np.max(abs(actual-expected))),
            versus_independent_assembly_relative_l2=float(error/norm),energy=energy))
    scaled,scaled_diagnostics=lag_kernel(shape,2*hx,2*hy)
    np.testing.assert_allclose(scaled,2*kernel,atol=2e-12,rtol=1e-10)
    np.testing.assert_allclose(mass_matrix(shape,2*hx,2*hy).toarray(),4*mass,atol=2e-16,rtol=3e-14)
    far,_=lag_kernel((13,26),hx,hy)
    far_value=far[12+12,25+25]
    radius=np.hypot(25*hx,12*hy)
    point_limit=-(hx*hy)**2/(2*np.pi*radius**3)
    assert abs(far_value/point_limit-1)<.02
    report.update(small_shape=shape,small_spacing=[hx,hy],kernel_diagnostics=patch.diagnostics,
        small_cost=benchmark(patch,cases[2][1]),
        isotropic_scale_max_difference=float(np.max(abs(scaled-2*kernel))),
        far_value=far_value,far_point_limit=point_limit,far_relative_difference=far_value/point_limit-1,
        plateau_energy=report['actions'][-1]['energy'])
    saved=dict(kernel=kernel,dense_reference=reference,mass=mass)
    if args.full_patch:
        full_shape=(len(y)-2,len(x)-2);spacing=[float(x[1]-x[0]),float(y[1]-y[0])]
        before=perf_counter();full=WeakPatch(full_shape,*spacing)
        cp.cuda.get_current_stream().synchronize()
        report['actual_patch']=dict(shape=full_shape,spacing=spacing,
            constructor_seconds=perf_counter()-before,diagnostics=full.diagnostics,
            application=benchmark(full,residual[1:-1,1:-1,0]),
            scope='Flat full-plane weak operator component; excludes physical hole/cross maps/mass solve/curved geometry')
        saved['actual_kernel']=full.kernel
    np.savez(args.output.with_suffix('.npz'),**saved)
    environment=measurement_environment(cp,'uncontrolled',
        'Serial numerical lease, BLAS/OMP1; component timings only, no implementation-speed ratio',
        'Synchronized component output; independent audits and archive outside application timing')
    environment.update(timing_scope='Weak-operator construction plus static verification; repeated apply timings include host upload, finite check and stream synchronization; no integration',
        startup_treatment='Actual-size constructor follows small-case shared-kernel/context warmup; includes new-size FFT plan work; repeated applies warmed; no cold-start claim')
    report.update(qualified=False,scope=__doc__,
        measurement_environment=environment,
        arrays_sha256=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest(),
        peak_process_rss_bytes=peak_rss_bytes(),ending_cuda_pool_reserved_bytes=cp.get_default_memory_pool().total_bytes(),
        memory_scope='Whole verification process peakRSS; ending CUDA pool is not peakVRAM. Dense references excluded from operator ownership figures.',
        limits=['K is weak stiffness, not a nodal normal velocity or M^-1K',
            'Full-plane compact patch only; no physical hull-hole/contact or global/local/body cross coupling',
            'Dense assembly shares inner scalar Rankine primitive; independent Poisson-work check remains separate',
            'No nonlinear trajectory or whole-simulation speed/accuracy qualification'])
    report['complete']=True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--residual',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--full-patch',action='store_true')
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output/source archive required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,residual_sha256=sha256(args.residual.read_bytes()).hexdigest())
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error),
            weak_diagnostics=getattr(error,'weak_diagnostics',None));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))

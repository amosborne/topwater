"""Fixed, graded triangular patches for validation-only geometry preflight."""
import argparse
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from .common import (prepare,read_owned,current_sources,PHYSICAL,VALIDATORS,
    digest,physical_identity,save,load_contract,study_preflight)
from .reference_mesh import (graded_mesh,verify_surface,farthest_first,facets,
    COORD_ATOL)


def lattice_triangles(n):
    """One fixed diagonal pattern on the integer barycentric triangular grid."""
    triangles=[]
    for i in range(n):
        for j in range(n-i):
            triangles.append(((i,j),(i+1,j),(i,j+1)))
            if i+j<n-1:
                triangles.append(((i+1,j),(i+1,j+1),(i,j+1)))
    return np.asarray(triangles)


def canonical_connectivity(triangles):
    return sorted(tuple(sorted(map(tuple,triangle))) for triangle in triangles)


def verify_patch(parent,children,n,spec):
    """Reject fallback remeshing, displaced nodes and changed diagonal patterns."""
    if np.shape(children)!=(n*n,3,3):
        raise ValueError('transfinite patch must contain exactly N squared triangles')
    a,b,c=parent
    basis=np.column_stack((b-a,c-a))
    uv=np.linalg.lstsq(basis,(children-a).reshape(-1,3).T,rcond=None)[0].T.reshape(-1,3,2)
    integer=np.rint(uv*n).astype(int)
    if np.max(abs(uv*n-integer))>spec['barycentric_lattice_absolute']:
        raise ValueError('patch nodes differ from the declared affine barycentric lattice')
    if canonical_connectivity(integer)!=canonical_connectivity(lattice_triangles(n)):
        raise ValueError('transfinite patch connectivity differs from the fixed diagonal pattern')
    parent_area=np.linalg.norm(np.cross(b-a,c-a))/2
    vector=np.cross(children[:,1]-children[:,0],children[:,2]-children[:,0])/2
    if np.any(vector@np.cross(b-a,c-a)<=0):
        raise ValueError('negative patch Jacobian')
    np.testing.assert_allclose(np.linalg.norm(vector,axis=1),parent_area/(n*n),
        rtol=spec['child_area_relative'],atol=spec['child_area_absolute_m2'])
    reconstructed=a+uv[...,0,None]*(b-a)+uv[...,1,None]*(c-a)
    if np.max(abs(reconstructed-children))>COORD_ATOL:
        raise ValueError('transfinite patch nodes leave the original physical plane')
    return float(np.max(abs(uv*n-integer)))


def shape_statistics(triangles):
    a,b,c=np.asarray(triangles).transpose(1,0,2)
    lengths=np.stack((np.linalg.norm(b-c,axis=1),np.linalg.norm(a-c,axis=1),np.linalg.norm(a-b,axis=1)),axis=1)
    area=np.linalg.norm(np.cross(b-a,c-a),axis=1)/2
    angles=[]
    for i in range(3):
        x,y,z=lengths[:,i],lengths[:,(i+1)%3],lengths[:,(i+2)%3]
        angles.append(np.rad2deg(np.arccos(np.clip((y*y+z*z-x*x)/(2*y*z),-1,1))))
    aspect=lengths.max(axis=1)**2/(2*area)
    if not np.isfinite(aspect).all() or np.any(area<=0):
        raise ValueError('nonfinite or degenerate panel shape')
    order=np.argsort(aspect)
    weighted=np.interp([0,.5,.95,.99,1],np.cumsum(area[order])/area.sum(),aspect[order])
    return dict(panels=len(triangles),area_m2=float(area.sum()),
        minimum_angle_degrees=float(np.min(angles)),
        edge_length_min_m=float(lengths.min()),edge_length_max_m=float(lengths.max()),
        longest_edge_over_altitude_quantiles=np.quantile(aspect,[0,.5,.95,.99,1]).tolist(),
        area_weighted_aspect_quantiles=weighted.tolist(),
        area_fraction_aspect_above_10_50_100=[float(area[aspect>limit].sum()/area.sum()) for limit in (10,50,100)])


def facet_statistics(original,triangles):
    rows=[]
    for i,facet in enumerate(facets(original)):
        n=facet['normal'];normal=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
        normal/=np.linalg.norm(normal,axis=1)[:,None]
        selected=(np.max(abs((triangles-facet['polygon'][0])@n),axis=1)<=COORD_ATOL)&(normal@n>1-1e-8)
        rows.append(dict(facet=i,normal=n.tolist(),**shape_statistics(triangles[selected])))
    if sum(row['panels'] for row in rows)!=len(triangles):
        raise ValueError('quality statistics do not cover every physical panel exactly once')
    return rows


def transfinite_refine(seed,n,spec):
    """Gmsh triangular transfinite patches, one affine map per frozen seed face."""
    import gmsh
    if gmsh.__version__!=spec['gmsh_version']:
        raise ValueError('pinned Gmsh version required')
    half=seed[:len(seed)//2]
    if len(seed)%2 or np.any(half[:,:,1]<-COORD_ATOL):
        raise ValueError('verified reflected seed ordering required')
    vertices,ids=np.unique(half.reshape(-1,3),axis=0,return_inverse=True)
    faces=ids.reshape(-1,3)
    gmsh.initialize(['reference-transfinite','-nopopup'])
    try:
        for key,value in {'General.Terminal':0,'General.NumThreads':1,
                'Mesh.MaxNumThreads1D':1,'Mesh.MaxNumThreads2D':1,
                'Mesh.TransfiniteTri':1,'Mesh.Smoothing':0,'Mesh.Optimize':0,
                'Mesh.OptimizeNetgen':0,'Mesh.MeshSizeFromPoints':0,
                'Mesh.MeshSizeFromCurvature':0,'Mesh.MeshSizeExtendFromBoundary':0}.items():
            gmsh.option.setNumber(key,value)
        gmsh.model.add('frozen-affine-patches')
        points=[gmsh.model.geo.addPoint(*point) for point in vertices]
        curves={};surfaces=[]
        for face in faces:
            sides=[]
            for a,b in zip(face,np.roll(face,-1),strict=True):
                key=tuple(sorted((int(a),int(b))))
                if key not in curves:
                    curves[key]=gmsh.model.geo.addLine(points[key[0]],points[key[1]])
                sides.append(curves[key] if a<b else -curves[key])
            loop=gmsh.model.geo.addCurveLoop(sides)
            surfaces.append(gmsh.model.geo.addPlaneSurface([loop]))
        gmsh.model.geo.synchronize()
        for curve in curves.values():
            gmsh.model.mesh.setTransfiniteCurve(curve,n+1,'Progression',1.)
        for surface,face in zip(surfaces,faces,strict=True):
            gmsh.model.mesh.setTransfiniteSurface(surface,'Left',[points[i] for i in face])
            gmsh.model.mesh.setSmoothing(2,surface,0)
        gmsh.model.mesh.generate(2)
        tags,xyz,_=gmsh.model.mesh.getNodes();xyz=np.asarray(xyz).reshape(-1,3)
        lookup={int(tag):i for i,tag in enumerate(tags)}
        pieces=[];maximum_error=0.
        for surface,parent in zip(surfaces,half,strict=True):
            types,_,nodes=gmsh.model.mesh.getElements(2,surface)
            if list(types)!=[2]:
                raise ValueError('structured linear triangles required')
            t=xyz[np.array([lookup[int(tag)] for tag in nodes[0]]).reshape(-1,3)]
            if np.mean(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])@np.cross(parent[1]-parent[0],parent[2]-parent[0]))<0:
                t=t[:,[0,2,1]]
            maximum_error=max(maximum_error,verify_patch(parent,t,n,spec))
            pieces.append(t)
        half=np.concatenate(pieces)
        half[:,:,1][abs(half[:,:,1])<COORD_ATOL]=0
        reflected=half[:,[0,2,1]].copy();reflected[:,:,1]*=-1
        result=farthest_first(np.concatenate((half,reflected)))
        curve_rows=[]
        for (a,b),tag in curves.items():
            _,coords,_=gmsh.model.mesh.getNodes(1,tag,includeBoundary=True)
            coords=np.asarray(coords).reshape(-1,3)
            delta=vertices[b]-vertices[a]
            position=np.sort((coords-vertices[a])@delta/(delta@delta))
            if len(position)!=n+1 or np.max(abs(position-np.linspace(0,1,n+1)))>spec['barycentric_lattice_absolute']:
                raise ValueError('shared curve subdivision differs from declared uniform affine map')
            curve_rows.append(dict(vertices=[a,b],positions=position.tolist()))
        return result,dict(subdivisions=n,half_patch_count=len(faces),
            maximum_integer_lattice_error=maximum_error,curves=curve_rows,
            seed_vertices=vertices.tolist(),seed_faces=faces.tolist())
    finally:
        gmsh.finalize()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol',type=Path,default=Path('validation/transfinite-preflight.json'))
    parser.add_argument('--contract',type=Path,default=Path('validation/canonical-boat.json'))
    parser.add_argument('--frozen',type=Path,required=True)
    parser.add_argument('--physical-reference',type=Path,default=Path('artifacts/validation-reset/frozen.json'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();spec=json.loads(args.protocol.read_text());contract=load_contract(args.contract)
    if (spec['format_version']!=1 or spec['gmsh_version']!='4.15.2'
            or spec['triangle_vertex_order']!='farthest-first'
            or spec['smoothing_steps']!=0 or spec['mesh_optimization'] is not False
            or spec['richardson_extrapolation'] is not False
            or spec['barycentric_lattice_absolute']!=1e-7
            or spec['child_area_relative']!=1e-8 or spec['child_area_absolute_m2']!=1e-14
            or spec['seed_maximum_size_to_facet_width']!=20
            or len(spec['subdivisions'])!=4 or spec['subdivisions']!=sorted(set(spec['subdivisions']))
            or any(type(n)!=int or n<2 for n in spec['subdivisions'])):
        raise ValueError('unsupported fixed-topology geometry protocol')
    report=prepare(args.output,args.protocol,'transfinite-geometry-preflight');arrays={};started=perf_counter()
    try:
        frozen,data=read_owned(args.frozen,digest(args.contract),'freeze',passed=True)
        current_sources(frozen,PHYSICAL+VALIDATORS)
        if physical_identity(frozen,data)!=frozen['physical_identity_sha256']:
            raise ValueError('frozen input array identity differs')
        input_file=args.frozen.with_suffix('.source')/'inputs.json'
        if digest(input_file)!=frozen['input_sha256']:
            raise ValueError('frozen input JSON identity differs')
        shutil.copyfile(input_file,args.output.with_suffix('.source')/'inputs.json')
        shutil.copyfile(args.contract,args.output.with_suffix('.source')/'canonical-contract.json')
        shutil.copyfile(args.frozen,args.output.with_suffix('.source')/'bound-freeze.json')
        physical=('wetted_triangles','world_cg','orientation','dry_mass','hydrostatic_stiffness','body_mass','hull_length','hull_beam')
        arrays.update({key:data[key] for key in physical})
        reference_header=json.loads(args.physical_reference.read_text())
        reference,original=read_owned(args.physical_reference,reference_header['contract_sha256'],'freeze',passed=True)
        reference_input=args.physical_reference.with_suffix('.source')/'inputs.json'
        if digest(reference_input)!=reference['input_sha256']:
            raise ValueError('physical reference input JSON identity differs')
        equality={key:bool(data[key].dtype==original[key].dtype and np.array_equal(data[key],original[key])) for key in physical}
        equality['input_json_bytes']=input_file.read_bytes()==reference_input.read_bytes()
        if not all(equality.values()):
            raise ValueError('original frozen physical problem changed')
        for source,label in ((args.frozen,'bound-freeze'),(args.physical_reference,'physical-reference')):
            for suffix in ('.json','.npz'):
                shutil.copyfile(source.with_suffix(suffix),args.output.with_suffix('.source')/(label+suffix))
            shutil.copytree(source.with_suffix('.source'),args.output.with_suffix('.source')/(label+'.source'))
        report.update(input_sha256=frozen['input_sha256'],frozen_manifest_sha256=digest(args.frozen),
            canonical_contract_sha256=digest(args.contract),physical_reference_manifest_sha256=digest(args.physical_reference),
            physical_equality_checks=equality,
            frozen_physical_identity_sha256=frozen['physical_identity_sha256'],
            qualification_status='unqualified; geometry only; no hydrodynamic calls',
            complete=False,levels=[])
        save(args.output,report,arrays)
        seed,settings=graded_mesh(data['wetted_triangles'],spec['seed_bulk_size_m'],spec['seed_edge_size_m'],growth=spec['seed_edge_growth'],
            maximum_size_to_facet_width=spec['seed_maximum_size_to_facet_width'])
        arrays['seed_triangles']=seed
        report['seed']=dict(settings=settings,quality=shape_statistics(seed),facets=facet_statistics(data['wetted_triangles'],seed))
        save(args.output,report,arrays)
        resources=contract['reference_numerics']
        report['study_plan']=study_preflight(len(seed),spec['subdivisions'],resources['maximum_estimated_workspace_bytes'],
            panel_counts=[len(seed)*n*n for n in spec['subdivisions']])
        for n in spec['subdivisions']:
            if perf_counter()-started>resources['maximum_elapsed_seconds']:
                raise TimeoutError('geometry preflight elapsed budget exhausted')
            t,mesh=transfinite_refine(seed,n,spec)
            arrays[f'triangles_{n}']=t
            proof=verify_surface(data['wetted_triangles'],t)
            if len(t)!=len(seed)*n*n:
                raise ValueError('full family count differs from frozen seed times N squared')
            report['levels'].append(dict(**mesh,geometry=proof,quality=shape_statistics(t),
                facets=facet_statistics(data['wetted_triangles'],t),complete=True))
            save(args.output,report,arrays)
            if perf_counter()-started>resources['maximum_elapsed_seconds']:
                raise TimeoutError('geometry preflight elapsed budget exceeded after native meshing')
        report.update(complete=True,passed=True)
    except (Exception,KeyboardInterrupt) as error:
        report.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    report['elapsed_seconds']=perf_counter()-started;save(args.output,report,arrays)
    print(json.dumps(dict(complete=report['complete'],passed=report['passed'],output=str(args.output))))
    return 0 if report['complete'] and report['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())

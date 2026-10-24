"""Validation-only planar-facet meshes; never change the frozen wet surface."""
import numpy as np
from scipy.spatial import ConvexHull


COORD_ATOL=1e-11
GEOM_ATOL=1e-12
GEOM_RTOL=1e-10


def farthest_first(triangles):
    """Keep winding/geometry, make Capytaine's first-vertex radius conservative."""
    triangles=np.asarray(triangles)
    centers=triangles.mean(axis=1)
    first=np.argmax(np.linalg.norm(triangles-centers[:,None],axis=2),axis=1)
    return np.take_along_axis(triangles,((first[:,None]+np.arange(3))%3)[:,:,None],axis=1)


def moments(triangles):
    a,b,c=np.asarray(triangles).transpose(1,0,2)
    vector=np.cross(b-a,c-a)/2
    area=np.linalg.norm(vector,axis=1)
    total=a+b+c
    second=(np.einsum('fi,fj->fij',total,total)
        +sum(np.einsum('fi,fj->fij',v,v) for v in (a,b,c)))/12
    return np.concatenate(([area.sum()],vector.sum(axis=0),
        np.sum(area[:,None]*total/3,axis=0),
        np.sum(area[:,None,None]*second,axis=0).ravel()))


def half_surface(triangles):
    """Clip to world y>=0; reflection is permitted only after area-moment proof."""
    result=[]
    for triangle in triangles:
        polygon=[]
        for a,b in zip(triangle,np.roll(triangle,-1,axis=0),strict=True):
            if a[1]>=0:
                polygon.append(a)
            if (a[1]<0<b[1]) or (b[1]<0<a[1]):
                q=a+(b-a)*(-a[1]/(b[1]-a[1]));q[1]=0
                polygon.append(q)
        for j in range(1,len(polygon)-1):
            t=np.array([polygon[0],polygon[j],polygon[j+1]])
            if np.linalg.norm(np.cross(t[1]-t[0],t[2]-t[0]))>GEOM_ATOL:
                result.append(t)
    return np.asarray(result)


def facets(triangles):
    """Merge convex coplanar facets only after proving their area moments."""
    triangles=np.asarray(triangles)
    normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    groups=[]
    for i,(t,n) in enumerate(zip(triangles,normals,strict=True)):
        for group in groups:
            first=group[0]
            if (np.max(abs(n-normals[first]))<COORD_ATOL
                    and np.max(abs((t-triangles[first,0])@normals[first]))<COORD_ATOL):
                group.append(i);break
        else:
            groups.append([i])
    result=[]
    for group in groups:
        source=triangles[group];n=normals[group[0]]
        points=source.reshape(-1,3)
        # Eliminate roundoff duplicates without changing a retained coordinate.
        _,indices=np.unique(np.round(points,11),axis=0,return_index=True)
        points=points[np.sort(indices)]
        origin=points[0];u=points[np.argmax(np.linalg.norm(points-origin,axis=1))]-origin
        u/=np.linalg.norm(u);v=np.cross(n,u)
        xy=np.column_stack(((points-origin)@u,(points-origin)@v))
        hull=ConvexHull(xy);polygon=points[hull.vertices]
        fan=np.array([[polygon[0],polygon[j],polygon[j+1]] for j in range(1,len(polygon)-1)])
        np.testing.assert_allclose(moments(fan),moments(source),atol=GEOM_ATOL,rtol=GEOM_RTOL)
        result.append(dict(polygon=polygon,normal=n,source=source,origin=origin,u=u,v=v))
    return result


def verify_surface(original,triangles):
    """Per-facet containment + moments and conforming edges; reject lost strips."""
    original=np.asarray(original);triangles=np.asarray(triangles)
    groups=facets(original);assignment=np.full(len(triangles),-1)
    for i,facet in enumerate(groups):
        polygon,n=facet['polygon'],facet['normal']
        planar=np.max(abs((triangles-polygon[0])@n),axis=1)<=COORD_ATOL
        inside=np.ones(len(triangles),bool)
        for a,b in zip(polygon,np.roll(polygon,-1,axis=0),strict=True):
            inward=np.cross(n,b-a);inward/=np.linalg.norm(inward)
            inside&=np.min((triangles-a)@inward,axis=1)>=-COORD_ATOL
        selected=planar&inside
        if not selected.any():
            raise ValueError('missing physical facet')
        if np.any(assignment[selected]>=0):
            raise ValueError('panel assigned to multiple physical facets')
        assignment[selected]=i
        t=triangles[selected]
        if np.any(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])@n<=0):
            raise ValueError('inverted/degenerate physical panel')
        np.testing.assert_allclose(moments(t),moments(facet['source']),atol=GEOM_ATOL,rtol=GEOM_RTOL)
    if np.any(assignment<0):
        raise ValueError('panel is not contained in original physical surface')
    # Every interior edge has exactly two opposite incidences. Every open edge
    # is on the original waterline, which is the only allowed open boundary.
    vertices,indices=np.unique(np.round(triangles.reshape(-1,3),11),axis=0,return_inverse=True)
    faces=indices.reshape(-1,3)
    directed=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    edges,inverse,counts=np.unique(np.sort(directed,axis=1),axis=0,return_inverse=True,return_counts=True)
    if np.any(counts>2):
        raise ValueError('overlapping/nonmanifold surface')
    open_edges=edges[counts==1]
    if np.max(abs(vertices[open_edges][...,2]))>COORD_ATOL:
        raise ValueError('gap or nonconforming physical feature edge')
    signed=np.where(directed[:,0]<directed[:,1],1,-1)
    orientation=np.bincount(inverse,weights=signed,minlength=len(edges))
    if np.any(orientation[counts==2]!=0):
        raise ValueError('inconsistent shared-edge winding')
    np.testing.assert_allclose(moments(triangles),moments(original),atol=GEOM_ATOL,rtol=GEOM_RTOL)
    return dict(facets=len(groups),panels=len(triangles),waterline_edges=len(open_edges),
        maximum_moment_difference=float(np.max(abs(moments(triangles)-moments(original)))))


def facet_width(polygon,normal):
    edges=np.roll(polygon,-1,axis=0)-polygon
    inward=np.cross(normal,edges)
    inward/=np.linalg.norm(inward,axis=1)[:,None]
    return float(np.min(np.ptp(polygon@inward.T,axis=0)))


def graded_mesh(original,bulk_size,edge_size,*,growth=.3,maximum_size_to_facet_width=None):
    """Pinned Gmsh, half-surface + reflected panels, full physical entry solves."""
    import gmsh
    if gmsh.__version__!='4.15.2':
        raise ValueError('pinned Gmsh 4.15.2 required')
    if not 0<edge_size<=bulk_size or not 0<growth<=1:
        raise ValueError('positive bulk/edge sizes and bounded growth required')
    if maximum_size_to_facet_width is not None and (
            not np.isfinite(maximum_size_to_facet_width) or maximum_size_to_facet_width<=0):
        raise ValueError('positive finite facet-width size factor required')
    half=half_surface(original)
    reflection=half[:,[0,2,1]].copy();reflection[:,:,1]*=-1
    verify_surface(original,np.concatenate((half,reflection)))
    polygons=facets(half)
    gmsh.initialize(['reference-mesh','-nopopup'])
    try:
        for name,value in {'General.Terminal':0,'General.NumThreads':1,
                'Mesh.MaxNumThreads1D':1,'Mesh.MaxNumThreads2D':1,
                'Mesh.Algorithm':6,'Mesh.MeshSizeFromPoints':0,
                'Mesh.MeshSizeFromCurvature':0,'Mesh.MeshSizeExtendFromBoundary':0,
                'Mesh.RandomFactor':1e-9}.items():
            gmsh.option.setNumber(name,value)
        gmsh.model.add('frozen-wet-half')
        points={};lines={};surfaces=[]
        for facet in polygons:
            ids=[]
            for point in facet['polygon']:
                key=tuple(np.round(point,11))
                if key not in points:
                    points[key]=gmsh.model.geo.addPoint(*point)
                ids.append(points[key])
            curves=[]
            for a,b in zip(ids,ids[1:]+ids[:1],strict=True):
                key=tuple(sorted((a,b)))
                if key not in lines:
                    lines[key]=gmsh.model.geo.addLine(*key)
                curves.append(lines[key] if a<b else -lines[key])
            loop=gmsh.model.geo.addCurveLoop(curves)
            surfaces.append(gmsh.model.geo.addPlaneSurface([loop]))
        gmsh.model.geo.synchronize()
        distance=gmsh.model.mesh.field.add('Distance')
        gmsh.model.mesh.field.setNumbers(distance,'CurvesList',list(lines.values()))
        gmsh.model.mesh.field.setNumber(distance,'Sampling',1000)
        size=gmsh.model.mesh.field.add('MathEval')
        gmsh.model.mesh.field.setString(size,'F',f'Min({bulk_size}, {edge_size}+{growth}*F{distance})')
        fields=[size];facet_limits=[]
        if maximum_size_to_facet_width is not None:
            for surface,facet in zip(surfaces,polygons,strict=True):
                width=facet_width(facet['polygon'],facet['normal'])
                limit=maximum_size_to_facet_width*width
                constant=gmsh.model.mesh.field.add('MathEval')
                gmsh.model.mesh.field.setString(constant,'F',str(limit))
                restrict=gmsh.model.mesh.field.add('Restrict')
                gmsh.model.mesh.field.setNumber(restrict,'InField',constant)
                gmsh.model.mesh.field.setNumbers(restrict,'SurfacesList',[surface])
                gmsh.model.mesh.field.setNumber(restrict,'IncludeBoundary',1)
                fields.append(restrict)
                facet_limits.append(dict(surface=surface,minimum_width_m=width,maximum_size_target_m=limit))
        background=gmsh.model.mesh.field.add('Min')
        gmsh.model.mesh.field.setNumbers(background,'FieldsList',fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(background)
        gmsh.model.mesh.generate(2)
        tags,xyz,_=gmsh.model.mesh.getNodes();xyz=np.asarray(xyz).reshape(-1,3)
        lookup={int(tag):i for i,tag in enumerate(tags)}
        pieces=[]
        for surface,facet in zip(surfaces,polygons,strict=True):
            types,_,nodes=gmsh.model.mesh.getElements(2,surface)
            if list(types)!=[2]:
                raise ValueError('linear triangles required')
            ids=np.array([lookup[int(tag)] for tag in nodes[0]]).reshape(-1,3)
            t=xyz[ids]
            if np.mean(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])@facet['normal'])<0:
                t=t[:,[0,2,1]]
            pieces.append(t)
        half=np.concatenate(pieces)
        # Enforce the intended shared reflection boundary at coordinate roundoff.
        half[:,:,1][abs(half[:,:,1])<COORD_ATOL]=0
        reflected=half[:,[0,2,1]].copy();reflected[:,:,1]*=-1
        result=farthest_first(np.concatenate((half,reflected)))
        quality=verify_surface(original,result)
        return result,dict(gmsh_version=gmsh.__version__,bulk_size_m=bulk_size,
            edge_size_m=edge_size,growth=growth,geometry=quality,
            maximum_size_to_facet_width=maximum_size_to_facet_width,facet_size_limits=facet_limits)
    finally:
        gmsh.finalize()

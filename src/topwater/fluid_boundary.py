"""Conforming planar fluid boundary from a convex hull and explicit surface graph."""
from dataclasses import dataclass
import numpy as np
from scipy.spatial import ConvexHull
from .inputs import HullMesh
from .frames import rotation_matrix


def _cross(a,b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def _ears(indices,vertices,normal):
    """Planar polygon triangulation retaining every boundary vertex."""
    axis=np.eye(3)[np.argmin(abs(normal))]
    u=np.cross(axis,normal);u/=np.linalg.norm(u);v=np.cross(normal,u)
    p=vertices@np.stack((u,v),axis=1)
    loop=list(indices)
    if sum(_cross(p[loop[i]],p[loop[(i+1)%len(loop)]]) for i in range(len(loop)))<0:loop.reverse()
    first=loop.index(min(loop));loop=loop[first:]+loop[:first]
    scale=max(np.ptp(p[loop],axis=0));eps=64*np.finfo(float).eps*scale**2
    triangles=[]
    while len(loop)>3:
        for i in range(len(loop)):
            a,b,c=loop[i-1],loop[i],loop[(i+1)%len(loop)]
            if _cross(p[b]-p[a],p[c]-p[a])<=eps:continue
            other=[j for j in loop if j not in (a,b,c)]
            q=p[other]
            inside=(_cross(p[b]-p[a],q-p[a])>=-eps)&(_cross(p[c]-p[b],q-p[b])>=-eps)&(_cross(p[a]-p[c],q-p[c])>=-eps)
            if np.any(inside):continue
            triangles.append([a,b,c]);loop.pop(i);break
        else:raise ValueError('wet facet polygon cannot be triangulated without losing a boundary vertex')
    if _cross(p[loop[1]]-p[loop[0]],p[loop[2]]-p[loop[0]])<=eps:
        raise ValueError('degenerate final wet-facet triangle')
    return triangles+[loop]


@dataclass(frozen=True)
class FluidBoundary:
    mesh: HullMesh
    tags: np.ndarray
    surface_vertex_indices: np.ndarray
    diagnostics: dict
    body: object
    state: object


def closed_boundary(body,state,surface_vertices,surface_faces,*,contact,outer,center_xy,depth,tolerance=1e-10):
    """Close a represented graph with actual wet facets, cylinder wall and bottom.

    contact/outer list the two CCW boundary loops in the supplied surface mesh.
    Depth and artificial-boundary conditions are separate physical declarations;
    this geometry constructor does not imply a seabed or qualify truncation.
    """
    sv=np.asarray(surface_vertices,float);sf=np.asarray(surface_faces,np.int64)
    contact=np.asarray(contact,np.int64);outer=np.asarray(outer,np.int64);origin=np.asarray(center_xy,float)
    if (sv.ndim!=2 or sv.shape[1]!=3 or not np.isfinite(sv).all() or sf.ndim!=2 or sf.shape[1]!=3
            or np.any(sf<0) or np.any(sf>=len(sv)) or len(contact)<3 or len(outer)<3
            or origin.shape!=(2,) or not np.isfinite(depth) or depth<=0):
        raise ValueError('finite surface graph, loops, origin and positive depth required')
    full=(body.mesh.vertices-body.cg)@rotation_matrix(state.orientation).T+state.position
    if np.min(full[:,2])<=-depth+tolerance:raise ValueError('computational bottom intersects the hull')
    t=full[body.mesh.faces];normals=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    planes=np.c_[normals,np.einsum('ni,ni->n',normals,t[:,0])]
    if np.max(full@planes[:,:3].T-planes[:,3])>tolerance:raise ValueError('convex supporting hull required')
    for triangle in sv[sf]:
        polygon=list(triangle)
        for plane in planes:
            if not polygon:break
            clipped=[];previous=polygon[-1];sp=previous@plane[:3]-plane[3]+tolerance
            for current in polygon:
                sc=current@plane[:3]-plane[3]+tolerance
                if (sp<=0)!=(sc<=0):clipped.append(previous+(current-previous)*sp/(sp-sc))
                if sc<=0:clipped.append(current)
                previous,sp=current,sc
            polygon=clipped
        if len(polygon)>=3:
            p=np.array(polygon)
            overlap=abs(np.sum(_cross(p[:,:2],np.roll(p[:,:2],-1,axis=0))))/2
            if overlap>tolerance**2:raise ValueError('surface triangle intersects hull interior')
    physical=set(body.mesh.physical_faces.tolist());groups=[]
    for i,p in enumerate(planes):
        match=next((g for g in groups if np.max(abs(p-g[0]))<=tolerance),None)
        if match is None:groups.append([p,[i]])
        else:match[1].append(i)
    pool=[p.copy() for p in sv]
    center=np.r_[origin,np.mean(sv[contact,2])]
    if np.max(center@planes[:,:3].T-planes[:,3])>=-tolerance:
        raise ValueError('contact fan center must lie strictly inside the hull')
    cap=np.array([[center,sv[contact[i]],sv[contact[(i+1)%len(contact)]]] for i in range(len(contact))])
    graph=np.r_[sv[sf],cap];xy=graph[:,:,:2]
    edge1=xy[:,1]-xy[:,0];edge2=xy[:,2]-xy[:,0];det=_cross(edge1,edge2)
    if np.any(det<=0):raise ValueError('surface and contact loops must give an upward single-valued graph')
    def height(point):
        d=point[:2]-xy[:,0];a=_cross(d,edge2)/det;b=_cross(edge1,d)/det
        valid=np.flatnonzero((a>=-1e-10)&(b>=-1e-10)&(a+b<=1+1e-10))
        if not len(valid):raise ValueError('hull point lies outside represented surface domain')
        j=valid[0]
        return graph[j,0,2]+a[j]*(graph[j,1,2]-graph[j,0,2])+b[j]*(graph[j,2,2]-graph[j,0,2])
    if len(body.mesh.opening_rim) and any(p[2]<=height(p)+tolerance for p in full[body.mesh.opening_rim]):
        raise ValueError('surface reaches opening rim')
    contact_distance=sv[contact]@planes[:,:3].T-planes[:,3]
    if np.max(contact_distance)>tolerance or np.any(np.min(abs(contact_distance),axis=1)>tolerance):
        raise ValueError('contact must lie on the actual hull')
    def vertex(point):
        distance=np.max(abs(np.array(pool)-point),axis=1);found=np.flatnonzero(distance<=tolerance)
        if len(found):return int(found[0])
        pool.append(point.copy());return len(pool)-1
    wet=[]
    for p,group in groups:
        keep=[i in physical for i in group]
        if not any(keep):continue
        if not all(keep):raise ValueError('coplanar physical and artificial closure need separate facet identity')
        ids=np.unique(body.mesh.faces[group]);points=full[ids]
        axis=np.eye(3)[np.argmin(abs(p[:3]))];u=np.cross(axis,p[:3]);u/=np.linalg.norm(u);v=np.cross(p[:3],u)
        order=ConvexHull(points@np.stack((u,v),axis=1)).vertices
        corners=points[order];edges=set()
        on=abs(sv[contact]@p[:3]-p[3])<=tolerance
        for i,j in zip(contact,np.roll(contact,-1)):
            if abs(sv[i]@p[:3]-p[3])<=tolerance and abs(sv[j]@p[:3]-p[3])<=tolerance:
                edges.add(tuple(sorted((int(i),int(j)))))
        for a,b in zip(corners,np.roll(corners,-1,axis=0)):
            d=b-a;length2=d@d;lam=(sv[contact]-a)@d/length2
            near=np.linalg.norm(sv[contact]-(a+lam[:,None]*d),axis=1)<=tolerance
            hits=np.flatnonzero(on & near & (lam>=-tolerance)&(lam<=1+tolerance))
            split=[(0.,vertex(a)),(1.,vertex(b))]+[(float(lam[k]),int(contact[k])) for k in hits]
            split=sorted(set(split))
            for (_,i),(_,j) in zip(split[:-1],split[1:]):
                if i==j:continue
                middle=(pool[i]+pool[j])/2
                if middle[2]<=height(middle)+tolerance:edges.add(tuple(sorted((i,j))))
        if not edges:continue
        nodes=sorted({i for edge in edges for i in edge});coords=np.array(pool)
        if len(nodes)<3 or np.linalg.matrix_rank((coords[nodes]-coords[nodes[0]])@np.stack((u,v),axis=1),tol=tolerance)<2:continue
        adjacent={i:[] for i in nodes}
        for i,j in edges:adjacent[i].append(j);adjacent[j].append(i)
        if any(len(q)!=2 for q in adjacent.values()):raise ValueError('wet facet boundary is not a collection of simple loops')
        remaining=set(nodes)
        while remaining:
            loop=[min(remaining)];previous=None;current=loop[0]
            while True:
                nxt=next(j for j in adjacent[current] if j!=previous)
                if nxt==loop[0]:break
                if nxt in loop:raise ValueError('self-intersecting wet facet graph')
                loop.append(nxt);previous,current=current,nxt
            remaining.difference_update(loop)
            wet.extend(_ears(loop,coords,p[:3]))
    if not wet:raise ValueError('no wetted hull facets')
    # Hull facets above were body-outward; the fluid domain uses their reverse.
    faces=[*sf.tolist(),*[q[::-1] for q in wet]];tags=['surface']*len(sf)+['body']*len(wet)
    bottom=[]
    for i in outer:
        if sv[i,2]<=-depth:raise ValueError('surface intersects bottom')
        bottom.append(vertex(np.r_[sv[i,:2],-depth]))
    center_bottom=vertex(np.r_[origin,-depth])
    for i in range(len(outer)):
        j=(i+1)%len(outer);a,b=int(outer[i]),int(outer[j]);c,d=bottom[i],bottom[j]
        faces.extend([[a,c,d],[a,d,b],[center_bottom,d,c]]);tags.extend(['wall','wall','bottom'])
    faces=np.asarray(faces,np.int64);used,inverse=np.unique(faces,return_inverse=True)
    reindex=np.full(len(pool),-1,np.int64);reindex[used]=np.arange(len(used))
    mesh=HullMesh(np.array(pool)[used],inverse.reshape(-1,3))
    tags=np.array(tags);tags.setflags(write=False);surface_ids=reindex[np.arange(len(sv))];surface_ids.setflags(write=False)
    if np.any(surface_ids<0):raise ValueError('unused supplied surface vertex')
    return FluidBoundary(mesh,tags,surface_ids,dict(panel_count=len(faces),vertices=len(used),
        panels_by_tag={name:int(np.sum(tags==name)) for name in ('surface','body','wall','bottom')}),body,state)

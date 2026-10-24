"""External mean-fluid sampling geometry and explicitly checked surface pullback."""
import numpy as np
from scipy.optimize import brentq


def hull_planes(vertices,faces,tolerance=1e-10):
    """Convex closed hull only; retain actual oriented support planes."""
    v=np.asarray(vertices);t=v[np.asarray(faces)]
    normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    planes=np.c_[normal,np.einsum('ni,ni->n',normal,t[:,0])]
    if not np.isfinite(planes).all() or np.max(v@planes[:,:3].T-planes[:,3])>tolerance:
        raise ValueError('radial reference map requires verified convex supporting hull facets')
    unique=[]
    for p in planes:
        if not any(np.max(abs(p-q))<=tolerance for q in unique):unique.append(p)
    return np.array(unique)


def radius(planes,angles,heights,origin):
    a=np.atleast_1d(angles);z=np.broadcast_to(heights,a.shape)
    d=np.c_[np.cos(a),np.sin(a)]
    numerator=planes[:,3]-np.asarray(origin)@planes[:,:2].T-z[:,None]*planes[:,2]
    denominator=d@planes[:,:2].T
    if not np.isfinite(numerator).all() or np.any(numerator<-1e-10):
        raise ValueError('polar origin must be inside every horizontal hull section')
    bounds=np.full(denominator.shape,np.inf)
    np.divide(numerator,denominator,out=bounds,where=denominator>1e-14)
    r=bounds.min(axis=1)
    if not np.isfinite(r).all() or np.any(r<=0):raise ValueError('missing positive bounded hull section')
    return r


def mean_angles(wetted,origin,divisions):
    corners=np.unique(wetted.reshape(-1,3),axis=0)
    corners=corners[abs(corners[:,2])<=1e-11]
    if len(corners)<3:raise ValueError('missing frozen mean-waterline vertices')
    angles=np.mod(np.arctan2(corners[:,1]-origin[1],corners[:,0]-origin[0]),2*np.pi)
    a=np.sort(np.r_[angles,np.arange(divisions)*2*np.pi/divisions])
    return a[np.r_[True,np.diff(a)>1e-12]]


def annulus(angles,inner,outer,offsets,origin,heights=None):
    angles=np.asarray(angles);inner=np.asarray(inner);offsets=np.asarray(offsets)
    r=np.r_[inner[None,:]+offsets[:,None],np.full((1,len(angles)),outer)]
    if np.any(np.diff(r,axis=0)<=0):raise ValueError('radial rows must remain ordered inside outer boundary')
    xy=np.asarray(origin)+r[:,:,None]*np.c_[np.cos(angles),np.sin(angles)][None,:,:]
    z=np.zeros(r.shape) if heights is None else np.broadcast_to(heights,r.shape)
    vertices=np.concatenate((xy,z[:,:,None]),axis=2).reshape(-1,3)
    n=len(angles);j,i=np.meshgrid(np.arange(len(r)-1),np.arange(n),indexing='ij')
    a=j*n+i;b=(j+1)*n+i;c=(j+1)*n+(i+1)%n;d=j*n+(i+1)%n
    faces=np.stack((np.stack((a,b,c),axis=-1),np.stack((a,c,d),axis=-1)),axis=-2).reshape(-1,3)
    return vertices,faces,r


def exterior_points(points,planes,tolerance=1e-10):
    p=np.asarray(points)
    if (p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()
            or np.any(p[:,2]>0) or np.any(np.max(p@planes[:,:3].T-planes[:,3],axis=1)<=tolerance)):
        raise ValueError('reference evaluations require strictly exterior mean-fluid points at z<=0')


def periodic_interpolate(angles,values,target):
    a=np.asarray(angles);v=np.asarray(values);q=np.mod(target,2*np.pi)
    return np.array([np.interp(q,np.r_[a,a[0]+2*np.pi],np.r_[row,row[0]])
                     for row in v.reshape(-1,len(a))]).reshape(v.shape[:-1]+np.shape(q))


def _overlap_area(triangle,planes,tolerance):
    # Clip to a strictly inset convex hull. Positive area is forbidden, including
    # overlap missed by testing just vertices/centroids of a free-surface panel.
    polygon=list(triangle)
    for p in planes:
        if not polygon:return 0.
        clipped=[];previous=polygon[-1];sp=previous@p[:3]-p[3]+tolerance
        for current in polygon:
            sc=current@p[:3]-p[3]+tolerance
            if (sp<=0)!=(sc<=0):clipped.append(previous+(current-previous)*sp/(sp-sc))
            if sc<=0:clipped.append(current)
            previous,sp=current,sc
        polygon=clipped
    if len(polygon)<3:return 0.
    p=np.array(polygon)
    return abs(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))/2


def check_surface(vertices,faces,planes,contact_count,tolerance=1e-10):
    triangles=vertices[faces]
    twice=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])[:,2]
    if np.any(twice<=tolerance**2):raise ValueError('folded or degenerate surface pullback')
    contact=vertices[:contact_count];distance=contact@planes[:,:3].T-planes[:,3]
    if np.max(distance)>tolerance or np.any(np.min(abs(distance),axis=1)>tolerance):
        raise ValueError('contact vertices do not lie on physical hull')
    on=abs(distance)<=tolerance
    if np.any(~np.any(on & np.roll(on,-1,axis=0),axis=1)):
        raise ValueError('contact edge crosses hull facets without an inserted feature vertex')
    overlap=max(_overlap_area(t,planes,tolerance) for t in triangles)
    if overlap>tolerance**2:raise ValueError('free-surface triangle overlaps actual hull interior')
    return dict(minimum_projected_twice_area_m2=float(twice.min()),maximum_interior_overlap_m2=float(overlap),
                contact_vertices=contact_count,triangles=len(faces))


def mapped_surface(angles,eta,mean_planes,actual_planes,origin,outer,offsets):
    """Pull back elevation by radial coordinates; explicitly insert facet switches.

    eta is a represented field on mean annulus rows, including boundary limit.
    This construction does not qualify the approximation used to obtain eta.
    """
    angles=np.asarray(angles);eta=np.asarray(eta)
    if eta.shape!=(len(offsets)+1,len(angles)) or not np.isfinite(eta).all():
        raise ValueError('finite elevation on every mean radial row required')
    extra=[];closed=np.r_[angles,angles[0]+2*np.pi]
    def height(x):return float(periodic_interpolate(angles,eta[:1],np.array([x]))[0,0])
    def cross_bound(x,p,q):
        d=np.array([np.cos(x),np.sin(x)]);z=height(x)
        return ((p[3]-origin@p[:2]-z*p[2])*(d@q[:2])
                -(q[3]-origin@q[:2]-z*q[2])*(d@p[:2]))
    for lo,hi in zip(closed[:-1],closed[1:]):
        for i,p in enumerate(actual_planes):
            for q in actual_planes[i+1:]:
                left,right=cross_bound(lo,p,q),cross_bound(hi,p,q)
                if left*right>=0:continue
                x=brentq(cross_bound,lo,hi,args=(p,q),xtol=1e-14)
                d=np.array([np.cos(x),np.sin(x)])
                if min(d@p[:2],d@q[:2])<=1e-14:continue
                r=radius(actual_planes,[x],height(x),origin)[0]
                point=np.r_[origin+r*d,height(x)]
                if max(abs(point@p[:3]-p[3]),abs(point@q[:3]-q[3]))<1e-9:extra.append(x%(2*np.pi))
    target=np.sort(np.r_[angles,extra]);target=target[np.r_[True,np.diff(target)>1e-12]]
    heights=periodic_interpolate(angles,eta,target)
    mean=radius(mean_planes,target,0.,origin);actual=radius(actual_planes,target,heights[0],origin)
    vertices,faces,r=annulus(target,mean,outer,offsets,origin,heights)
    blend=(outer-r)/(outer-mean[None,:])
    mapped_r=r+blend*(actual-mean)[None,:]
    if np.any(np.diff(mapped_r,axis=0)<=0):raise ValueError('non-injective radial map')
    vertices[:,:2]=(origin+mapped_r[:,:,None]*np.c_[np.cos(target),np.sin(target)][None,:,:]).reshape(-1,2)
    checks=check_surface(vertices,faces,actual_planes,len(target))
    checks['inserted_angles']=len(target)-len(angles)
    return target,vertices,faces,checks

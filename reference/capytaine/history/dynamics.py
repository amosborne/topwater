"""Superseded pre-WEC-Sim prototype; preserved only, not the reference runtime."""
import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import lu_factor,lu_solve


def ramp(times,duration):
    u=np.clip(np.asarray(times)/duration,0.,1.)
    return u**3*(10.+u*(-15.+6.*u))


def radiation_kernel(omega,damping,times):
    """Finite-band cosine transform, preserving all couplings without fitting."""
    omega=np.asarray(omega);damping=np.asarray(damping);times=np.asarray(times)
    if damping.shape!=(len(omega),6,6) or not np.isfinite(damping).all() or np.any(np.diff(omega)<=0):raise ValueError('complete finite ordered radiation data required')
    K=np.empty((len(times),6,6))
    for start in range(0,len(times),32):
        cosine=np.cos(times[start:start+32,None]*omega)
        K[start:start+32]=2/np.pi*trapezoid(cosine[:,:,None,None]*damping[None],omega,axis=1)
    return K


def integrate(times,mass,added_infinite,restoring,kernel,force):
    """Newmark average acceleration with trapezoidal causal memory, from rest."""
    t=np.asarray(times);dt=t[1]-t[0];N=len(t)
    if t[0]!=0 or dt<=0 or not np.allclose(np.diff(t),dt,rtol=1e-10,atol=1e-12):raise ValueError('uniform times from zero required')
    mass=np.asarray(mass);added=np.asarray(added_infinite);C=np.asarray(restoring);K=np.asarray(kernel);F=np.asarray(force)
    if any(x.shape!=(6,6) for x in (mass,added,C)) or K.shape!=(N,6,6) or F.shape!=(N,6):raise ValueError('full six-mode arrays required')
    if not all(np.isfinite(x).all() for x in (mass,added,C,K,F)):raise ValueError('finite dynamics inputs required')
    q=np.zeros((N,6));v=q.copy();a=q.copy();effective=mass+added
    a[0]=np.linalg.solve(effective,F[0])
    lu=lu_factor(effective+.25*dt**2*C+.25*dt**2*K[0])
    for n in range(1,N):
        qp=q[n-1]+dt*v[n-1]+.25*dt**2*a[n-1]
        vp=v[n-1]+.5*dt*a[n-1]
        memory=dt*np.einsum('kij,kj->i',K[1:n][::-1],v[1:n])+.5*dt*K[0]@vp
        a[n]=lu_solve(lu,F[n]-C@qp-memory)
        q[n]=qp+.25*dt**2*a[n];v[n]=vp+.5*dt*a[n]
    if not all(np.isfinite(x).all() for x in (q,v,a)):raise ValueError('nonfinite motion history')
    return q,v,a

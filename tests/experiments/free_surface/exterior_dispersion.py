"""Spatial wave-dispersion falsifier for the actual P1 water exterior.

One six-tetrahedron cell is periodic horizontally with a Bloch phase. Interior
vertical potentials minimize the same assembled energy as the canoe exterior;
the free-surface mass is its consistent P1 mass. No hull, absorber or time
integrator participates. Five-period phase is predicted from the semidiscrete
frequency, not measured from a nonlinear trajectory.
"""
from itertools import permutations
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix

from .pocket import stiffness
from .surface_damping import matrices


def symbol(dx,dy,z,kx,ky):
    """Discrete DtN eigenvalue using current tetra energy and surface mass."""
    z=np.asarray(z,dtype=float)
    if (z.ndim!=1 or len(z)<2 or not np.all(np.isfinite(z)) or np.any(np.diff(z)<=0)
            or not np.all(np.isfinite([dx,dy,kx,ky])) or dx<=0 or dy<=0):
        raise ValueError('positive cell widths, increasing finite layers and finite wavenumbers required')
    layers=len(z)
    points=np.array([[i*dx,j*dy,value] for i in range(2) for j in range(2) for value in z])
    numbering=np.arange(len(points)).reshape(2,2,layers)
    paths=np.zeros((6,4,3),dtype=int)
    paths[:,1:]=np.cumsum(np.eye(3,dtype=int)[list(permutations(range(3)))],axis=1)
    tets=np.array([numbering[tuple((path+[0,0,level]).T)]
        for level in range(layers-1) for path in paths])
    phase=np.exp(1j*(kx*points[:,0]+ky*points[:,1]))
    levels=np.tile(np.arange(layers),4)
    periodic=coo_matrix((phase,(np.arange(len(points)),levels)),
        shape=(len(points),layers)).tocsr()
    reduced=(periodic.conj().T@stiffness(points,tets)@periodic).toarray()
    profile=np.r_[np.linalg.solve(reduced[:-1,:-1],-reduced[:-1,-1]),1.]
    energy=profile.conj()@reduced@profile
    top=numbering[:,:,-1]
    faces=np.array([[top[0,0],top[1,0],top[1,1]],[top[0,0],top[1,1],top[0,1]]])
    free=top.ravel()
    mass,_=matrices(points,faces,free,np.zeros(4))
    surface_mass=phase[free].conj()@mass@phase[free]
    if abs(energy.imag)>1e-11*max(1.,abs(energy)) or abs(surface_mass.imag)>1e-11*abs(surface_mass):
        raise ArithmeticError('periodic fluid energy or surface mass is not real')
    return float(energy.real/surface_mass.real)


def uniform_symbol(dx,dy,nz,kx,ky,depth=4.):
    """Independent closed form for this six-tetrahedron uniform-layer stencil."""
    ax,ay=kx*dx,ky*dy
    horizontal=4*np.sin(ax/2)**2/dx**2+4*np.sin(ay/2)**2/dy**2
    dz=depth/nz
    gamma=2*np.arcsinh(.5*dz*np.sqrt(horizontal))
    mass=.5+(np.cos(ax)+np.cos(ay)+np.cos(ax+ay))/6
    return float(np.sinh(gamma)*np.tanh(nz*gamma)/(dz*mass))


def measure(n,omega,angle,*,horizontal_refinement=1,vertical_refinement=1):
    """Compare semidiscrete frequency with finite-depth gravity-wave dispersion."""
    dx,dy=4*np.pi/n/horizontal_refinement,2*np.pi/n/horizontal_refinement
    nz=n//3*vertical_refinement
    k=brentq(lambda value:9.81*value*np.tanh(4*value)-omega**2,1e-10,omega**2/9.81+10)
    kx,ky=k*np.cos(angle),k*np.sin(angle)
    actual=symbol(dx,dy,np.linspace(-4.,0.,nz+1),kx,ky)
    independent=uniform_symbol(dx,dy,nz,kx,ky)
    if not np.isclose(actual,independent,rtol=2e-11,atol=2e-12):
        raise ArithmeticError('assembled exterior differs from its independently derived stencil')
    frequency=np.sqrt(9.81*actual)
    phase=10*np.pi*(frequency/omega-1)
    return dict(n=n,omega=omega,angle=angle,dx=dx,dy=dy,dz=4/nz,
        horizontal_refinement=horizontal_refinement,vertical_refinement=vertical_refinement,
        discrete_dtn=actual,exact_dtn=k*np.tanh(4*k),
        relative_frequency_error=float(frequency/omega-1),
        predicted_five_period_phase_error_radians=float(phase),
        phase_gate_passed=bool(abs(phase)<=.01))


if __name__=='__main__':
    cases=[]
    for n in (24,36,48):
        for omega in (2.,4.,6.):
            for angle in (0.,np.pi/2,np.pi/4,-np.pi/4):
                cases.append(measure(n,omega,angle))
    report=dict(qualified=False,phase_gate_radians=.01,
        phase_gate_passed=all(case['phase_gate_passed'] for case in cases),cases=cases,
        scope='Actual P1 exterior spatial dispersion only. Phase predicted after five exact '
            'continuum periods; no time-integration error. Does not qualify nonlinear wave, '
            'interface or canoe loads. A failed gate is retained, not corrected or fitted.')
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/exterior-dispersion.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(phase_gate_passed=report['phase_gate_passed'],
        maximum_phase_error_radians=max(abs(case['predicted_five_period_phase_error_radians']) for case in cases))))
    raise SystemExit(0 if report['phase_gate_passed'] else 1)

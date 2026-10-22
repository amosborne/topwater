"""Fresh-operator nonlinear near/far fluid baseline; body coupling is unfinished."""
from dataclasses import dataclass
from time import perf_counter
import numpy as np
from scipy.linalg import cho_solve
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from topwater import HullMesh
from topwater.rankine import TrianglePanels
from .evolution import Tank,make_tank as make_surface,bernoulli_load
from .operator import assemble_energy
from .pocket import grid,boundary,stiffness


@dataclass(frozen=True)
class PocketTank:
    surface: Tank
    points: np.ndarray
    exterior_tets: np.ndarray
    surface_nodes: np.ndarray
    local_nodes: np.ndarray
    local_faces: np.ndarray
    column_index: np.ndarray
    gpu: bool


def make_tank(n=12,shift=0,gpu=True):
    points,tets,near,fs=grid(n,0.,shift)
    # The pocket ends one metre below the reference surface; exterior FEM
    # remains below it as well as alongside it. All interface nodes are shared.
    near &= np.all(points[tets][:,:,2]>=-1.,axis=1)
    if not np.any(near):
        raise ValueError('the shallow pocket needs at least four vertical cells')
    faces=boundary(points,tets[near])
    nodes,inverse=np.unique(faces,return_inverse=True)
    local=HullMesh(points[nodes],inverse.reshape(-1,3))
    surface=make_surface((n,n,n//3))
    if not np.array_equal(points[fs],surface.mesh.vertices[surface.surface_nodes]):
        raise ValueError('surface and volume grid nodes do not match')
    return PocketTank(surface,points,tets[~near],fs,nodes,local.faces,
        np.repeat(np.arange(len(fs)),n//3+1),gpu)


def rhs(data,eta,psi,*,timing=None):
    start=perf_counter()
    if not np.all(np.isfinite(eta)) or not np.all(np.isfinite(psi)) or np.min(eta)<=-4:
        raise ValueError('nonfinite state or graph touching the tank bottom')
    points=data.points.copy()
    points[:,2]+=(points[:,2]+4)/4*eta[data.column_index]
    panels=TrianglePanels(points[data.local_nodes][data.local_faces],1)
    geometry_seconds=perf_counter()-start
    stage=perf_counter()
    if data.gpu:
        from .cuda_operator import assemble_energy as assemble
    else:
        assemble=assemble_energy
    local,diagnostics=assemble(data.local_faces,panels,3)
    local_seconds=perf_counter()-stage
    stage=perf_counter()
    energy=stiffness(points,data.exterior_tets)
    nodes=data.local_nodes
    energy+=coo_matrix((local.ravel(),(np.repeat(nodes,len(nodes)),np.tile(nodes,len(nodes)))),
        shape=energy.shape).tocsr()
    exterior_seconds=perf_counter()-stage
    active=np.flatnonzero(energy.diagonal()>0)
    fs=data.surface_nodes
    unknown=np.setdiff1d(active,fs)
    potential=np.zeros(len(points));potential[fs]=psi
    stage=perf_counter()
    operator=energy[unknown][:,unknown].tocsc()
    right=-energy[unknown][:,fs]@psi
    factor=splu(operator)
    potential[unknown]=factor.solve(right)
    flux=energy[fs]@potential
    solve_seconds=perf_counter()-stage
    surface=data.surface
    eta_t=cho_solve(surface.projected_factor,flux)
    nonlinear=bernoulli_load(surface,eta,psi,eta_t)
    psi_t=-9.81*eta+cho_solve(surface.projected_factor,nonlinear)
    if not np.all(np.isfinite(eta_t)) or not np.all(np.isfinite(psi_t)):
        raise FloatingPointError('nonfinite near/far fluid stage')
    mass=surface.projected_mass
    kinetic=.5*psi@flux
    total=kinetic+.5*9.81*eta@mass@eta
    if timing is not None:
        timing.update(geometry_seconds=geometry_seconds,local_seconds=local_seconds,
            exterior_seconds=exterior_seconds,solve_seconds=solve_seconds,
            total_seconds=perf_counter()-start,local_panels=len(panels),local_nodes=len(nodes),
            total_energy=float(total),
            global_unknowns=len(unknown),matrix_bytes=energy.data.nbytes+energy.indices.nbytes+energy.indptr.nbytes,
            original_algebraic_absolute_residual=float(np.max(np.abs(operator@potential[unknown]-right))),
            local_diagnostics=diagnostics)
    return eta_t,psi_t,float(np.ones(len(fs))@mass@eta_t),float(total)

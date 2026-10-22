"""Actual canonical-canoe instantaneous reaction through a shallow BEM pocket.

The surrounding finite-domain water is sparse tetrahedral FEM. Free-surface
potential is zero (infinite-frequency limit); this is not finite-frequency wave
correlation, nonlinear evolution or a production backend.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits
from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from topwater.rankine import TrianglePanels
from tests.reference_hull import input_data
from tests.experiments.free_surface.operator import assemble_energy
from tests.experiments.free_surface.canoe import body_load,join_wetted_hull
from tests.benchmark_rankine import peak_rss_bytes
from tests.experiments.capytaine_mesh_cost import subdivide
from .pocket import grid,boundary,stiffness


def run(n,subdivisions,expected):
    start=perf_counter()
    data=SimulationInputs.from_json(json.dumps(input_data(2)))
    rest,wet=equilibrate_heave(data.body,data.initial_state)
    hull=wet.wetted_triangles
    for _ in range(subdivisions):
        hull=subdivide(hull)
    result=run_hull(n,hull,rest.position,expected,start=start)
    result['subdivisions']=subdivisions
    return result


def run_hull(n,hull,center,expected,*,start=None):
    """Solve the same pocket for supplied physical wet triangles and body CG.

    The caller supplies body-outward hull winding, as in hydrostatics. This remains
    the canonical convex, horizontal-waterline experiment, not a mesh API.
    """
    if start is None:
        start=perf_counter()
    points,tets,near,fs=grid(n,0.,0)
    # An actual local pocket may end below the hull, rather than at the seabed.
    near &= np.all(points[tets][:,:,2]>=-1.,axis=1)
    faces=boundary(points,tets[near])
    top=np.all(points[faces,2]==0,axis=1)
    walls=points[faces[~top]]
    edges=np.stack((hull,np.roll(hull,-1,axis=1)),axis=2).reshape(-1,2,3)
    segments=edges[np.all(edges[:,:,2]==0,axis=1)]
    mesh,body,free=join_wetted_hull(hull,segments,walls,points[np.unique(faces[top])])
    panels=TrianglePanels.from_mesh(mesh)
    mapping=[];lookup={tuple(point):i for i,point in enumerate(points)}
    additional=[]
    for vertex in mesh.vertices:
        key=tuple(vertex)
        if key not in lookup:
            lookup[key]=len(points)+len(additional);additional.append(vertex)
        mapping.append(lookup[key])
    mapping=np.asarray(mapping)
    points=np.vstack((points,additional))
    free_nodes=np.unique(np.r_[fs,mapping[np.unique(mesh.faces[free])]])
    geometry_seconds=perf_counter()-start
    assembly_start=perf_counter()
    energy=stiffness(points,tets[~near])
    local,diagnostics=assemble_energy(mesh.faces,panels,3)
    energy+=coo_matrix((local.ravel(),(np.repeat(mapping,len(mapping)),
        np.tile(mapping,len(mapping)))),shape=energy.shape).tocsr()
    load=np.zeros((len(points),6));load[mapping]=body_load(mesh,panels,body,center)
    active=np.flatnonzero(energy.diagonal()>0)
    unknown=np.setdiff1d(active,free_nodes)
    operator=energy[unknown][:,unknown].tocsc()
    assembly_seconds=perf_counter()-assembly_start
    factor_start=perf_counter();factor=splu(operator)
    factor_seconds=perf_counter()-factor_start
    times=[]
    for _ in range(6):
        solve_start=perf_counter();answer=factor.solve(load[unknown]);times.append(perf_counter()-solve_start)
    potential=np.zeros_like(load);potential[unknown]=answer
    actual=1000*load.T@potential
    candidate_seconds=perf_counter()-start
    scale=np.array([1.,1.,1.,3.5,3.5,3.5])
    normalized=(actual-expected)/scale[:,None]/scale[None,:]/(1000*3.5**3)
    factor_bytes=sum(m.data.nbytes+m.indices.nbytes+m.indptr.nbytes for m in (factor.L,factor.U))
    return dict(qualified=False,n=n,pocket_bottom=-1.,hull_panels=len(body),
        local_panels=len(panels),local_nodes=len(mapping),unknowns=len(unknown),
        geometry_seconds=geometry_seconds,assembly_seconds=assembly_seconds,
        factor_seconds=factor_seconds,solve_first_seconds=times[0],solve_warm_seconds=times[1:],
        complete_candidate_seconds=candidate_seconds,
        matrix_bytes=energy.data.nbytes+energy.indices.nbytes+energy.indptr.nbytes,
        factor_bytes=factor_bytes,peak_rss_bytes=peak_rss_bytes(),
        actual=actual.tolist(),reference=expected.tolist(),
        diagonal_relative_error=(np.diag(actual-expected)/np.diag(expected)).tolist(),
        maximum_normalized_matrix_error=float(np.max(np.abs(normalized))),
        original_algebraic_relative_residual=float(np.linalg.norm(operator@answer-load[unknown])/np.linalg.norm(load[unknown])),
        **diagnostics,scope='Finite box/bed versus infinite-depth Capytaine; instantaneous acceleration reaction only, not a wave-body qualification.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    args=parser.parse_args()
    from tests.capytaine_reference import instantaneous_inertia
    with threadpool_limits(1):
        data=SimulationInputs.from_json(json.dumps(input_data(2)))
        rest,wet=equilibrate_heave(data.body,data.initial_state)
        reference=subdivide(subdivide(wet.wetted_triangles))
        expected,damping=instantaneous_inertia(reference,rest.position,data.body.mass)
        if np.any(damping):
            raise RuntimeError('infinite-frequency reference has nonzero damping')
        rows=[]
        for n,subdivisions in ((12,0),(24,0),(12,1),(12,2)):
            row=run(n,subdivisions,expected);rows.append(row)
            print(json.dumps(row),flush=True)
    report=dict(qualified=False,capytaine_reference_panels=len(reference),results=rows,
        scope='Finite-domain instantaneous reaction comparison, not finite-frequency or nonlinear wave-body qualification.')
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-pocket-canoe.json').write_text(json.dumps(report,indent=2)+'\n')

"""Material-first P1 mesh comparison on the canonical canoe; not a solver mode.

All layouts retain the same planar physical hull. A smaller mesh must pass the
same six-mode check at held-out poses; a resting-pose improvement is insufficient.
"""
from dataclasses import replace
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits
from topwater import SimulationInputs
from topwater._clipping import clip_flat_water
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave, hydrostatics
from tests.reference_hull import input_data
from tests.experiments.capytaine_mesh_cost import subdivide
from .pocket_canoe import run_hull


def material_mesh(body,transverse=False):
    """Fixed numerical subdivision of the canonical hull's physical faces.

    The optional bisection adds cross-sectional potential resolution. Shared
    edges make the same geometric decision on both incident faces. This uses
    the reference hull's longitudinal x axis, not a general adaptive mesher.
    """
    triangles=body.mesh.vertices[body.mesh.faces[body.mesh.physical_faces]]
    triangles=subdivide(subdivide(triangles))
    if transverse:
        children=[]
        for triangle in triangles:
            selected=triangle[:,0]==np.roll(triangle[:,0],-1)
            count=int(selected.sum())
            if count==0:
                children.append(triangle)
            elif count==1:
                edge=int(np.flatnonzero(selected)[0])
                a,b,c=np.roll(triangle,-edge,axis=0)
                midpoint=(a+b)/2
                children.extend(([a,midpoint,c],[midpoint,b,c]))
            elif count==3:
                children.extend(subdivide(triangle[None]))
            else:
                raise AssertionError('two equal-x edges imply three')
        triangles=np.asarray(children)
    vertices,inverse=np.unique(triangles.reshape(-1,3),axis=0,return_inverse=True)
    return vertices,inverse.reshape(-1,3)


def run():
    from tests.capytaine_reference import instantaneous_inertia

    data=SimulationInputs.from_json(json.dumps(input_data(2)))
    rest,_=equilibrate_heave(data.body,data.initial_state)
    materials={};preparation={}
    for name,transverse in (('material',False),('material_transverse',True)):
        start=perf_counter()
        materials[name]=material_mesh(data.body,transverse)
        preparation[name]=perf_counter()-start
    # Existing pocket engineering screen, not a physical accuracy requirement.
    screen=.034
    rows=[]
    for heave,roll in ((0.,0.),(.002,.01),(-.003,-.015)):
        state=replace(rest,position=rest.position+np.array([0,0,heave]),
            orientation=np.roll(Rotation.from_rotvec([roll,0.,0.]).as_quat(),1))
        wet=hydrostatics(data.body,state)
        expected,damping=instantaneous_inertia(subdivide(subdivide(wet.wetted_triangles)),
            state.position,data.body.mass)
        if np.any(damping):
            raise RuntimeError('infinite-frequency reference has nonzero damping')
        for name in ('postclip','material','material_transverse'):
            start=perf_counter()
            if name=='postclip':
                current=hydrostatics(data.body,state)
                hull=subdivide(subdivide(current.wetted_triangles))
            else:
                vertices,faces=materials[name]
                world=(vertices-data.body.cg)@rotation_matrix(state.orientation).T+state.position
                hull,*_=clip_flat_water(world,faces,np.arange(len(faces)),0.,True)
            result=run_hull(12,hull,state.position,expected,start=start)
            passed=bool(np.max(np.abs(result['diagonal_relative_error']))<screen)
            if result['original_algebraic_relative_residual']>1e-10:
                raise RuntimeError('current pocket equations were not solved accurately')
            if name=='postclip' and not passed:
                raise RuntimeError('reference layout failed the existing pocket screen')
            row=dict(layout=name,heave_m=heave,roll_rad=roll,
                engineering_screen_passed=passed,**result)
            rows.append(row)
            print(json.dumps({key:row[key] for key in ('layout','heave_m','roll_rad',
                'hull_panels','local_nodes','diagonal_relative_error','engineering_screen_passed')}),flush=True)
    return dict(qualified=False,diagonal_engineering_screen=screen,
        layout_screens={name:all(row['engineering_screen_passed'] for row in rows if row['layout']==name)
            for name in ('postclip','material','material_transverse')},
        material_preparation_seconds=preparation,results=rows,
        scope='Same physical canoe at rest and two prescribed heave/heel poses; fresh Capytaine per pose. '
            'Per-case time includes clipping and the complete pocket reaction solve, excludes one-time '
            'material preparation (reported separately) and Capytaine. First timings are not warm benchmarks. '
            'No pressure, evolved motion, finite-frequency or real-time qualification.')


if __name__=='__main__':
    with threadpool_limits(1):
        report=run()
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-material-layout.json').write_text(json.dumps(report,indent=2)+'\n')

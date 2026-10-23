"""Independent finite-difference incident-load derivative for momentum audits."""
import numpy as np
from .boat_world_field import boundary_data


def incident_load_rate(base,plus,minus,values,twist,time,step):
    """Cheap h-stencil of the SAME prescribed incident Neumann quadrature."""
    if not np.isfinite(step) or step<=0:
        raise ValueError('positive finite incident derivative step required')
    basis=('faces','body_faces','free_faces','free_nodes','body_source_group',
           'free_source_parent','free_parent_indices','node_chart_key')
    if 'node_chart_key' not in base:
        raise ValueError('explicit current chart ownership required')
    for moved in (plus,minus):
        if any(not np.array_equal(base[name],moved[name]) for name in basis):
            raise ValueError('incident derivative must preserve current chart/basis')
    loads=[]
    for moved,offset in ((plus,step),(minus,-step)):
        boundary=boundary_data(moved,values,twist,time=time+offset)
        loads.append(boundary['right']-boundary['body_load']@twist)
    return (loads[0]-loads[1])/(2*step)



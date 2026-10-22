"""Regress actual moving-canoe image interactions that stopped adaptive probes."""

from dataclasses import replace
import json

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation

from topwater import SimulationInputs
from topwater._clipping import clip_flat_water
from topwater.bem import assemble_dirichlet_halfspace
from topwater.frames import rotation_matrix
from topwater.hydrostatics import equilibrate_heave
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence, QuadratureError
from .experiments.capytaine_mesh_cost import subdivide
from .rankine_reference import edge_aligned_values
from .reference_hull import input_data


@pytest.fixture(scope='module')
def material_canoe():
    data = SimulationInputs.from_json(json.dumps(input_data(2)))
    rest, _ = equilibrate_heave(data.body, data.initial_state)
    triangles = data.body.mesh.vertices[data.body.mesh.faces[data.body.mesh.physical_faces]]
    # Refine coplanar material faces before moving/clipping; physical hull fixed.
    triangles = subdivide(subdivide(triangles))
    vertices, inverse = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    return data.body, rest, vertices, inverse.reshape(-1, 3)


@pytest.mark.parametrize('frame,index', [(1, 239), (48, 57)])
@pytest.mark.parametrize('sign', [-1, 1])
def test_moving_canoe_thin_waterline_images(material_canoe, frame, index, sign):
    body, rest, vertices, faces = material_canoe
    time = .02*frame
    state = replace(rest, position=rest.position+np.array([0, 0, .02*np.sin(1.6*time)]),
        orientation=np.roll(Rotation.from_rotvec([.05*np.sin(2*time+.4), 0, 0]).as_quat(), 1))
    world = (vertices-body.cg)@rotation_matrix(state.orientation).T+state.position
    wet, *_ = clip_flat_water(world, faces, np.arange(len(faces)), 0., True)
    panels = TrianglePanels(wet, sign)
    # These two actual replay frames previously failed complete assembly.
    operators = assemble_dirichlet_halfspace(panels)
    target = panels.centroids[index].copy(); target[2] *= -1
    image = panel_influence(target[None], TrianglePanels(wet[index:index+1], sign))
    expected = edge_aligned_values(target, wet[index], sign)
    assert_allclose([image.single_layer[0, 0], image.double_layer[0, 0]],
                    expected, rtol=1e-10, atol=1e-14)
    # Direct planar self K is zero; the entire diagonal is the reflected image.
    assert_allclose(operators.double_layer[index, index], -expected[1], rtol=1e-10)
    with pytest.raises(QuadratureError, match='max_refinement'):
        panel_influence(target[None], TrianglePanels(wet[index:index+1], sign),
                        controls=LaplaceControls(max_refinement=0))

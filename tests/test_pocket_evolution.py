"""Physical invariants and independent wave terms of the integrated fluid stage."""
import numpy as np
from threadpoolctl import threadpool_limits
from .experiments.free_surface.pocket_evolution import make_tank,rhs
from .experiments.free_surface.evolution import check_quadratic
from .experiments.free_surface.shape_check import check as check_shape


def test_near_far_calm_water_gauge_and_volume():
    with threadpool_limits(1):
        data=make_tank(12,gpu=False)
        size=len(data.surface_nodes)
        still=rhs(data,np.zeros(size),np.ones(size))
        assert np.max(np.abs(still[0]))<1e-11
        assert np.max(np.abs(still[1]))<1e-11
        x=data.points[data.surface_nodes,0]+2*np.pi
        moving=rhs(data,.02*np.cos(.25*x),.01*np.cos(.25*x))
        assert abs(moving[2])<1e-11
        assert moving[3]>0


def test_near_far_nonlinear_wave_terms_refine():
    results=[]
    with threadpool_limits(1):
        for n in (12,24):
            data=make_tank(n,gpu=False)
            result=check_quadratic(data.surface,
                lambda surface,eta,psi:rhs(data,eta,psi),.02,.25)
            assert max(abs(value) for value in result['integrated_elevation_rates'])<1e-11
            results.append(result)
    # Independent second-order finite-depth terms, integrated over the complete
    # free surface. This proves refinement, not that physical accuracy is enough.
    for metric in ('continuous_quadratic_kinematic_relative_rms',
                   'continuous_quadratic_dynamic_relative_rms'):
        assert results[1][metric]<results[0][metric]/2


def test_near_far_surface_work_refines():
    errors=[]
    with threadpool_limits(1):
        for n in (12,24):
            data=make_tank(n,gpu=False)
            result=check_shape(data.surface,lambda eta,psi:rhs(data,eta,psi),.25)
            coarse,fine=result['results']
            # A directional derivative of actual fluid energy independently
            # checks the dynamic surface work; the residual is spatial error.
            assert abs(coarse['finite_difference']/fine['finite_difference']-1)<1e-6
            errors.append(fine['relative_discrepancy'])
    assert errors[1]<errors[0]/2

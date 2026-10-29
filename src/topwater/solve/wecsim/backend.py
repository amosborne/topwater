"""Optional WEC-Sim motion operations used by the shared solve workflow."""
from pathlib import Path
import numpy as np
from topwater.solve.settings import DEFAULTS
from topwater.solve.types import PreparedMotion
from topwater.solve.fourier import surface_radiation_frequencies
from . import engine, exchange


def preflight(options, influence):
    if options.time_step_s != DEFAULTS.time_step_s:
        raise ValueError('Custom motion timestep requires the native solver')
    if influence is not None:
        raise ValueError('External influence integration requires the native solver')
    return engine.preflight()


def prepare(setup, workspace, options, runtime):
    h = setup.hydrostatic
    excitation, radiation, headings = exchange.reference_coefficient_grids(setup.train.omega[-1], options)
    hydro = exchange.build_hydrodynamics(setup.wet_triangles_m, setup.mean_cg_m,
        setup.mass_matrix, workspace/'hydrodynamics', hydrostatic_stiffness=h.stiffness,
        cb=h.center_of_buoyancy, volume=h.displaced_volume,
        excitation_omega=excitation, radiation_omega=radiation, headings=headings,
        field_omega=np.unique(np.r_[setup.train.omega[1:],
            surface_radiation_frequencies(setup.time_s, options.radiation_settings(1))[1:]]),
        green_cache_dir=Path.home()/'.cache/topwater/capytaine-3.0.0')
    return PreparedMotion(hydro)


def integrate(setup, prepared, workspace, options, influence, runtime):
    return engine.integrate_motion(setup, prepared.hydrodynamics, workspace/'body', runtime)

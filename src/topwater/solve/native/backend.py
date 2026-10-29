"""Native motion operations used by the shared solve workflow."""
from time import perf_counter
import sys
from topwater.solve.hydrodynamics import prepare_hydrodynamics
from topwater.solve.types import PreparedMotion
from .dynamics import integrate_motion
from .excitation import Excitation


def preflight(options, influence):
    """The native motion backend needs no external runtime."""
    return None


def prepare(setup, workspace, options, runtime):
    hydro = prepare_hydrodynamics(setup, workspace/'hydrodynamics', options)
    return PreparedMotion(hydro, Excitation.from_case(setup, hydro))


def integrate(setup, prepared, workspace, options, influence, runtime):
    started = perf_counter()
    print('STAGE native Python integration started', file=sys.stderr)
    motion = integrate_motion(setup, prepared.hydrodynamics, influence=influence,
        dt_s=options.time_step_s, excitation=prepared.inputs)
    print(f'STAGE native Python integration complete elapsed_s={perf_counter()-started:.1f}', file=sys.stderr)
    return motion

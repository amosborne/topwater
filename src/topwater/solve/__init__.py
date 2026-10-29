"""Shared solve workflow with independently selected motion and field backends."""
from .settings import SolveOptions, DEFAULT_DWELL_S
from .native.influence import BoatState, Influence, LoadEvaluation


def simulate(boat, sea, dwell_s=DEFAULT_DWELL_S, *, rest_s=10., ramp_s=5.,
             solver='native', influence=None, device='auto', options=None):
    """Return motion and total water geometry; external influences remain caller-owned."""
    from .workflow import solve
    return solve(boat, sea, dwell_s, rest_s=rest_s, ramp_s=ramp_s,
        solver=solver, influence=influence, device=device, options=options)


__all__ = ['BoatState', 'Influence', 'LoadEvaluation', 'SolveOptions', 'simulate']

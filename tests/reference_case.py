"""Shared boat/observed-sea inputs; not an initialized fluid state or a solver.

The boat starts at its calm-water equilibrium with bow east. A moving-wave
solver must initialize/ramp its fluid and contact state consistently. Density
and gravity are owned by reference_boat, and the incident model is deep water.
"""
from topwater import SimulationInputs
from topwater.sea_states import realize
from .reference_boat import BEAM, GRAVITY, equilibrium, reference_body


BASELINE_PHASE_SEED = 42


def reference_inputs(*, seed=BASELINE_PHASE_SEED, refinement=0):
    """Fixed baseline by default; seed=None requests a fresh recorded realization."""
    body = reference_body(refinement)
    rest, _ = equilibrium(refinement)
    waves = realize(beam_m=BEAM, gravity=GRAVITY, seed=seed)
    return SimulationInputs(body, waves, rest)

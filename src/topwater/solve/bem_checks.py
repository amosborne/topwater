"""Capytaine 3.0 frequency-bank checks, owned by solve preparation."""
import sys

import capytaine as cpt
import numpy as np
from capytaine.bem.problems_and_results import LinearPotentialFlowProblem
from capytaine.bem.problems_checks import _check_wavelength_and_mesh_resolution, _check_wavelength_and_irregular_frequencies, _check_wavelength_and_water_depth


def check_frequency_bank(body, frequencies, *, density, gravity, depth):
    """Check the complete bank before tabulation, BEM solves, or cache lookup."""
    omega = np.asarray(frequencies, dtype=float)
    if omega.ndim != 1 or not len(omega) or np.any(np.isnan(omega)) or np.any(omega <= 0):
        raise ValueError('BEM preflight requires positive frequencies, optionally including infinity')
    omega = np.unique(omega)
    problems = [LinearPotentialFlowProblem(body=body, omega=float(w), rho=density,
                g=gravity, water_depth=depth) for w in omega]
    print(f'BEM_PREFLIGHT frequencies={len(omega)} mesh_resolution=once irregular_frequencies=once',
          file=sys.stderr)
    _check_wavelength_and_mesh_resolution(problems)
    _check_wavelength_and_irregular_frequencies(problems)
    if depth != np.inf:
        _check_wavelength_and_water_depth(problems)
    return frozenset(float(w) for w in omega)


class PrecheckedBEMSolver(cpt.BEMSolver):
    """Retain ordinary numerical checks; skip the completed wavelength preflight."""
    def __init__(self, *, checked_frequencies, **kwargs):
        super().__init__(**kwargs)
        self.checked_frequencies = checked_frequencies

    def solve(self, problem, method=None, keep_details=True, n_threads=None, _check_wavelength=True):
        if float(problem.omega) not in self.checked_frequencies:
            raise ValueError('BEM frequency was not declared in the initial preflight bank; '
                             'include it in build_hydrodynamics(field_omega=...)')
        return super().solve(problem, method=method, keep_details=keep_details,
                             n_threads=n_threads, _check_wavelength=False)

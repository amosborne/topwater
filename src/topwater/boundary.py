"""Explicit reusable Laplace boundary-system factorization.

The owner holds operators, boundary assignments and length scaling together.
Changing any of them requires constructing a new system, not a cache lookup.
"""

from dataclasses import dataclass
from time import perf_counter
import warnings

import numpy as np
from numba import njit
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from scipy.linalg.lapack import dgecon, dlange
from threadpoolctl import ThreadpoolController

from ._validation import array
from .bem import LaplaceOperators

# NumPy/SciPy BLAS libraries are loaded above. Discover their controls once;
# rescanning all process libraries per factorization dominated small systems.
_blas = ThreadpoolController()


@njit(cache=True, nogil=True)
def _fill_mixed_matrix(single, double, dirichlet, length, matrix):
    """Tile C-order operators into the one Fortran-order factor buffer."""
    n = len(single)
    for first in range(0, n, 16):
        for start in range(0, n, 16):
            for col in range(start, min(start+16, n)):
                if dirichlet[col]:
                    for row in range(first, min(first+16, n)):
                        matrix[row, col] = -single[row, col]/length
                else:
                    for row in range(first, min(first+16, n)):
                        matrix[row, col] = double[row, col]
    for j in range(n):
        if not dirichlet[j]:
            matrix[j, j] += .5


@dataclass(frozen=True, eq=False)
class BoundarySolution:
    potential: np.ndarray
    normal_velocity: np.ndarray
    diagnostics: dict


@dataclass(frozen=True, eq=False, init=False)
class BoundarySystem:
    """Factor H*phi=S*q for one geometry and Dirichlet/Neumann assignment.

    True entries in dirichlet prescribe potential; False entries prescribe
    normal velocity. Supply panel areas and a positive physical length scale L.
    Unknowns are phi and L*q, so columns have comparable units. A common
    potential scale cancels from the relative diagnostics and needs no API knob.
    A bounded all-Neumann problem uses an area-weighted zero-mean gauge;
    exterior problems instead require a decaying field at infinity.

    Factorization defaults to one BLAS thread, based on measured dense setup
    costs. None preserves the caller's limits. The scoped limit is temporarily
    process-wide; concurrent applications should coordinate BLAS externally.
    Repeated solves retain the caller's ambient threading.
    """

    operators: LaplaceOperators
    dirichlet: np.ndarray
    areas: np.ndarray
    length_scale: float
    diagnostics: dict
    factorization_threads: int | None

    def __init__(self, operators, *, dirichlet, areas, length_scale, bounded,
                 factorization_threads=1):
        start = perf_counter()
        if not isinstance(operators, LaplaceOperators):
            raise ValueError("operators must be assembled LaplaceOperators")
        n = len(operators.single_layer)
        if n == 0:
            raise ValueError("at least one panel is required")
        for operator in (operators.single_layer, operators.double_layer):
            if (operator.shape != (n, n) or operator.dtype != np.dtype('float64')
                    or operator.flags.writeable or not np.all(np.isfinite(operator))):
                raise ValueError("operators must be finite read-only square float64 arrays")
        mask = np.asarray(dirichlet)
        if mask.dtype != np.bool_ or mask.shape != (n,):
            raise ValueError("dirichlet must be one boolean per panel")
        mask = mask.copy()
        mask.setflags(write=False)
        areas = array(areas, (n,), "areas")
        if np.any(areas <= 0):
            raise ValueError("panel areas must be positive")
        if isinstance(length_scale, (bool, np.bool_)) or not np.isscalar(length_scale):
            raise ValueError("length_scale must be positive and finite")
        length = float(length_scale)
        if not np.isfinite(length) or length <= 0:
            raise ValueError("length_scale must be positive and finite")
        if type(bounded) is not bool:
            raise ValueError("bounded must explicitly identify the fluid domain")
        if factorization_threads is not None and (
            isinstance(factorization_threads, (bool, np.bool_))
            or not isinstance(factorization_threads, (int, np.integer))
            or factorization_threads < 1
        ):
            raise ValueError("factorization_threads must be a positive integer or None")
        if np.any(np.diag(operators.double_layer) != 0):
            raise ValueError("flat-panel principal-value K must have zero diagonal")
        if factorization_threads is not None:
            factorization_threads = int(factorization_threads)
        weights = areas/areas.sum()
        weights.setflags(write=False)
        for name, value in dict(operators=operators, dirichlet=mask, areas=areas,
                                length_scale=length, _weights=weights,
                                factorization_threads=factorization_threads,
                                _has_dirichlet=bool(np.any(mask)),
                                _has_neumann=bool(np.any(~mask)),
                                _bordered=bounded and not np.any(mask)).items():
            object.__setattr__(self, name, value)
        size = n+int(self._bordered)
        # One Fortran-order system buffer becomes its LU factors in place.
        matrix = np.empty((size, size), order="F")
        _fill_mixed_matrix(operators.single_layer, operators.double_layer, mask, length, matrix)
        if self._bordered:
            matrix[:n, n] = 1
            matrix[n, :n] = self._weights
            matrix[n, n] = 0
        object.__setattr__(self, '_matrix_norm', dlange('1', matrix))
        with _blas.limit(limits=factorization_threads, user_api='blas'), warnings.catch_warnings():
            factor_pools = [pool for pool in _blas.info() if pool['user_api'] == 'blas']
            warnings.simplefilter("error", LinAlgWarning)
            try:
                lu, pivots = lu_factor(matrix, overwrite_a=True, check_finite=False)
            except LinAlgWarning as error:
                raise np.linalg.LinAlgError("singular boundary system") from error
            rcond, info = dgecon(lu, self._matrix_norm)
        if info != 0 or not np.isfinite(rcond) or rcond < 1e-12:
            raise np.linalg.LinAlgError(f"ill-conditioned boundary system (rcond={rcond:g})")
        lu.setflags(write=False)
        pivots.setflags(write=False)
        object.__setattr__(self, '_lu', lu)
        object.__setattr__(self, '_pivots', pivots)
        # Transposed C-order operators are already Fortran views. The infinity
        # norm of that view is the original 1-norm, without copies/abs matrices.
        # K has zero self diagonal; hence ||H||_1 = .5 + ||K||_1.
        h_norm = .5+dlange('I', operators.double_layer.T)
        s_norm = dlange('I', operators.single_layer.T)/length
        object.__setattr__(self, '_h_norm', h_norm)
        object.__setattr__(self, '_s_norm', s_norm)
        object.__setattr__(self, 'diagnostics', dict(
            panel_count=n, bordered=bool(self._bordered), rcond=float(rcond),
            setup_seconds=perf_counter()-start, factor_bytes=lu.nbytes+pivots.nbytes,
            factorization_blas=factor_pools))

    def solve(self, *, potential=None, normal_velocity=None):
        """Solve one or several RHS using physical full-panel prescribed fields.

        Arrays have shape (N,) or (N,R). Entries on panels where that field is
        unknown are ignored; supplied arrays must still be finite. Omit a field
        only when none of its entries are prescribed. Factors are reused.
        """
        start = perf_counter()
        n = len(self.areas)
        supplied = [np.asarray(value) for value in (potential, normal_velocity) if value is not None]
        if not supplied:
            raise ValueError("supply the prescribed boundary fields")
        shape = supplied[0].shape
        if len(shape) not in (1, 2) or shape[0] != n or (len(shape) == 2 and shape[1] == 0):
            raise ValueError("boundary fields must have shape (N,) or (N,R)")
        vector = len(shape) == 1
        width = 1 if vector else shape[1]

        # Copy once into mutable RHS fields; no panel/operator revalidation.
        def field(value, required, name):
            if value is None:
                if np.any(required):
                    raise ValueError(f"missing prescribed {name}")
                return np.zeros((n, width))
            values = array(value, shape, name)
            values.setflags(write=True)
            values = values.reshape(n, width)
            values[~required] = 0
            return values

        phi = field(potential, self.dirichlet, "potential")
        q = field(normal_velocity, ~self.dirichlet, "normal_velocity")
        compatibility = None
        if self._bordered:
            flux = self.areas @ q
            magnitude = self.areas @ np.abs(q)
            if np.any(np.abs(flux) > 1e-12*magnitude):
                raise ValueError("bounded all-Neumann data have incompatible net flux")
            compatibility = np.divide(np.abs(flux), magnitude, out=np.zeros(width),
                                       where=magnitude != 0).tolist()
        single, double = self.operators.single_layer, self.operators.double_layer
        rhs = single @ q if self._has_neumann else np.zeros_like(phi)
        if self._has_dirichlet:
            rhs -= double @ phi+.5*phi
        if self._bordered:
            rhs = np.vstack((rhs, np.zeros((1, width))))
        if not np.all(np.isfinite(rhs)):
            raise np.linalg.LinAlgError("non-finite boundary right-hand side")
        rhs_norm = np.sum(np.abs(rhs), axis=0)
        rhs = np.asfortranarray(rhs)
        unknown = lu_solve((self._lu, self._pivots), rhs, overwrite_b=True, check_finite=False)
        phi[~self.dirichlet] = unknown[:n][~self.dirichlet]
        q[self.dirichlet] = unknown[:n][self.dirichlet]/self.length_scale
        original = double @ phi+.5*phi-single @ q
        multiplier = unknown[n] if self._bordered else np.zeros(width)
        gauge = self._weights @ phi if self._bordered else np.zeros(width)
        residual = np.sum(np.abs(original+multiplier), axis=0)+np.abs(gauge)
        denominator = self._matrix_norm*np.sum(np.abs(unknown), axis=0)+rhs_norm
        backward = np.divide(residual, denominator, out=np.zeros(width), where=denominator != 0)
        if (not np.all(np.isfinite(unknown)) or not np.all(np.isfinite(backward))
                or np.any(backward > 1e-12)):
            raise np.linalg.LinAlgError("boundary solve failed the backward-error gate")
        reference = self._h_norm*np.sum(np.abs(phi), axis=0)+self._s_norm*self.length_scale*np.sum(np.abs(q), axis=0)
        bie = np.divide(np.sum(np.abs(original), axis=0), reference,
                        out=np.zeros(width), where=reference != 0)
        if not np.all(np.isfinite(bie)) or not np.all(np.isfinite(gauge)):
            raise np.linalg.LinAlgError("non-finite boundary diagnostics")
        net_flux = self.areas @ q
        if vector:
            phi, q = phi[:, 0], q[:, 0]
        phi.setflags(write=False)
        q.setflags(write=False)
        diagnostics = dict(solve_seconds=perf_counter()-start, rhs_count=width,
                           backward_error=backward.tolist(), original_bie_residual=bie.tolist(),
                           gauge_residual=gauge.tolist(), gauge_multiplier=multiplier.tolist(),
                           net_flux=net_flux.tolist(), compatibility_residual=compatibility)
        return BoundarySolution(phi, q, diagnostics)

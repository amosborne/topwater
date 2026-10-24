"""Quadrature-defined symmetric Green energy, fixed flat-domain experiment.

The inputs are RAW boundary-work blocks, not an evolution operator. The real
scalar Dirichlet energy defines its derivatives, including both Green terms.
Raw reciprocity defects remain diagnostics; this construction establishes
neither continuum accuracy nor positivity of the reduced water energy.
All energies/functionals here are divided by water density.
"""
import numpy as np
from scipy.linalg import cholesky, eigvalsh, solve_triangular


def _matrix(value, shape, name):
    value = np.asarray(value)
    if (value.shape != shape or np.iscomplexobj(value)
            or value.dtype.kind not in 'fiu' or not np.isfinite(value).all()):
        raise ValueError(name+' requires a finite real matrix of shape '+str(shape))
    return np.asarray(value, dtype=np.float64)


def _defect(numerator, denominator):
    return float(numerator/denominator) if denominator else (0. if numerator == 0 else None)


class OddRitzEnergy:
    """One mixed energy owner for water coefficients and six body velocities.

    E = integral grad(phi_a)*grad(phi_b) is defined by the symmetric Green
    integrand on the COMPLETE physical boundary. Its finite-quadrature blocks
    are (rawHH+rawHH.T)/2, (rawHR+rawRH.T)/2 and (rawRR+rawRR.T)/2.
    This is a new variational discretization, not a repair of a pre-existing D.

    The mixed functional is E(phi,phi)/2 - u.T*(BH*psi+BR*sigma).
    It is NOT the kinetic energy. Source-density stationarity uses its actual
    derivative, not an independently simplified sampled Neumann equation.
    The reference-area scaling changes coordinates only; no penalty is added.
    """
    def __init__(self, hh, hr, rh, rr, bh, br, source_areas):
        hh, rr = np.asarray(hh), np.asarray(rr)
        if hh.ndim != 2 or rr.ndim != 2 or min(hh.shape+rr.shape) <= 0:
            raise ValueError('nonempty square wave/source work blocks required')
        n, b = hh.shape[0], rr.shape[0]
        hh = _matrix(hh, (n, n), 'HH')
        hr = _matrix(hr, (n, b), 'HR')
        rh = _matrix(rh, (b, n), 'RH')
        rr = _matrix(rr, (b, b), 'RR')
        bh = _matrix(bh, (6, n), 'BH')
        br = _matrix(br, (6, b), 'BR')
        area = _matrix(source_areas, (b,), 'source areas')
        if np.any(area <= 0):
            raise ValueError('every source panel must have positive physical area')
        self.wave_count, self.source_count = n, b
        self.raw_diagnostics = dict(
            HH_adjoint_relative=_defect(np.linalg.norm(hh-hh.T), 2*np.linalg.norm(hh)),
            cross_adjoint_relative=_defect(np.linalg.norm(hr-rh.T), np.linalg.norm(hr)+np.linalg.norm(rh)),
            RR_adjoint_relative=_defect(np.linalg.norm(rr-rr.T), 2*np.linalg.norm(rr)))
        self.wave_energy = .5*(hh+hh.T)
        self.cross_energy = .5*(hr+rh.T)
        self.source_energy = .5*(rr+rr.T)
        self.BH, self.BR = bh.copy(), br.copy()
        self._area_scale = np.sqrt(area)
        scaled = self.source_energy/self._area_scale[:, None]/self._area_scale[None, :]
        self.source_eigenvalues = eigvalsh(scaled, check_finite=False)
        if self.source_eigenvalues[0] <= 0:
            raise ValueError('source energy is not positive definite; no rank or eigenvalue repair')
        self._source_factor = cholesky(scaled, lower=True, check_finite=False)
        self._wave_response = solve_triangular(self._source_factor,
            self.cross_energy.T/self._area_scale[:, None], lower=True, check_finite=False)
        self._body_response = solve_triangular(self._source_factor,
            self.BR.T/self._area_scale[:, None], lower=True, check_finite=False)
        # These are Schur complements/Gram products of the SAME scalar energy.
        # Do not symmetrize an evolution matrix or change any eigenvalues here.
        self.G = self.wave_energy-self._wave_response.T@self._wave_response
        self.C = self.BH-self._body_response.T@self._wave_response
        self.A = self._body_response.T@self._body_response
        self.reduction_symmetry_max = dict(
            G=float(abs(self.G-self.G.T).max()), A=float(abs(self.A-self.A.T).max()))
        for value in (self.wave_energy, self.cross_energy, self.source_energy,
                      self.BH, self.BR, self._area_scale, self._source_factor,
                      self._wave_response, self._body_response, self.G, self.C,
                      self.A, self.source_eigenvalues):
            value.setflags(write=False)

    def _fields(self, psi, u):
        psi, u = np.asarray(psi), np.asarray(u)
        if (psi.ndim not in (1, 2) or u.ndim != psi.ndim
                or psi.shape[0] != self.wave_count or u.shape[0] != 6
                or psi.shape[1:] != u.shape[1:]
                or not np.isfinite(psi).all() or not np.isfinite(u).all()):
            raise ValueError('matching finite wave fields and six body velocities required')
        return psi, u

    def solve_density(self, psi, u):
        psi, u = self._fields(psi, u)
        result = solve_triangular(self._source_factor.T,
            self._body_response@u-self._wave_response@psi,
            lower=False, check_finite=False)
        return result/self._area_scale if result.ndim == 1 else result/self._area_scale[:, None]

    def wave_flux(self, psi, u):
        """Canonical exposed-water MOMENT flux; eta_dot = Mphysical^-1 flux."""
        psi, u = self._fields(psi, u)
        return self.G@psi-self.C.T@u

    def body_port(self, psi, u):
        """Body potential port; at fixed geometry pressure load = -rho*z_dot."""
        psi, u = self._fields(psi, u)
        return self.C@psi+self.A@u

    def gradients(self, psi, sigma, u):
        """Derivatives of the unreduced mixed functional, all terms retained."""
        psi, u = self._fields(psi, u)
        sigma = np.asarray(sigma)
        if sigma.shape != (self.source_count,)+psi.shape[1:] or not np.isfinite(sigma).all():
            raise ValueError('source field must match wave/body columns')
        return (self.wave_energy@psi+self.cross_energy@sigma-self.BH.T@u,
                self.cross_energy.T@psi+self.source_energy@sigma-self.BR.T@u,
                -self.BH@psi-self.BR@sigma)

    def kinetic_energy_per_density(self, psi, sigma):
        """Quadrature-defined real kinetic functional, not time averaging."""
        psi, sigma = np.asarray(psi), np.asarray(sigma)
        if (psi.shape != (self.wave_count,) or sigma.shape != (self.source_count,)
                or np.iscomplexobj(psi) or np.iscomplexobj(sigma)
                or not np.isfinite(psi).all() or not np.isfinite(sigma).all()):
            raise ValueError('real physical wave/source vectors required for scalar energy')
        return float(.5*psi@(self.wave_energy@psi)
                     +psi@(self.cross_energy@sigma)
                     +.5*sigma@(self.source_energy@sigma))

    def mixed_functional(self, psi, sigma, u):
        """Kinetic energy per density minus prescribed body-potential work."""
        psi, u = self._fields(psi, u)
        if u.shape != (6,) or np.iscomplexobj(u):
            raise ValueError('one real six-component body velocity required')
        return self.kinetic_energy_per_density(psi, sigma)-float(u@(self.BH@psi+self.BR@sigma))

"""Causal linear water-wave PML falsifier; never a nonlinear solver boundary.

The SAME pointwise q12 sigma as the steady strip is realized with zero-history
auxiliaries. Implicit midpoint is a correctness reference, not an RK23 claim.
See Cohen/Imperiale, doi:10.4208/cicp.201109.261110s, sections4.2--4.5.
The first declared case is omega2/45deg, depth4, uncompressed width7.87758,
sigma8, prop+lowest-evanescent left trace, two time steps. No automatic family.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .pml_strip import _triangle_rule, stretching
from .sponge_strip import (DEPTH, GRAVITY, _line, make_strip, surface_matrices,
                           vertical_modes)

WIDTH = 7.877580409572782


def drive(t, omega, *, ramp_periods=1., plateau_periods=4., dc=.2):
    """C2 finite pulse with nonzero DC integral, then exactly zero drive."""
    period = 2*np.pi/omega
    ramp = ramp_periods*period; end = (2*ramp_periods+plateau_periods)*period
    def smooth(x):
        x = np.clip(x, 0., 1.)
        return x**3*(10-15*x+6*x*x)
    envelope = smooth(np.asarray(t)/ramp)*smooth((end-np.asarray(t))/ramp)
    return envelope*(np.cos(omega*np.asarray(t))+dc)


def left_profile(strip):
    """Fixed propagating plus .25 lowest-evanescent mode; no inferred port."""
    _, mu, _ = vertical_modes(strip.omega, 1)
    z = strip.points[:, 1]
    return np.cosh(strip.k*(z+DEPTH))/np.cosh(strip.k*DEPTH)+.25*np.cos(mu[0]*(z+DEPTH))


def stationary_field(strip, points, strength=8.):
    """Independent continuum harmonic field for THIS left Dirichlet/wall BC.

    The left amplitude is one, unlike the earlier unit-INCOMING modal port.
    DC pulse response is not represented by this nonzero-frequency formula.
    """
    points = np.asarray(points);left = strip.xlevels[0];width = strip.xlevels[-1]
    f, _, _ = stretching(points[..., 0], width, strength, strip.omega)
    end, _, _ = stretching(width, width, strength, strip.omega)
    _, mu, _ = vertical_modes(strip.omega, 1)
    alpha = np.sqrt(mu[0]**2+strip.ky**2)
    prop = (np.exp(1j*strip.kappa*(f-left))+np.exp(1j*strip.kappa*(2*end-f-left)))
    prop /= 1+np.exp(2j*strip.kappa*(end-left))
    evan = (np.exp(-alpha*(f-left))+np.exp(-alpha*(2*end-f-left)))
    evan /= 1+np.exp(-2*alpha*(end-left))
    z = points[..., 1]
    return (prop*np.cosh(strip.k*(z+DEPTH))/np.cosh(strip.k*DEPTH)
            +.25*evan*np.cos(mu[0]*(z+DEPTH)))


class MidpointStrip:
    """One fixed real factor per dt, three q-point histories, no projection."""
    def __init__(self, strip, dt, strength=8., order=12):
        if not np.isfinite(dt) or dt <= 0 or not np.isfinite(strength) or strength < 0:
            raise ValueError('positive finite dt and nonnegative sigma required')
        if len(strip.points) > 2000 or len(strip.triangles) > 250:
            raise ValueError('bounded strip requires <=2000nodes/250triangles')
        self.strip = strip;self.dt = dt;self.h = dt/2
        self.value, derivative, weight = _triangle_rule(order)
        corners = strip.points[strip.triangles]
        jac = corners[:, 1:]-corners[:, :1];det = np.linalg.det(jac)
        if np.any(det <= 0):
            raise ValueError('positive affine cells required')
        self.gradient = np.einsum('qaj,tij->tqai', derivative, np.linalg.inv(jac))
        self.measure = det[:, None]*weight
        x = np.einsum('qa,ta->tq', self.value, strip.points[strip.elements, 0])
        _, _, self.sigma = stretching(x, strip.xlevels[-1], strength, strip.omega)
        h = self.h;s = self.sigma;g = self.gradient;n = self.value
        local = np.einsum('tqa,tqb,tq->tab', g[..., 0], g[..., 0], self.measure/(1+h*s))
        local += np.einsum('tqa,tqb,tq->tab', g[..., 1], g[..., 1], self.measure*(1+h*s))
        local += strip.ky**2*np.einsum('qa,qb,tq->tab', n, n, self.measure*(1+h*s))
        rows = np.repeat(strip.elements, 15, axis=1).ravel()
        cols = np.tile(strip.elements, (1, 15)).ravel()
        size = len(strip.points)
        matrix = coo_matrix((local.ravel(), (rows, cols)), shape=(size, size)).tocsc()
        self.free, self.mass, self.weighted = surface_matrices(strip, strength, order=order)
        surface = (self.mass/h**2+self.weighted/h)/GRAVITY
        entries = surface.tocoo()
        matrix += coo_matrix((entries.data, (self.free[entries.row], self.free[entries.col])),
                             shape=matrix.shape).tocsc()
        self.left = np.unique(strip.left)
        self.unknown = np.setdiff1d(np.arange(size), self.left)
        self.boundary_block = matrix[self.unknown][:, self.left]
        self.factor = splu(matrix[self.unknown][:, self.unknown])
        self.profile = left_profile(strip)[self.left]
        self.chi = np.zeros_like(s);self.z = np.zeros_like(s);self.q = np.zeros_like(s)
        self.psi = np.zeros(len(self.free));self.v = self.psi.copy()
        self.potential = np.zeros(size);self.mid_v = self.v.copy();self.time = 0.

    def step(self):
        strip = self.strip;h = self.h;s = self.sigma;g = self.gradient
        local = np.einsum('tqa,tq->ta', g[..., 0], self.measure*s*self.chi/(1+h*s))
        local -= np.einsum('tqa,tq->ta', g[..., 1], self.measure*s*self.z)
        local -= strip.ky**2*np.einsum('qa,tq->ta', self.value, self.measure*s*self.q)
        rhs = np.zeros(len(strip.points));np.add.at(rhs, strip.elements.ravel(), local.ravel())
        rhs[self.free] += ((self.mass@self.psi)/h**2+(self.weighted@self.psi)/h
                          +(self.mass@self.v)/h)/GRAVITY
        end = self.time+self.dt
        boundary = .5*(drive(self.time, strip.omega)+drive(end, strip.omega))*self.profile
        mid = np.empty_like(rhs);mid[self.left] = boundary
        mid[self.unknown] = self.factor.solve(rhs[self.unknown]-self.boundary_block@boundary)
        if not np.all(np.isfinite(mid)):
            raise FloatingPointError('nonfinite time PML solution')
        grad = np.einsum('tqai,ta->tqi', g, mid[strip.elements])
        field = np.einsum('qa,ta->tq', self.value, mid[strip.elements])
        self.chi = 2*(self.chi+h*grad[..., 0])/(1+h*s)-self.chi
        self.z += 2*h*grad[..., 1];self.q += 2*h*field
        vmid = (mid[self.free]-self.psi)/h
        self.v = 2*vmid-self.v;self.psi = 2*mid[self.free]-self.psi
        self.potential = mid
        self.mid_v = vmid.copy()
        self.time = end
        return self.potential.copy()

    def diagnostics(self):
        """Quiet PHYSICAL energy only; auxiliary norms are not physical energy."""
        strip = self.strip
        grad = np.einsum('tqai,ta->tqi', self.gradient, self.potential[strip.elements])
        phi = np.einsum('qa,ta->tq', self.value, self.potential[strip.elements])
        x = np.einsum('qa,ta->tq', self.value, strip.points[strip.elements, 0])
        kinetic = .5*np.sum(self.measure*(x < 0)*(np.sum(grad**2, axis=-1)+strip.ky**2*phi**2))
        t, weight, shape = _line(12);length = np.diff(strip.xlevels)
        lookup = np.full(len(strip.points), -1);lookup[self.free] = np.arange(len(self.free))
        v = np.einsum('qa,ta->tq', shape, self.mid_v[lookup[strip.top]])
        surface_x = strip.xlevels[:-1, None]+length[:, None]*t
        gravity = .5/GRAVITY*np.sum(length[:, None]*weight*(surface_x < 0)*v**2)
        return np.r_[kinetic+gravity, [np.max(abs(a)) for a in (self.chi, self.z, self.q)]]


def run(dt, output):
    """One declared pulse/tail, not an automatically qualified boundary."""
    output = Path(output)
    if output.exists() or output.with_suffix('.npz').exists():
        raise FileExistsError('preserve existing strip outputs')
    strip = make_strip(2., np.pi/4, WIDTH, spacing=np.pi/6, layers=4, quiet_wavelengths=.25)
    period = 2*np.pi/strip.omega;duration = 12*period
    if not np.isfinite(dt) or dt <= 0 or duration/dt > 20000 or dt > period/20:
        raise ValueError('declared strip requires20..1666steps/period, bounded saved histories')
    solver = MidpointStrip(strip, dt)
    count = int(np.ceil(duration/dt));times = [];fields = [];diagnostics = []
    endpoint_times = [0.];surface_v = [solver.v.copy()];surface_psi = [solver.psi.copy()]
    for _ in range(count):
        fields.append(solver.step());times.append(solver.time-dt/2);diagnostics.append(solver.diagnostics())
        endpoint_times.append(solver.time);surface_v.append(solver.v.copy());surface_psi.append(solver.psi.copy())
    times = np.asarray(times);fields = np.asarray(fields)
    fit = (times >= 3*period)&(times <= 5*period)
    design = np.column_stack((np.cos(strip.omega*times[fit]), np.sin(strip.omega*times[fit]),
                              np.ones(np.count_nonzero(fit))))
    coefficients = np.linalg.lstsq(design, fields[fit], rcond=None)[0]
    harmonic = coefficients[0]+1j*coefficients[1]
    exact = stationary_field(strip, strip.points)
    # A finite-duration fit can retain startup transients. Preserve rather than
    # explain away that error; this is not a steady-reference qualification.
    report = {'status': 'UNQUALIFIED time-domain strip falsifier', 'dt': dt,
              'omega': strip.omega, 'angle': strip.angle, 'depth': DEPTH,
              'width': WIDTH, 'sigma_terminal': 8., 'volume_quadrature': 12,
              'nodes': len(strip.points), 'triangles': len(strip.triangles),
              'pulse_ramp_periods': 1., 'pulse_plateau_periods': 4., 'pulse_dc': .2,
              'duration_periods': 12., 'fit_period_range': [3., 5.],
              'harmonic_nodal_error_including_transients': float(np.linalg.norm(harmonic-exact)/np.linalg.norm(exact)),
              'diagnostic_columns': ['quiet_energy_per_density', 'chi_max', 'z_max', 'q_max'],
              'zero_history': True, 'no_mean_or_history_reset': True,
              'potential_and_energy_location': 'actual midpoint solve',
              'surface_state_location': 'accepted fixed-step endpoints',
              'boundary_v_checkerboard_is_retained': True}
    np.savez(output.with_suffix('.npz'), times=times, potential=fields, diagnostics=diagnostics,
             points=strip.points, triangles=strip.triangles, harmonic=harmonic, exact=exact,
             endpoint_times=endpoint_times, surface_psi=surface_psi, surface_v=surface_v,
             free_nodes=solver.free, left_nodes=solver.left,
             final_chi=solver.chi, final_z=solver.z, final_q=solver.q)
    output.write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        print(json.dumps(run(args.dt, args.output), indent=2))

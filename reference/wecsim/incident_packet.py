"""One shared dispersive incident field for native excitation and rendering."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat
import xarray as xr
from threadpoolctl import threadpool_limits


class IncidentPacket:
    def __init__(self, case, controls):
        sea = case['sea']
        self.a = np.asarray(sea['amplitudes_m'])
        self.w0 = 2*np.pi*np.asarray(sea['frequencies_hz'])
        self.phi = np.asarray(sea['phases_rad'])
        self.theta = np.deg2rad(sea['propagation_headings_degrees'])
        self.g = case['water']['gravity_m_s2']
        self.B = np.asarray(controls['half_bandwidth_rad_s'])
        self.tc = np.asarray(controls['focus_time_s'])
        n = controls['spectral_samples_per_component']
        nu = self.B[:, None]*np.linspace(-1, 1, n)
        self.w = self.w0[:, None]+nu
        self.k = self.w**2/self.g
        weights = np.cos(np.pi*nu/(2*self.B[:, None]))**4
        weights[:, [0, -1]] = 0
        self.weights = weights/weights.sum(axis=1)[:, None]
        self.c = self.a[:, None]*self.weights*np.exp(1j*(self.phi[:, None]+nu*self.tc[:, None]))
        self.normal = np.c_[np.cos(self.theta), np.sin(self.theta)]

    def elevation(self, xy, t):
        xy = np.atleast_2d(xy)
        result = np.zeros((len(xy), 3))
        for j in range(3):
            s = xy@self.normal[j]
            for start in range(0, len(xy), 512):
                end = start+512
                z = np.exp(1j*(s[start:end, None]*self.k[j]-self.w[j]*t))
                result[start:end, j] = (z@self.c[j]).real
        return result

    def initial_bounds(self, xyz, spacing, deadline):
        max_eta = max_speed = max_grad_eta = max_grad_velocity = 0.
        for start in range(0, len(xyz), 256):
            if time.monotonic() > deadline:
                raise TimeoutError('Declared incident diagnostic runtime exceeded')
            pts = xyz[start:start+256]
            eta = np.zeros(len(pts))
            velocity = np.zeros((len(pts), 3))
            de = np.zeros((len(pts), 2))
            dv = np.zeros((len(pts), 3, 3))
            for j in range(3):
                phase = np.exp(1j*((pts[:, :2]@self.normal[j])[:, None]*self.k[j]))
                amp = phase*self.c[j]
                eta += amp.sum(axis=1).real
                for d in range(2):
                    de[:, d] += (amp*(1j*self.k[j]*self.normal[j, d])).sum(axis=1).real
                water_amp = amp*np.exp(pts[:, 2, None]*self.k[j])
                component = np.array([*self.normal[j], -1j])
                factors = [1j*self.k[j]*self.normal[j, 0],
                           1j*self.k[j]*self.normal[j, 1], self.k[j]]
                for q in range(3):
                    velocity[:, q] += (water_amp*self.w[j]*component[q]).sum(axis=1).real
                    for d in range(3):
                        dv[:, q, d] += (water_amp*self.w[j]*component[q]*factors[d]).sum(axis=1).real
            max_eta = max(max_eta, abs(eta).max())
            max_speed = max(max_speed, np.linalg.norm(velocity, axis=1).max())
            max_grad_eta = max(max_grad_eta, np.linalg.norm(de, axis=1).max())
            max_grad_velocity = max(max_grad_velocity, np.linalg.norm(dv, axis=(1, 2)).max())
        r2 = np.sqrt(2)*spacing/2
        r3 = np.sqrt(3)*spacing/2
        ehess = np.sum(self.a[:, None]*self.weights*self.k**2)
        vhess = 2*np.sqrt(2)*np.sum(self.a[:, None]*self.weights*self.w*self.k**2)
        return dict(samples=len(xyz), grid_spacing_m=spacing,
                    sampled_max_abs_eta_m=float(max_eta), sampled_max_speed_m_s=float(max_speed),
                    eta_region_bound_m=float(max_eta+r2*max_grad_eta+0.5*r2*r2*ehess),
                    velocity_region_bound_m_s=float(max_speed+r3*max_grad_velocity+0.5*r3*r3*vhess),
                    bound='sample maximum + local gradient radius + global spectral Hessian remainder')


def excitation(dataset):
    value = dataset['diffraction_force']+dataset['Froude_Krylov_force']
    value = value.sel(influenced_dof=['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])
    value = value.transpose('influenced_dof', 'wave_direction', 'omega', 'complex')
    return value.sel(complex='re').values+1j*value.sel(complex='im').values


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--controls', type=Path, default=Path('reference/wecsim/packet-controls.json'))
    args = p.parse_args()
    started = time.monotonic()
    case_path = Path('reference/capytaine/canonical-boat.json')
    native = Path('artifacts/reference/minimum-20261024/native-gl2')
    consumer = Path('artifacts/reference/minimum-20261024/consumer')
    case = json.loads(case_path.read_text())
    controls = json.loads(args.controls.read_text())
    field = IncidentPacket(case, controls)
    args.output.mkdir(parents=True, exist_ok=False)
    dt = controls['dt_s']
    tf = np.arange(round((controls['forcing_time_s'][1]-controls['forcing_time_s'][0])/dt)+1)*dt+controls['forcing_time_s'][0]
    to = np.arange(round(10/dt)+1)*dt
    eta = np.array([field.elevation([[0, 0]], t)[0] for t in tf])
    for j in range(3):
        np.savetxt(args.output/f'wave_{j+1}.txt', np.c_[tf, eta[:, j]], fmt='%.16e')
    np.savez(args.output/'packet.npz', w=field.w, k=field.k, c=field.c,
             weights=field.weights, direction=field.theta, times=tf, elevation=eta,
             carrier_omega=field.w0, carrier_amplitude=field.a, focus_time=field.tc)
    wet = np.asarray(case['minimum_wetted_triangles_world_m']).reshape(-1, 3)
    lo = wet.min(axis=0)
    hi = wet.max(axis=0)
    hi[2] = 0
    h = controls['initial_grid_spacing_m']
    axes = [np.linspace(lo[d], hi[d], int(np.ceil((hi[d]-lo[d])/h))+1) for d in range(3)]
    xyz = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    bounds = field.initial_bounds(xyz, h, started+controls['maximum_elapsed_seconds'])
    bounds['region_min_xyz_m'] = lo.tolist()
    bounds['region_max_xyz_m'] = hi.tolist()
    ds = xr.open_dataset(native/'canonical.nc')
    H = excitation(ds)
    omega = ds.omega.values
    directions = ds.wave_direction.values
    direct = np.zeros((len(to), 6))
    scale = np.zeros(6)
    for j in range(3):
        idir = int(np.argmin(abs(directions-field.theta[j])))
        assert abs(directions[idir]-field.theta[j]) < 1e-12
        for q in range(6):
            hq = np.interp(field.w[j], omega, H[q, idir].real)+1j*np.interp(field.w[j], omega, H[q, idir].imag)
            direct[:, q] += (np.exp(-1j*to[:, None]*field.w[j])@(field.c[j]*hq)).real
            hcarrier = np.interp(field.w0[j], omega, H[q, idir].real)+1j*np.interp(field.w0[j], omega, H[q, idir].imag)
            scale[q] += field.a[j]*abs(hcarrier)
    bemio = loadmat(consumer/'bemio.mat', simplify_cells=True)
    kernel = bemio['h']['ex_K']*case['water']['density_kg_m3']*case['water']['gravity_m_s2']
    ktime = bemio['h']['ex_t']
    assert np.allclose(ktime, np.arange(-20, 20+dt/2, dt))
    actual = sum(np.stack([np.convolve(eta[:, j], kernel[q, j], mode='valid')*dt
                          for q in range(6)], axis=1) for j in range(3))
    assert actual.shape == direct.shape
    initial_actual = abs(actual[0])
    initial_direct = abs(direct[0])
    tiny = scale < 1e-12
    limits = np.minimum(controls['maximum_initial_force_N_or_Nm'], controls['calm_fraction']*scale)
    arrival = np.max(abs(eta[(tf >= 0) & (tf <= 10)]), axis=0)/field.a
    freqmax = min(.6, np.sqrt(field.g/(2*np.pi*(2*9*.3048))))
    frac = controls['calm_fraction']
    eta_limit = frac*field.a.sum()
    velocity_limit = frac*np.sum(field.a*field.w0)
    gates = dict(eta=bool(bounds['eta_region_bound_m'] <= eta_limit),
                 velocity=bool(bounds['velocity_region_bound_m_s'] <= velocity_limit),
                 force_actual=bool(np.all(initial_actual <= limits)),
                 force_direct=bool(np.all(initial_direct <= limits)),
                 arrival=bool(np.any(arrival >= controls['minimum_arrival_fraction'])),
                 allowed_band=bool(field.w.min() >= omega.min() and field.w.max() <= 2*np.pi*freqmax))
    report = dict(complete=True, accepted=bool(all(gates.values())), gates=gates,
                  initial=bounds, initial_eta_limit_m=float(eta_limit),
                  initial_speed_limit_m_s=float(velocity_limit),
                  initial_actual_force=initial_actual.tolist(), initial_spectral_force=initial_direct.tolist(),
                  force_limits=limits.tolist(), carrier_force_scale=scale.tolist(), degenerate_axes=tiny.tolist(),
                  arrival_fraction_by_component=arrival.tolist(),
                  spectral_frequency_range_hz=(np.array([field.w.min(), field.w.max()])/(2*np.pi)).tolist(),
                  allowed_maximum_frequency_hz=float(freqmax),
                  convolution_vs_spectral_peak_error=np.max(abs(actual-direct), axis=0).tolist(),
                  carrier_amplitudes_m=field.a.tolist(), carrier_frequencies_hz=(field.w0/(2*np.pi)).tolist(),
                  body_initial_state='q=v=radiation memory exactly zero; no body simulated',
                  elapsed_s=time.monotonic()-started,
                  input_sha256={str(v):hashlib.sha256(v.read_bytes()).hexdigest() for v in [case_path,args.controls]})
    np.savez(args.output/'startup-forces.npz', time=to, force_convolution=actual, force_spectral=direct)
    (args.output/'diagnostic.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
    return 0 if report['accepted'] else 2


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())

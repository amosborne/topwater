"""Bounded canonical-sea bandwidth and incident-only microbenchmark; no PDE.

Query clouds are workload controls, not hull meshes or fluid discretizations.
Nine outputs: eta, eta_x, eta_y, eta_t, phi, phi_t, u_x, u_y, u_z.
"""
import hashlib
import json
import math
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
from numba import njit
from threadpoolctl import threadpool_info, threadpool_limits

from topwater import sea_states
from tests.benchmark_rankine import peak_rss_bytes

G, L, B, RHO = 9.81, 24*.3048, 9*.3048, 1025.


def coefficients(waves):
    a, omega, direction = waves.amplitudes, waves.angular_frequencies, waves.directions
    k = omega**2/G
    return np.column_stack((a, omega, k, np.cos(direction), np.sin(direction), waves.phases))


def vectorized(points, values, time):
    a, omega, k, dx, dy, phase0 = values.T
    phase = points[:, :1]*(k*dx)+points[:, 1:2]*(k*dy)-time*omega+phase0
    sine, cosine = np.sin(phase), np.cos(phase)
    decay = np.exp(points[:, 2:3]*k)
    out = np.column_stack((cosine@a, sine@(-a*k*dx), sine@(-a*k*dy), sine@(a*omega),
        (sine*decay)@(G*a/omega), (cosine*decay)@(-G*a),
        (cosine*decay)@(a*omega*dx), (cosine*decay)@(a*omega*dy), (sine*decay)@(a*omega)))
    return out


@njit(fastmath=False)
def fused_cpu(points, values, time):
    out = np.zeros((len(points), 9))
    for i in range(len(points)):
        for j in range(len(values)):
            a, omega, k, dx, dy, phase0 = values[j]
            phase = k*(dx*points[i, 0]+dy*points[i, 1])-time*omega+phase0
            sine, cosine, decay = math.sin(phase), math.cos(phase), math.exp(k*points[i, 2])
            out[i, 0] += a*cosine
            out[i, 1] -= a*k*dx*sine
            out[i, 2] -= a*k*dy*sine
            out[i, 3] += a*omega*sine
            out[i, 4] += G*a/omega*decay*sine
            out[i, 5] -= G*a*decay*cosine
            out[i, 6] += a*omega*dx*decay*cosine
            out[i, 7] += a*omega*dy*decay*cosine
            out[i, 8] += a*omega*decay*sine
    return out


CUDA = r'''
extern "C" __global__ void incident(const double* points, const double* v,
                                   int np, int nc, double time, double* out) {
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=np) return;
  double accum[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double x=points[3*i],y=points[3*i+1],z=points[3*i+2];
  for(int j=0;j<nc;j++) {
    double a=v[6*j],w=v[6*j+1],k=v[6*j+2],dx=v[6*j+3],dy=v[6*j+4];
    double phase=k*(dx*x+dy*y)-time*w+v[6*j+5],s,c;
    sincos(phase,&s,&c); double e=exp(k*z);
    accum[0]+=a*c; accum[1]-=a*k*dx*s; accum[2]-=a*k*dy*s; accum[3]+=a*w*s;
    accum[4]+=9.81*a/w*e*s; accum[5]-=9.81*a*e*c;
    accum[6]+=a*w*dx*e*c; accum[7]+=a*w*dy*e*c; accum[8]+=a*w*e*s;
  }
  for(int j=0;j<9;j++) out[9*i+j]=accum[j];
}'''


def clock(function, repeats=7):
    function()
    seconds = []
    for i in range(repeats):
        start = perf_counter()
        function()
        seconds.append(perf_counter()-start)
    return dict(median_seconds=float(np.median(seconds)), minimum_seconds=float(min(seconds)),
                repeats=repeats)


def weighted_quantiles(value, weight, quantiles=(.005, .05, .5, .95, .995)):
    order = np.argsort(value, kind='stable')
    cumulative = np.cumsum(weight[order])/np.sum(weight)
    return [float(value[order[min(np.searchsorted(cumulative, q), len(value)-1)]]) for q in quantiles]


def shortest_arc(weight, fraction):
    weight = weight/weight.sum()
    prefix = np.r_[0., np.cumsum(np.tile(weight, 2))]
    choices = []
    for start in range(len(weight)):
        stop = int(np.searchsorted(prefix, prefix[start]+fraction, side='left'))-1
        choices.append((stop-start, start, stop))
    width, start, stop = min(choices)
    return dict(start_from_deg=start*.5, end_from_deg=(stop*.5) % 360, span_deg=width*.5,
                crossing_north=stop >= len(weight), fraction=fraction,
                resolution_deg=.5)


def main():
    report_path = Path('artifacts/canonical-sea-incident-benchmark.json')
    if report_path.exists():
        raise FileExistsError('preserve existing benchmark')
    start = perf_counter()
    state = sea_states.observed()
    waves = sea_states.realize(state, beam_m=B, seed=42, gravity=G)
    variance = state.band_amplitudes_m**2/2
    edges = state.frequency_edges_hz
    cdf = np.r_[0., np.cumsum(variance)]/variance.sum()
    q = np.array([.005, .05, .5, .95, .995])
    source_frequency = np.interp(q, cdf, edges)
    directional, scales = zip(*(sea_states._directional_weights(m) for m in state.direction_moments))
    directional = np.asarray(directional)
    angular_energy = variance@directional
    bearing = np.arange(720)*np.pi/360
    z = np.sum(angular_energy*np.exp(1j*bearing))/variance.sum()
    report = dict(qualified=False,
        scope='Exact analytic linear incident deep-water fields only; no total fluid pressure, '
              'body response, near/far coupling, trajectory or end-to-end throughput qualification.',
        source_time_utc=state.time_utc, observation_id=state.id, seed=42,
        source_library_sha256=state.provenance['source_library_sha256'],
        catalog_sha256=hashlib.sha256(Path('src/topwater/data/socal_sea_states.json').read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        boat=dict(length_m=L, beam_m=B, mass_kg=2500., heading='bow east; no heading optimization'),
        gravity=G, density=RHO, hm0_m=float(state.hm0_m), tm02_s=state.tm02_s,
        quantile_probabilities=q.tolist(),
        measured_band_frequency_quantiles_hz=source_frequency.tolist(),
        measured_band_wavelengths_at_quantiles_m=(G/(2*np.pi*source_frequency**2)).tolist(),
        measured_band_kL_at_quantiles=(source_frequency**2*4*np.pi**2/G*L).tolist(),
        measured_band_kB_at_quantiles=(source_frequency**2*4*np.pi**2/G*B).tolist(),
        frequency_quantile_model='Exact overlap integral under the catalog flat-density within-band model.',
        measured_coverage_hz=[float(edges[0]), float(edges[-1])],
        positive_energy_upper_hz=float(edges[1:][variance > 0].max()), unmeasured_tails='unknown',
        directional_reconstruction=dict(mean_from_deg=float(np.angle(z) % (2*np.pi)*180/np.pi),
            first_harmonic_concentration=float(abs(z)),
            shortest_arcs=[shortest_arc(angular_energy, p) for p in (.9, .99)],
            regularization_min_scale=float(min(scales)),
            scope='MEM reconstructed direction PDF from four observed moments; not uniquely measured PDF.'),
        components=[], timings=[], runtime=dict(python=platform.python_version(), numpy=np.__version__,
            native_threadpools=threadpool_info()),
        timing_scope='Single process, BLAS1; seven warmed repeats; cold setup separate; no concurrent PDE jobs.')
    # Derivative moments reveal why elevation-energy cropping alone is not an
    # accuracy argument for slope, acceleration or hull force.
    for name, power in (('elevation', 0), ('surface_slope_magnitude', 4), ('eta_time_derivative', 2)):
        density = variance/np.diff(edges)
        moment = density*(edges[1:]**(power+1)-edges[:-1]**(power+1))/(power+1)
        cumulative = np.r_[0., np.cumsum(moment)]
        target = q*cumulative[-1]
        indices = np.clip(np.searchsorted(cumulative, target, side='right')-1, 0, len(variance)-1)
        within = edges[indices]**(power+1)+(target-cumulative[indices])*(power+1)/np.maximum(density[indices], 1e-300)
        frequency = within**(1/(power+1))
        report[name+'_frequency_quantiles_hz'] = frequency.tolist()
    print('bandwidth ready', flush=True)
    # CUDA runtime discovery is already owned by the existing experimental helper.
    try:
        before = perf_counter()
        from tests.experiments.free_surface.cuda_operator import cp
        kernel = cp.RawKernel(CUDA, 'incident', options=('--std=c++11',))
        report['runtime']['gpu_name'] = cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
        report['runtime']['cupy'] = cp.__version__
        report['gpu_import_seconds'] = perf_counter()-before
    except Exception as error:
        cp, kernel = None, None
        report['gpu_unavailable'] = dict(type=type(error).__name__, message=str(error))
    rng = np.random.default_rng(20261022)
    for per_band in (1, 2, 4, 8):
        case = sea_states.realize(state, beam_m=B, seed=42, gravity=G, components_per_band=per_band)
        values = coefficients(case)
        theta = (1.5*np.pi-case.directions.reshape(48, per_band)) % (2*np.pi)
        realized_moments = np.stack((np.cos(theta), np.sin(theta), np.cos(2*theta), np.sin(2*theta)), axis=-1).mean(axis=1)
        target_moments = np.column_stack([directional@function(bearing) for function in
            (np.cos, np.sin, lambda x: np.cos(2*x), lambda x: np.sin(2*x))])
        report['components'].append(dict(per_band=per_band, components=len(values),
            variance_relative_error=float(abs(np.sum(case.amplitudes**2)/2-state.variance_m2)/state.variance_m2),
            directional_moment_maximum_absolute_error=float(abs(realized_moments-target_moments).max()),
            directional_moment_energy_weighted_mean_absolute_error=float(np.sum(variance[:, None]*abs(realized_moments-target_moments))/(4*variance.sum())),
            scope='Different component counts are different realizations, not pointwise approximations of seed42/8perband.'))
        device_values = cp.asarray(values) if cp is not None else None
        for npoints in (256, 1024, 4096):
            if perf_counter()-start > 100:
                report['budget_exhausted'] = True
                break
            points = rng.uniform([-L, -B, -2.], [L, B, 0.], size=(npoints, 3))
            reference = vectorized(points, values, .37)
            before = perf_counter()
            actual = fused_cpu(points, values, .37)
            cold = perf_counter()-before
            error = float(np.max(abs(actual-reference)))
            if not np.allclose(actual, reference, rtol=3e-12, atol=3e-13):
                raise ValueError('CPU fused incident field differs')
            row = dict(points=npoints, components=len(values),
                cpu_fused_cold_seconds=cold, cpu_fused_maximum_absolute_error=error,
                numpy=clock(lambda: vectorized(points, values, .37)),
                cpu_fused=clock(lambda: fused_cpu(points, values, .37)),
                numpy_largest_pair_array_bytes=npoints*len(values)*8,
                cpu_fused_output_bytes=npoints*9*8)
            if cp is not None:
                device_points, device_out = cp.asarray(points), cp.empty((npoints, 9), dtype=cp.float64)
                def gpu():
                    kernel(((npoints+127)//128,), (128,), (device_points, device_values,
                        np.int32(npoints), np.int32(len(values)), np.float64(.37), device_out))
                    cp.cuda.get_current_stream().synchronize()
                before = perf_counter()
                gpu()
                row['gpu_cold_seconds'] = perf_counter()-before
                device_actual = device_out.get()
                if not np.allclose(device_actual, reference, rtol=3e-12, atol=3e-13):
                    raise ValueError('GPU fused incident field differs')
                row['gpu_maximum_absolute_error'] = float(abs(device_actual-reference).max())
                row['gpu_resident'] = clock(gpu)
                def transfer_inclusive():
                    device_points.set(points)
                    gpu()
                    device_out.get()
                row['gpu_points_upload_output_download'] = clock(transfer_inclusive)
            report['timings'].append(row)
            print(json.dumps(dict(points=npoints, components=len(values),
                cpu=row['cpu_fused']['median_seconds'], gpu=row.get('gpu_resident'))), flush=True)
            report.update(elapsed_seconds=perf_counter()-start, peak_rss_bytes=peak_rss_bytes())
            report_path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    report.update(complete=not report.get('budget_exhausted', False), elapsed_seconds=perf_counter()-start,
                  peak_rss_bytes=peak_rss_bytes())
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    with threadpool_limits(1):
        main()

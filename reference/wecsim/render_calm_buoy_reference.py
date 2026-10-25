"""Preview/encode an accepted native60s case and its exact incoming train.

Source preparation does not authorize execution. Preview and encoding require
separate coordinator allocations. The visual cache never supplies body forcing.
"""
import os
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import argparse
import ctypes
import importlib.metadata
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation, cm, colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from calm_buoy_train import CalmBuoyTrain
from render_generated_reference import ROOT, read, write, sha, resolved, check_deadline, motion_plot


def load_inputs(config):
    if config['status'] != 'ready_for_coordinator_render_allocation':
        raise ValueError('Accepted packet/motion identities have not been bound')
    required = ('physical_case', 'heading_file', 'packet_npz', 'incident_report',
                'motion_mat', 'summary_json', 'native_output_mat', 'motion_acceptance_receipt')
    identities = config['input_sha256']
    for key in required:
        assert config[key] and config[key] in identities, f'Unbound input: {key}'
    for name in ('render_calm_buoy_reference.py', 'render_generated_reference.py', 'calm_buoy_train.py'):
        assert 'reference/wecsim/'+name in identities
    assert all(digest and sha(resolved(path)) == digest for path, digest in identities.items())
    receipt = read(resolved(config['motion_acceptance_receipt']))
    assert receipt['decision'] == 'accepted_for_visualization' and receipt['reviewed_by']
    assert receipt['motion_sha256'] == identities[config['motion_mat']]
    assert receipt['summary_sha256'] == identities[config['summary_json']]
    assert receipt['native_output_sha256'] == identities[config['native_output_mat']]
    assert receipt['heading_sha256'] == identities[config['heading_file']]
    assert receipt['packet_sha256'] == identities[config['packet_npz']]
    assert receipt['incident_report_sha256'] == identities[config['incident_report']]
    assert receipt['trajectory_status'] in ('provisional_model_trajectory', 'accepted_reference_case')
    case = read(resolved(config['physical_case']))
    heading = read(resolved(config['heading_file']))
    summary = read(resolved(config['summary_json']))
    incident = read(resolved(config['incident_report']))
    assert incident['complete'] and incident['accepted']
    assert summary['complete'] and summary['samples'] == 3001 and summary['duration_s'] == 60
    assert heading['observation_id'] == config['observation_id'] == receipt['observation_id']
    assert heading['seed'] == config['seed'] == 42
    train = CalmBuoyTrain(resolved(config['packet_npz']))
    assert np.isfinite(train.C).all() and np.isfinite(train.omega).all()
    assert train.omega[0] >= 0 and np.all(train.omega[1:] > 0)
    assert np.all(np.diff(train.omega) > 0)
    if train.omega[0] == 0:
        assert np.array_equal(train.C[:, :, 0], np.zeros_like(train.C[:, :, 0])), 'Nonzero DC wave coefficient'
    assert np.array_equal(train.data['signs'], [1., -1.])
    assert np.array_equal(train.data['component_id'], np.arange(384))
    waves = heading['boat_frame_waves']
    for key, wave_key in (('carrier_amplitudes_m', 'amplitudes'),
                          ('carrier_omega_rad_s', 'angular_frequencies'),
                          ('carrier_phases_body_rad', 'phases'),
                          ('directions_body_rad', 'directions')):
        assert np.array_equal(train.data[key], np.asarray(waves[wave_key]))
    heading_rotation = np.asarray(heading['global_from_boat_rotation'])
    translation = np.asarray(heading['global_from_boat_translation_m'])
    assert np.array_equal(train.data['global_from_boat_rotation'], heading_rotation)
    assert np.array_equal(train.data['global_from_boat_translation_m'], translation)
    yaw = float(heading['mean_box_yaw_global_rad'])
    np.testing.assert_allclose(heading_rotation, Rotation.from_euler('z', yaw).as_matrix(), rtol=0, atol=2e-14)
    assert translation.shape == (3,) and translation[2] == 0
    motion = loadmat(resolved(config['motion_mat']), simplify_cells=True)
    raw = loadmat(resolved(config['native_output_mat']), simplify_cells=True)['native_body']
    t = np.asarray(motion['time']).ravel()
    q, position = np.asarray(motion['q']), np.asarray(motion['position_world'])
    assert t.shape == (3001,) and q.shape == position.shape == (3001, 6)
    assert np.isfinite(t).all() and np.isfinite(q).all() and np.isfinite(position).all()
    np.testing.assert_allclose(t, np.arange(3001)*.02, rtol=0, atol=1e-9)
    assert np.array_equal(np.asarray(raw['time']).ravel(), t)
    values = np.asarray(raw['signals']['values'])
    assert values.shape == (3001, 60) and np.isfinite(values).all()
    assert np.array_equal(values[:, :6], position)
    cg = np.asarray(case['world_cg_m'])
    np.testing.assert_allclose(np.asarray(motion['mean_world_cg']).ravel(), cg, rtol=0, atol=2e-12)
    np.testing.assert_allclose(position[:, :3]-cg, q[:, :3], rtol=0, atol=2e-12)
    np.testing.assert_allclose(position[:, 3:], q[:, 3:], rtol=0, atol=2e-12)
    hull = case['minimum_hull']
    eq = hull['equilibrium_quaternion_wxyz']
    equilibrium = Rotation.from_quat([*eq[1:], eq[0]]).as_matrix()
    local = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@equilibrium.T
    faces = np.asarray([f for i, f in enumerate(hull['triangular_faces'])
                        if i not in hull['nonphysical_closure_faces']], dtype=int)
    return dict(case=case, heading=heading, summary=summary, receipt=receipt,
        train=train, time=t, q=q, position=position, cg=cg, local=local, faces=faces,
        heading_rotation=heading_rotation, translation=translation, identities=identities)


def prepare_frames(data, config, started, budget):
    indices = np.arange(0, 3001, 4)
    times = data['time'][indices]
    np.testing.assert_allclose(times, np.arange(751)/12.5, rtol=0, atol=1e-9)
    rotations = Rotation.from_euler('xyz', data['position'][indices, 3:]).as_matrix()
    native = np.einsum('fij,vj->fvi', rotations, data['local'])+data['position'][indices, None, :3]
    hull = native@data['heading_rotation'].T+data['translation']
    centers = (hull[:, :, :2].min(axis=1)+hull[:, :, :2].max(axis=1))/2
    width, nodes = float(config['viewport_width_xy_m']), int(config['surface_grid_nodes'])
    assert width == 12 and nodes == 65
    dx, dy = np.meshgrid(np.linspace(-width/2, width/2, nodes), np.linspace(-width/2, width/2, nodes))
    offsets = np.c_[dx.ravel(), dy.ravel()]
    # Spatial field coordinates are in the mean-boat frame, not instantaneous
    # body coordinates. Applying body rotations to waves would move the ocean.
    boat_centers = (np.c_[centers, np.zeros(len(centers))]-data['translation'])@data['heading_rotation']
    boat_offsets = offsets@data['heading_rotation'][:2, :2]
    train = data['train']
    directions = np.c_[np.cos(train.theta), np.sin(train.theta)]
    radius = width/np.sqrt(2)
    lo = float(np.min(boat_centers[:, :2]@directions.T)-radius)
    hi = float(np.max(boat_centers[:, :2]@directions.T)+radius)
    requested_h = float(config['profile_grid_spacing_m'])
    s_grid = np.linspace(lo, hi, int(np.ceil((hi-lo)/requested_h))+1)
    h = float(np.max(np.diff(s_grid)))
    eta = np.zeros((751, nodes*nodes))
    active = np.flatnonzero(train.data['carrier_amplitudes_m'])
    batch = int(config['profile_component_batch'])
    for first in range(0, len(active), batch):
        check_deadline(started, budget)
        ids = active[first:first+batch]
        profiles = train.profiles(ids, s_grid, times)
        assert profiles.shape == (len(ids), len(times), len(s_grid)) and np.isfinite(profiles).all()
        for jj, component in enumerate(ids):
            projected_offset = boat_offsets@directions[component]
            projected_center = boat_centers[:, :2]@directions[component]
            for frame in range(len(times)):
                projected = projected_offset+projected_center[frame]
                assert projected.min() >= s_grid[0] and projected.max() <= s_grid[-1]
                eta[frame] += np.interp(projected, s_grid, profiles[jj, frame])
        check_deadline(started, budget)
    # A global analytic linear-interpolation error bound. This is a display
    # approximation report, not a new physical/sea-state acceptance threshold.
    curvature = float(np.sum(abs(train.C)*train.k[None, None, :]**2))
    interpolation_bound = h*h*curvature/8
    direct_error, rotation_error = 0., 0.
    sample_ids = np.array([0, nodes*nodes//2, nodes*nodes-1])
    for frame in (0, 375, 750):
        global_points = centers[frame]+offsets[sample_ids]
        boat_points = (np.c_[global_points, np.zeros(len(sample_ids))]-data['translation'])@data['heading_rotation']
        direct = train.elevation(boat_points[:, :2], [times[frame]])[0]
        direct_error = max(direct_error, float(np.max(abs(direct-eta[frame, sample_ids]))))
        # Separate scalar phase expansion verifies world->meanboat mapping.
        for point, expected in zip(global_points, direct):
            boat_point = data['heading_rotation'][:2, :2].T@(point-data['translation'][:2])
            actual = 0j
            for component in active:
                s = float(directions[component]@boat_point)
                phase = np.exp(-1j*train.omega*times[frame])
                actual += np.sum(phase*(train.C[component, 0]*np.exp(1j*train.k*s)+
                                       train.C[component, 1]*np.exp(-1j*train.k*s)))
            rotation_error = max(rotation_error, abs(float(actual.real)-float(expected)))
    assert direct_error <= interpolation_bound+1e-9
    assert rotation_error < 3e-11
    lower, upper = config['z_bounds_m']
    assert hull[:, :, 2].min() > lower and hull[:, :, 2].max() < upper, 'Hull z clipping'
    assert eta.min() > lower and eta.max() < upper, 'Incident surface z clipping'
    horizontal_extent = float(np.max(abs(hull[:, :, :2]-centers[:, None, :])))
    assert horizontal_extent < width/2, 'Hull x/y clipping'
    return dict(indices=indices, time=times, hull=hull, center=centers, dx=dx, dy=dy,
        eta=eta.reshape(751, nodes, nodes),
        check=dict(hull_z_extrema_m=[float(hull[:, :, 2].min()), float(hull[:, :, 2].max())],
            incident_grid_z_extrema_m=[float(eta.min()), float(eta.max())],
            hull_max_horizontal_distance_from_camera_m=horizontal_extent,
            profile_grid_points=len(s_grid), profile_grid_spacing_m=h,
            profile_interval_m=[lo, hi], profile_component_batch=batch,
            analytic_surface_interpolation_error_bound_m=interpolation_bound,
            direct_spectral_vs_cached_surface_max_error_m=direct_error,
            world_meanboat_scalar_phase_max_error_m=rotation_error,
            coverage='all751displayed native-derived frames;65x65surface; allphysical nonzero parent components and both signed branches'))


def save_cache(path, frames):
    np.savez_compressed(path, **{key:frames[key] for key in ('indices', 'time', 'hull', 'center', 'dx', 'dy', 'eta')})


def load_reviewed_cache(review_path, identity):
    review = read(review_path)
    assert review['decision'] == 'preview_pass' and review['reviewed_by']
    preflight_path, cache_path = (resolved(review[key]) for key in ('preflight_json', 'surface_cache_npz'))
    assert sha(preflight_path) == review['preflight_sha256']
    assert sha(cache_path) == review['surface_cache_sha256']
    preflight = read(preflight_path)
    assert preflight['complete'] and preflight['mode'] == 'preview' and preflight['identity'] == identity
    assert preflight['surface_cache_sha256'] == review['surface_cache_sha256']
    with np.load(cache_path) as stored:
        frames = {key:stored[key].copy() for key in stored.files}
    assert frames['eta'].shape == (751, 65, 65) and np.isfinite(frames['eta']).all()
    assert np.array_equal(frames['indices'], np.arange(0, 3001, 4))
    frames['check'] = preflight['preflight']
    bound = {str(review_path.resolve()):sha(review_path), str(preflight_path):sha(preflight_path), str(cache_path):sha(cache_path)}
    return frames, bound


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mode', choices=('preview', 'encode'), default='preview')
    parser.add_argument('--preview-review', type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = read(config_path)
    assert config['display_fps'] == 12.5 and config['native_sample_stride'] == 4 and config['display_frames'] == 751
    data = load_inputs(config)
    identity = dict(source_sha256=sha(__file__), config_sha256=sha(config_path), input_sha256=data['identities'])
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError('Use a fresh output path; preserve prior preview/encode attempts')
    out.mkdir(parents=True)
    budget = config['preview_maximum_seconds'] if args.mode == 'preview' else config['encode_maximum_seconds']
    initial = dict(complete=False, mode=args.mode, identity=identity, command=sys.argv,
        dependencies={name:importlib.metadata.version(name) for name in ('numpy', 'scipy', 'matplotlib', 'threadpoolctl')},
        trajectory_status=data['receipt']['trajectory_status'], compute_cpu_threads=1, maximum_elapsed_seconds=budget)
    write(out/'render-preflight.json', initial)
    review_identity = None
    if args.mode == 'preview':
        frames = prepare_frames(data, config, started, budget)
        save_cache(out/'surface-cache.npz', frames)
        cache_sha = sha(out/'surface-cache.npz')
    else:
        assert args.preview_review, 'Independent literal/provenance preview review required'
        frames, review_identity = load_reviewed_cache(args.preview_review.resolve(), identity)
        cache_sha = next(v for k, v in review_identity.items() if k.endswith('surface-cache.npz'))
    motion_plot(data, out)
    fig = plt.figure(figsize=(12.8, 9))
    ax = fig.add_axes([.025, .035, .82, .70], projection='3d', computed_zorder=False)
    width = config['viewport_width_xy_m']
    lower, upper = config['z_bounds_m']
    ax.set(xlabel='Global east (m)', ylabel='Global north (m)', zlabel='Global up (m)', zlim=[lower, upper])
    ax.set_box_aspect((width, width, upper-lower))
    ax.set_proj_type('ortho')
    yaw = float(data['heading']['mean_box_yaw_global_rad'])
    azimuth = float(np.rad2deg(yaw)-135)
    ax.view_init(elev=35.26438968, azim=azimuth)
    ax.grid(True, alpha=.2)
    fig.suptitle(config.get('title','Native WEC-Sim - incoming buoy sea'), fontsize=16, y=.985)
    fig.text(.5, .947, 'Prescribed incident surface; physical motion and elevation in metres', ha='center', fontsize=11)
    fig.text(.5, .918, 'Equal actual metres; camera follows hull XY; no motion or wave-height amplification', ha='center', fontsize=10)
    fig.text(.5, .890, f'Stationary source Hm0 {config["source_hm0_m"]:.3f} m; finite train differs; mean bow {np.rad2deg(yaw)%360:.3f} deg CCW from east', ha='center', fontsize=10)
    stamp = fig.text(.07, .852, '', fontsize=12)
    telemetry = fig.text(.07, .787, '', fontsize=10, linespacing=1.4)
    color_extent = float(config['surface_color_extent_m'])
    norm = colors.Normalize(vmin=-color_extent, vmax=color_extent)
    colorbar = fig.add_axes([.91, .20, .02, .45])
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='viridis'), cax=colorbar, label='Prescribed incident elevation (m)')
    hull_artist = Poly3DCollection(frames['hull'][0, data['faces']], facecolor='#e9e5d8', edgecolor='#152d3a', linewidth=1., zorder=5)
    ax.add_collection3d(hull_artist)
    surface = None
    def draw(index):
        nonlocal surface
        check_deadline(started, budget)
        center = frames['center'][index]
        ax.set_xlim(center[0]-width/2, center[0]+width/2)
        ax.set_ylim(center[1]-width/2, center[1]+width/2)
        if surface is not None:
            surface.remove()
        elevation = frames['eta'][index]
        surface = ax.plot_surface(frames['dx']+center[0], frames['dy']+center[1], elevation,
            facecolors=cm.viridis(norm(elevation)), rstride=1, cstride=1, linewidth=0,
            antialiased=False, shade=False, alpha=.40, zorder=1)
        hull_artist.set_verts(frames['hull'][index, data['faces']])
        native_index = frames['indices'][index]
        q = data['q'][native_index]
        angles = np.rad2deg(q[3:])
        stamp.set_text(f't = {frames["time"][index]:5.2f} s')
        telemetry.set_text(f'Surge {q[0]:+7.3f} m | Sway {q[1]:+7.3f} m | Heave {q[2]:+7.3f} m\n'
            f'Roll {angles[0]:+7.3f} deg | Pitch {angles[1]:+7.3f} deg | Yaw perturbation {angles[2]:+7.3f} deg')
    for index in (0, 375, 750):
        draw(index)
        fig.savefig(out/f'frame-{index:04d}.png', dpi=100)
    if args.mode == 'encode':
        writer = animation.FFMpegWriter(fps=12.5, codec='libx264', bitrate=3500, extra_args=['-pix_fmt', 'yuv420p', '-threads', '1'])
        with writer.saving(fig, str(out/'reference-isometric.mp4'), dpi=100):
            encoder_pid = writer._proc.pid
            print(dict(encoder_pid=encoder_pid, parent_pid=os.getpid()), flush=True)
            for index in range(751):
                draw(index)
                writer.grab_frame()
        assert writer._proc.returncode == 0
        print(dict(encoder_pid=encoder_pid, exact_process_terminal=True, exit_code=0), flush=True)
    plt.close(fig)
    assert all(sha(resolved(path)) == digest for path, digest in data['identities'].items())
    assert sha(config_path) == identity['config_sha256']
    if review_identity:
        assert all(sha(path) == digest for path, digest in review_identity.items())
    final = dict(initial, complete=True, elapsed_seconds=time.monotonic()-started,
        preflight=frames['check'], surface_cache_sha256=cache_sha,
        native_motion_seconds=60, display_frames=751, native_sample_stride=4, fps=12.5,
        encoded_duration_seconds=60.08 if args.mode == 'encode' else None,
        endpoint_note='Native t=60 endpoint held for one12.5fps frame; no temporal interpolation',
        camera=dict(projection='orthographic', elevation_deg=35.26438968, azimuth_deg=azimuth,
            follows='actual full-hull midpoint XY; no vertical/yaw following', width_xy_m=width,
            z_bounds_m=[lower, upper], azimuth_basis='FIXED mean global bow minus135deg'),
        mapping='equilibrium full hull relative meanCG once; nativeRzRyRx+actualCG; globalR/t once',
        field='exact accepted signed-kappa dispersive train; batchedCZT profiles plus bounded spatial display interpolation; no forcing from visual cache',
        preview_review_identity=review_identity,
        color_saturated_grid_points=int(np.count_nonzero(abs(frames['eta']) > color_extent)),
        previews_seconds=[0, 30, 60], physical_acceptance='Rendering does not establish mesh/time/nonlinear-motion accuracy')
    write(out/'render-preflight.json', final)
    write(out/'output-sha256.json', {path.name:sha(path) for path in out.iterdir() if path.is_file()})
    print(dict(complete=True, mode=args.mode, elapsed_seconds=final['elapsed_seconds'], output=str(out)), flush=True)


if __name__ == '__main__':
    if os.name == 'nt':
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.GetCurrentProcess.restype = ctypes.c_void_p
        kernel.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
        kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        handle = kernel.GetCurrentProcess()
        allowed, system = ctypes.c_size_t(), ctypes.c_size_t()
        assert kernel.GetProcessAffinityMask(handle, ctypes.byref(allowed), ctypes.byref(system))
        assert kernel.SetProcessAffinityMask(handle, allowed.value & -allowed.value)
    with threadpool_limits(1):
        main()

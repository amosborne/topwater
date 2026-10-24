"""Render a bound, reviewed native trajectory and its unchanged buoy sea.

Source preparation does not authorize execution. Preview and encoding are
separate serial artifact jobs; all input identities must be bound first.
"""
import os
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import argparse
import ctypes
import hashlib
import importlib.metadata
import json
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

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def resolved(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT/path


def check_deadline(started, seconds):
    if time.monotonic()-started > seconds:
        raise TimeoutError('Declared render budget exhausted; partial artifacts are retained')


def sea_arrays(bank, gravity):
    a, omega, theta, phase = (np.asarray(bank[key], dtype=float) for key in
        ('amplitudes', 'angular_frequencies', 'directions', 'phases'))
    assert all(v.shape == (384,) and np.isfinite(v).all() for v in (a, omega, theta, phase))
    assert np.all(a >= 0) and np.all(omega > 0)
    k = omega[:, None]**2/gravity*np.column_stack((np.cos(theta), np.sin(theta)))
    return a, omega, phase, k


def eta_direct(bank, gravity, points, tt):
    a, omega, phase, k = sea_arrays(bank, gravity)
    return a@np.cos(k@np.asarray(points).T+phase[:, None]-omega[:, None]*tt)


def load_inputs(config):
    if config['status'] != 'ready_for_coordinator_render_allocation':
        raise ValueError('Source-only configuration: accepted motion has not been bound')
    required = ('case', 'heading_file', 'global_bank', 'motion_mat', 'summary_json',
                'native_output_mat', 'motion_acceptance_receipt')
    identities = config['input_sha256']
    for key in required:
        assert config[key] and config[key] in identities, f'Unbound required input: {key}'
    assert Path(__file__).resolve().relative_to(ROOT).as_posix() in identities
    for path, digest in identities.items():
        assert digest and sha(resolved(path)) == digest, f'Input identity differs: {path}'
    motion_sha = identities[config['motion_mat']]
    assert motion_sha not in config['rejected_motion_sha256'], 'Known rejected trajectory'
    receipt = read(resolved(config['motion_acceptance_receipt']))
    assert receipt['decision'] == 'accepted_for_visualization'
    assert receipt['motion_sha256'] == motion_sha
    assert receipt['summary_sha256'] == identities[config['summary_json']]
    assert receipt['native_output_sha256'] == identities[config['native_output_mat']]
    assert receipt['heading_sha256'] == identities[config['heading_file']]
    assert receipt['trajectory_status'] in ('provisional_model_trajectory', 'accepted_reference_case')
    assert isinstance(receipt['reviewed_by'], str) and receipt['reviewed_by']
    case, heading, bank, summary = (read(resolved(config[key])) for key in
        ('case', 'heading_file', 'global_bank', 'summary_json'))
    assert summary['complete'] and summary['samples'] == 3001 and summary['duration_s'] == 60
    assert 'reject' not in str(receipt['trajectory_status']).lower()
    assert heading['observation_id'] == bank['provenance']['observation_id'] == receipt['observation_id']
    assert heading['physical_input_sha256'] == identities[config['case']]
    motion = loadmat(resolved(config['motion_mat']), simplify_cells=True)
    raw = loadmat(resolved(config['native_output_mat']), simplify_cells=True)['native_body']
    t = np.asarray(motion['time']).ravel()
    q, position = np.asarray(motion['q']), np.asarray(motion['position_world'])
    assert t.shape == (3001,) and q.shape == position.shape == (3001, 6)
    assert np.isfinite(q).all() and np.isfinite(position).all() and np.isfinite(t).all()
    np.testing.assert_allclose(t, np.arange(3001)*.02, rtol=0, atol=1e-9)
    assert np.array_equal(np.asarray(raw['time']).ravel(), t)
    assert np.asarray(raw['signals']['values']).shape == (3001, 60)
    assert np.array_equal(np.asarray(raw['signals']['values'])[:, :6], position)
    cg = np.asarray(case['world_cg_m'])
    np.testing.assert_allclose(np.asarray(motion['mean_world_cg']).ravel(), cg, rtol=0, atol=2e-12)
    np.testing.assert_allclose(position[:, :3]-cg, q[:, :3], rtol=0, atol=2e-12)
    np.testing.assert_allclose(position[:, 3:], q[:, 3:], rtol=0, atol=2e-12)
    hull = case['minimum_hull']
    eq = hull['equilibrium_quaternion_wxyz']
    equilibrium = Rotation.from_quat([*eq[1:], eq[0]]).as_matrix()
    mean_vertices = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@equilibrium.T+cg
    faces = np.asarray([f for i, f in enumerate(hull['triangular_faces'])
                        if i not in hull['nonphysical_closure_faces']], dtype=int)
    heading_rotation = np.asarray(heading['global_from_boat_rotation'])
    translation = np.asarray(heading['global_from_boat_translation_m'])
    yaw = float(heading['mean_box_yaw_global_rad'])
    np.testing.assert_allclose(heading_rotation, Rotation.from_euler('z', yaw).as_matrix(), rtol=0, atol=2e-14)
    assert translation.shape == (3,) and translation[2] == 0
    a, omega, phase, k = sea_arrays(bank['waves'], case['water']['gravity_m_s2'])
    ab, wb, pb, kb = sea_arrays(heading['boat_frame_waves'], case['water']['gravity_m_s2'])
    assert np.array_equal(a, ab) and np.array_equal(omega, wb)
    np.testing.assert_allclose(kb, k@heading_rotation[:2, :2], rtol=0, atol=3e-15)
    np.testing.assert_allclose(np.exp(1j*pb), np.exp(1j*(phase+k@translation[:2])), rtol=0, atol=3e-14)
    assert abs(summary['initial_global_heading_rad']-yaw) < 1e-12
    return dict(case=case, heading=heading, bank=bank, summary=summary, receipt=receipt,
        time=t, q=q, position=position, cg=cg, local=mean_vertices-cg, faces=faces,
        heading_rotation=heading_rotation, translation=translation, a=a, omega=omega,
        phase=phase, k=k, identities=identities)


def prepare_frames(data, config, started, budget):
    # Every fourth native sample: no temporal interpolation or sample rebasing.
    indices = np.arange(0, 3001, 4)
    times = data['time'][indices]
    np.testing.assert_allclose(times, np.arange(751)/12.5, rtol=0, atol=1e-9)
    r_native = Rotation.from_euler('xyz', data['position'][indices, 3:]).as_matrix()
    posed_native = np.einsum('fij,vj->fvi', r_native, data['local'])+data['position'][indices, None, :3]
    posed_global = posed_native@data['heading_rotation'].T+data['translation']
    # Follow the geometric hull midpoint horizontally; global z remains fixed,
    # preserving visible heave. Camera azimuth never follows instantaneous yaw.
    centers = (posed_global[:, :, :2].min(axis=1)+posed_global[:, :, :2].max(axis=1))/2
    width, nodes = float(config['viewport_width_xy_m']), int(config['surface_grid_nodes'])
    assert width == 12 and nodes == 65
    dx, dy = np.meshgrid(np.linspace(-width/2, width/2, nodes), np.linspace(-width/2, width/2, nodes))
    offsets = np.c_[dx.ravel(), dy.ravel()]
    spatial = data['k']@offsets.T
    cosine, sine = np.cos(spatial), np.sin(spatial)
    eta = np.empty((751, nodes*nodes))
    for first in range(0, len(times), 64):
        check_deadline(started, budget)
        sl = slice(first, first+64)
        argument = centers[sl]@data['k'].T+data['phase']-times[sl, None]*data['omega']
        eta[sl] = (np.cos(argument)*data['a'])@cosine-(np.sin(argument)*data['a'])@sine
    lower, upper = config['z_bounds_m']
    assert lower == -3 and upper == 3
    assert np.min(posed_global[:, :, 2]) > lower and np.max(posed_global[:, :, 2]) < upper, 'Hull z clipping'
    assert np.min(eta) > lower and np.max(eta) < upper, 'Incident surface z clipping'
    relative_xy = posed_global[:, :, :2]-centers[:, None, :]
    assert np.max(abs(relative_xy)) < width/2, 'Hull x/y clipping'
    scalar_error, closure_error = 0., 0.
    gravity = data['case']['water']['gravity_m_s2']
    for frame in (0, 375, 750):
        points = offsets[[0, 2112, 4224]]+centers[frame]
        # Independent scalar all-component sum checks optimized spatial factoring.
        for point, value in zip(points, eta[frame, [0, 2112, 4224]]):
            direct = sum(float(a)*np.cos(float(k@point)+float(phi)-float(w*times[frame]))
                         for a, k, phi, w in zip(data['a'], data['k'], data['phase'], data['omega']))
            scalar_error = max(scalar_error, abs(float(value)-direct))
        native_points = (np.c_[points, np.zeros(3)]-data['translation'])@data['heading_rotation']
        local_eta = eta_direct(data['heading']['boat_frame_waves'], gravity, native_points[:, :2], times[frame])
        global_eta = eta_direct(data['bank']['waves'], gravity, points, times[frame])
        closure_error = max(closure_error, float(np.max(abs(local_eta-global_eta))))
    assert scalar_error < 3e-12 and closure_error < 3e-12
    return dict(indices=indices, time=times, hull=posed_global, center=centers,
        dx=dx, dy=dy, eta=eta.reshape(751, nodes, nodes),
        check=dict(hull_z_extrema_m=[float(posed_global[:, :, 2].min()), float(posed_global[:, :, 2].max())],
            incident_grid_z_extrema_m=[float(eta.min()), float(eta.max())],
            hull_max_horizontal_distance_from_camera_m=float(np.max(abs(relative_xy))),
            scalar_vs_factored_eta_max_error_m=float(scalar_error),
            global_native_eta_closure_max_error_m=float(closure_error),
            coverage='all751displayed native-derived frames; surface grid65x65, not continuous-region bounds'))


def motion_plot(data, out):
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True, layout='constrained')
    labels = ('Surge (m)', 'Sway (m)', 'Heave (m)', 'Roll (deg)', 'Pitch (deg)', 'Yaw perturbation (deg)')
    for j, ax in enumerate(axes.flat):
        values = data['q'][:, j] if j < 3 else np.rad2deg(data['q'][:, j])
        ax.plot(data['time'], values, color='#075c80', linewidth=1.1)
        ax.set(title=labels[j], xlim=(0, 60))
        ax.grid(True, alpha=.3)
    for ax in axes[-1]:
        ax.set_xlabel('Time from declared initial state (s)')
    fig.suptitle('Native WEC-Sim: complete six motion histories; known neutral pose retained')
    fig.savefig(out/'six-motions.png', dpi=180)
    plt.close(fig)


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
    preview_review_identity = None
    if args.mode == 'encode':
        assert args.preview_review, 'Independent preview review receipt required before encoding'
        review = read(args.preview_review)
        assert review['decision'] == 'preview_pass' and review['reviewed_by']
        preflight_path = resolved(review['preflight_json'])
        assert sha(preflight_path) == review['preflight_sha256']
        preview = read(preflight_path)
        assert preview['identity'] == identity and preview['mode'] == 'preview' and preview['complete']
        preview_review_identity = {str(args.preview_review.resolve()): sha(args.preview_review),
                                   str(preflight_path): sha(preflight_path)}
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError('Use a fresh output directory; preserve all prior renders and previews')
    out.mkdir(parents=True)
    budget = config['preview_maximum_seconds'] if args.mode == 'preview' else config['encode_maximum_seconds']
    initial = dict(complete=False, mode=args.mode, identity=identity, command=sys.argv,
        dependencies={name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'matplotlib', 'threadpoolctl')},
        trajectory_status=data['receipt']['trajectory_status'], compute_cpu_threads=1,
        maximum_elapsed_seconds=budget)
    write(out/'render-preflight.json', initial)
    frames = prepare_frames(data, config, started, budget)
    motion_plot(data, out)
    fig = plt.figure(figsize=(12.8, 9))
    ax = fig.add_axes([.025, .04, .80, .69], projection='3d', computed_zorder=False)
    ax.set(xlabel='Global east (m)', ylabel='Global north (m)', zlabel='Global up (m)', zlim=config['z_bounds_m'])
    ax.set_box_aspect((12, 12, 6))
    ax.set_proj_type('ortho')
    yaw = data['heading']['mean_box_yaw_global_rad']
    azimuth = np.rad2deg(yaw)-135
    ax.view_init(elev=35.26438968, azim=azimuth)
    ax.grid(True, alpha=.2)
    fig.suptitle('Native WEC-Sim — collected-buoy sea', fontsize=16, y=.985)
    qualifier = 'Provisional model trajectory' if data['receipt']['trajectory_status'] == 'provisional_model_trajectory' else 'Accepted reference case'
    fig.text(.5, .947, f'{qualifier}; prescribed incident surface only', ha='center', fontsize=11)
    fig.text(.5, .918, 'Equal actual metres; no motion or wave-height amplification; camera follows hull XY', ha='center', fontsize=10)
    fig.text(.5, .890, f'Mean bow {np.rad2deg(yaw)%360:.3f}° CCW from east; heading offset {data["heading"]["relative_heading_offset_deg"]:+g}°', ha='center', fontsize=10)
    stamp = fig.text(.07, .852, '', fontsize=12)
    telemetry = fig.text(.07, .787, '', fontsize=10, linespacing=1.4)
    color_extent = 3*np.sqrt(np.sum(data['a']**2)/2)
    assert color_extent > 0
    norm = colors.Normalize(vmin=-color_extent, vmax=color_extent)
    colorbar = fig.add_axes([.89, .20, .02, .46])
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='viridis'), cax=colorbar, label='Incident elevation (m); colors ±3σ')
    hull_artist = Poly3DCollection(frames['hull'][0, data['faces']], facecolor='#e9e5d8',
        edgecolor='#152d3a', linewidth=1.0, zorder=5)
    ax.add_collection3d(hull_artist)
    surface = None
    def draw(index):
        nonlocal surface
        check_deadline(started, budget)
        center = frames['center'][index]
        ax.set_xlim(center[0]-6, center[0]+6)
        ax.set_ylim(center[1]-6, center[1]+6)
        x, y = frames['dx']+center[0], frames['dy']+center[1]
        if surface is not None:
            surface.remove()
        elevation = frames['eta'][index]
        surface = ax.plot_surface(x, y, elevation, facecolors=cm.viridis(norm(elevation)),
            rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False, alpha=.40, zorder=1)
        hull_artist.set_verts(frames['hull'][index, data['faces']])
        native = frames['indices'][index]
        q = data['q'][native]
        angles = np.rad2deg(q[3:])
        cg_global = data['position'][native, :3]@data['heading_rotation'].T+data['translation']
        stamp.set_text(f't = {frames["time"][index]:5.2f} s')
        telemetry.set_text(f'Surge {q[0]:+7.3f} m | Sway {q[1]:+7.3f} m | Heave {q[2]:+7.3f} m\n'
            f'Roll {angles[0]:+7.3f}° | Pitch {angles[1]:+7.3f}° | Yaw perturbation {angles[2]:+7.3f}°\n'
            f'Global CG: E {cg_global[0]:+.3f} m, N {cg_global[1]:+.3f} m, Z {cg_global[2]:+.3f} m')
    for index in (0, 375, 750):
        draw(index)
        fig.savefig(out/f'frame-{index:04d}.png', dpi=100)
    if args.mode == 'encode':
        writer = animation.FFMpegWriter(fps=12.5, codec='libx264', bitrate=3500,
            extra_args=['-pix_fmt', 'yuv420p', '-threads', '1'])
        with writer.saving(fig, str(out/'reference-isometric.mp4'), dpi=100):
            encoder_pid = writer._proc.pid
            print(json.dumps(dict(encoder_pid=encoder_pid, parent_pid=os.getpid())), flush=True)
            for index in range(751):
                draw(index)
                writer.grab_frame()
        assert writer._proc.returncode == 0
        print(json.dumps(dict(encoder_pid=encoder_pid, exact_process_terminal=True, exit_code=0)), flush=True)
    plt.close(fig)
    # Producers/data and the separately reviewed preview stay immutable.
    assert all(sha(resolved(path)) == digest for path, digest in data['identities'].items())
    assert sha(config_path) == identity['config_sha256']
    if preview_review_identity:
        assert all(sha(path) == digest for path, digest in preview_review_identity.items())
    final = dict(initial, complete=True, elapsed_seconds=time.monotonic()-started,
        preflight=frames['check'], native_motion_seconds=60, display_frames=751,
        native_sample_stride=4, fps=12.5, encoded_duration_seconds=60.08 if args.mode == 'encode' else None,
        endpoint_note='Native t=60 endpoint held for one12.5fps frame; no temporal interpolation',
        camera=dict(projection='orthographic', elevation_deg=35.26438968, azimuth_deg=float(azimuth),
            follows='actual physical full-hull midpoint XY; no vertical/yaw following', width_xy_m=12, z_bounds_m=[-3, 3]),
        mapping='equilibrium full hull minus known meanCG once; nativeRzRyRx + actualCG; globalheadingR/t once',
        field='all384unchanged global harmonics at actual global coordinates, Re(a*exp(i(k.x+phase-omega*t)))',
        color_scale='±3standard deviations only affects colors; actual elevations are unmodified',
        preview_review_identity=preview_review_identity,
        color_saturated_grid_points=int(np.count_nonzero(abs(frames['eta']) > color_extent)),
        previews_seconds=[0, 30, 60], physical_acceptance='rendering does not establish mesh/geometry/linearity accuracy')
    write(out/'render-preflight.json', final)
    outputs = {path.name: sha(path) for path in out.iterdir() if path.is_file()}
    write(out/'output-sha256.json', outputs)
    print(json.dumps(dict(complete=True, mode=args.mode, elapsed_seconds=final['elapsed_seconds'], output=str(out))), flush=True)


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

"""Render native motion and the exact accepted prescribed incident field."""
import argparse
import ctypes
import json
import os
from pathlib import Path
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


def main():
    started = time.monotonic()
    p = argparse.ArgumentParser()
    p.add_argument('--case', type=Path, required=True)
    p.add_argument('--output', type=Path)
    p.add_argument('--preview-only', action='store_true')
    args = p.parse_args()
    out = args.output or args.case
    out.mkdir(parents=True, exist_ok=True)
    assert(not (out/'reference-isometric.mp4').exists())
    m = loadmat(args.case/'motion.mat', simplify_cells=True)
    t, q = np.asarray(m['time']), np.asarray(m['q'])
    pos = np.asarray(m['position_world'])
    assert t[0] == 0 and abs(t[-1]-10) < 1e-9
    assert np.max(abs(q[0])) < 1e-9
    packet = np.load('artifacts/reference/minimum-20261024/packet-attempt3/packet.npz')
    boat = json.loads(Path('reference/capytaine/canonical-boat.json').read_text())
    hull = boat['minimum_hull']
    eq = hull['equilibrium_quaternion_wxyz']
    R0 = Rotation.from_quat([*eq[1:], eq[0]]).as_matrix()
    cg = np.asarray(boat['world_cg_m'])
    vertices = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@R0.T+cg
    faces = np.array([f for i, f in enumerate(hull['triangular_faces']) if i not in hull['nonphysical_closure_faces']])
    local = vertices-cg
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), sharex=True, layout='constrained')
    labels = ['Surge (m)', 'Sway (m)', 'Heave (m)', 'Roll (deg)', 'Pitch (deg)', 'Yaw (deg)']
    for i, ax in enumerate(axes.flat):
        y = q[:, i] if i < 3 else np.rad2deg(q[:, i])
        ax.plot(t, y, color='#075c80')
        ax.set_ylabel(labels[i])
        ax.grid(True, alpha=.3)
        ax.set_xlim(0, 10)
    axes[-1, 0].set_xlabel('Time from actual rest (s)')
    axes[-1, 1].set_xlabel('Time from actual rest (s)')
    fig.suptitle('Minimum-mesh WEC-Sim: native six motion histories')
    fig.savefig(out/'six-motions.png', dpi=180)
    plt.close(fig)
    view_center = (vertices[:, :2].min(axis=0)+vertices[:, :2].max(axis=0))/2
    view_width = 12.
    xlim = (view_center[0]-view_width/2, view_center[0]+view_width/2)
    ylim = (view_center[1]-view_width/2, view_center[1]+view_width/2)
    x, y = np.meshgrid(np.linspace(*xlim, 65), np.linspace(*ylim, 65))
    xy = np.c_[x.ravel(), y.ravel()]
    theta, k, omega, c = packet['direction'], packet['k'], packet['w'], packet['c']
    spatial = [np.exp(1j*np.outer(xy@np.array([np.cos(th), np.sin(th)]), k[j])) for j, th in enumerate(theta)]
    norm = colors.Normalize(vmin=-sum(boat['sea']['amplitudes_m']), vmax=sum(boat['sea']['amplitudes_m']))
    fig = plt.figure(figsize=(12.8, 9))
    ax = fig.add_axes([.03, .04, .79, .80], projection='3d', computed_zorder=False)
    ax.set(xlim=xlim, ylim=ylim, zlim=(-2, 2), xlabel='World x (m)', ylabel='World y (m)', zlabel='World z (m)')
    ax.set_box_aspect((view_width, view_width, 4))
    ax.set_proj_type('ortho')
    ax.view_init(elev=35.26438968, azim=-135)
    ax.grid(True, alpha=.2)
    fig.suptitle('Freely floating 24 ft boat - minimum solver mesh', fontsize=16, y=.98)
    fig.text(.5, .935, 'Transparent prescribed incident surface; no computed scattered/radiated field', ha='center', fontsize=11)
    fig.text(.5, .905, 'Equal physical spatial units; no motion or height exaggeration', ha='center', fontsize=11)
    stamp = fig.text(.07, .855, '', fontsize=13)
    telemetry = fig.text(.07, .805, '', fontsize=11, linespacing=1.5)
    hull_artist = Poly3DCollection(vertices[faces], facecolor='#e9e5d8', edgecolor='#152d3a', linewidth=1.1, zorder=5)
    ax.add_collection3d(hull_artist)
    bar = fig.add_axes([.89, .24, .02, .50])
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap='viridis'), cax=bar, label='Prescribed incident elevation (m)')
    surface = None
    fps = 25
    times = np.arange(10*fps+1)/fps
    def draw_frame(index, tt):
        nonlocal surface
        eta = sum((spatial[j]@(c[j]*np.exp(-1j*omega[j]*tt))).real for j in range(3)).reshape(x.shape)
        if surface is not None:
            surface.remove()
        surface = ax.plot_surface(x, y, eta, facecolors=cm.viridis(norm(eta)), rstride=1, cstride=1,
                                  linewidth=0, antialiased=False, shade=False, alpha=.65, zorder=1)
        angles = np.array([np.interp(tt, t, pos[:, j]) for j in range(3, 6)])
        world_cg = np.array([np.interp(tt, t, pos[:, j]) for j in range(3)])
        # Native quaternion2EulXYZ is the standard roll/pitch/yaw inverse
        # of Rz(yaw)*Ry(pitch)*Rx(roll), scipy's extrinsic lower-case xyz.
        R = Rotation.from_euler('xyz', angles).as_matrix()
        moved = local@R.T+world_cg
        hull_artist.set_verts(moved[faces])
        stamp.set_text(f't = {tt:5.2f} s')
        displacement_mm = 1000*(world_cg-cg)
        rotation_deg = np.rad2deg(angles)
        telemetry.set_text(
            f'Surge {displacement_mm[0]:+7.3f} mm | Sway {displacement_mm[1]:+7.3f} mm | Heave {displacement_mm[2]:+7.3f} mm\n'
            f'Roll  {rotation_deg[0]:+8.4f} deg | Pitch {rotation_deg[1]:+8.4f} deg | Yaw {rotation_deg[2]:+8.4f} deg')
        if index in [0, 125, 250]:
            fig.savefig(out/f'frame-{index:03d}.png', dpi=100)
    if args.preview_only:
        for index in [0, 125, 250]:
            draw_frame(index, times[index])
        plt.close(fig)
        report = dict(preview_only=True, timestamps_s=[0, 5, 10], elapsed_s=time.monotonic()-started)
        (out/'preview-summary.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(report), flush=True)
        return
    writer = animation.FFMpegWriter(fps=fps, codec='libx264', bitrate=3500,
                                   extra_args=['-pix_fmt', 'yuv420p', '-threads', '1'])
    with writer.saving(fig, str(out/'reference-isometric.mp4'), dpi=100):
        for index, tt in enumerate(times):
            draw_frame(index, tt)
            writer.grab_frame()
    plt.close(fig)
    result = dict(native_motion_seconds=float(t[-1]-t[0]), encoded_frames=len(times), fps=fps,
                  elapsed_s=time.monotonic()-started, compute_cpu_affinity='one allowed logical CPU; child encoder inherits',
                  encoded_duration_s=len(times)/fps, first_timestamp_s=0, last_timestamp_s=10,
                  endpoint_note='Last native t=10 sample is held for one encoded frame',
                  camera=dict(elevation_deg=35.26438968, azimuth_deg=-135, projection='orthographic', fixed=True),
                  viewport=dict(center_xy_m=view_center.tolist(), width_xy_m=view_width, z_bounds_m=[-2, 2]),
                  artist_visibility='hull overlay on transparent incident-only surface; alpha0.65',
                  spatial_units='equal actual metres in x/y/z; no exaggerated motion or incident elevation',
                  motion_mapping='mean-oriented native hull centered once at mean CG; native RzRyRx about instantaneous world CG',
                  prescribed_surface='exact frozen packet.npz complex spectral coefficients; incident only')
    (out/'render-summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


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

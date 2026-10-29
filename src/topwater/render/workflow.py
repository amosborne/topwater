"""Shared saved-geometry video orchestration for native CPU and hardware GPU."""
from pathlib import Path
import sys
import tempfile
import time
import numpy as np
from .scene import RenderSettings
from ..workflow_data import Result
from .mesh import mesh_edges
from .backend import DrawingBackend
from .hud import HUD,water_colors


def render_result(result,output_path,settings=None,*,device='auto'):
    started=time.monotonic()
    from . import geometry
    from .frames import VideoFrames
    if isinstance(result,Result):
        frames=VideoFrames.from_result(result,RenderSettings() if settings is None else settings)
    elif isinstance(result,VideoFrames):
        frames=result
    else:
        raise TypeError('render_result requires a validated workflow_data.Result')
    size,stamps,indices,fps=frames.size,frames.stamps,frames.indices,frames.fps
    count=len(indices)
    lower,upper=frames.lower,frames.upper
    output=Path(output_path)
    if output.suffix.lower()!='.mp4':raise ValueError('output must name a complete .mp4 file')
    if frames.ambiguous:
        print('WARNING render: Euler HUD angles are nonunique at gimbal lock; geometry uses the saved quaternions',file=sys.stderr)
    magnitude=float(np.ceil(max(float(abs(frames.water_vertices[:,2]).max()),1e-12)/.1)*.1)
    edges,rim=mesh_edges(frames.hull_faces)
    reference=geometry._flat_reference_waterline(frames.hull_at(0),frames.hull_faces,lower,upper)
    ticks=geometry._edge_graticule_ticks(lower,upper)
    preparation=time.monotonic()-started
    backend=None
    completed=changes=0
    validation=draw_seconds=send_seconds=close_seconds=verify_seconds=geometry_seconds=0.
    timings=dict(mask_seconds=0.,edge_seconds=0.,contact_seconds=0.)
    previous_faces=previous_xy=previous_count=projected=loops=None
    phase='drawing backend startup'
    phase_started=time.monotonic()
    try:
        backend=DrawingBackend(size,lower,upper,device=device)
        print('render backend: '+str(backend.info),flush=True)
        hud=HUD(size,lower,upper,magnitude,ticks)
        backend.texture('graticule',hud.graticule())
        setup=time.monotonic()-started-preparation
        import imageio_ffmpeg
        from matplotlib import tri
        with tempfile.TemporaryDirectory(prefix='topwater-render-') as temporary:
            encoded=Path(temporary)/'video.mp4'
            writer=imageio_ffmpeg.write_frames(str(encoded),size,fps=fps,codec='libx264',
                pix_fmt_in='rgba',pix_fmt_out='yuv420p',bitrate='3000k',macro_block_size=2,
                ffmpeg_timeout=0,output_params=['-threads','1'])
            phase='encoder startup';phase_started=time.monotonic()
            writer.send(None)
            encoder_start=time.monotonic()-phase_started
            try:
                for index in indices:
                    phase='geometry validation';phase_started=time.monotonic()
                    points,faces=frames.surface_at(index)
                    hull=frames.hull_at(index)
                    changed=previous_count!=len(points) or previous_faces is None or not np.array_equal(previous_faces,faces)
                    if changed:
                        loops=geometry._topology(faces,len(points))
                        if previous_faces is not None:changes+=1
                    fixed=not changed and len(loops)==1 and np.array_equal(previous_xy,points[:,:2])
                    if not fixed:
                        projected=geometry._projected_frame(points,faces,hull,lower,upper,loops,tri)
                        previous_xy=points[:,:2].copy()
                    previous_faces,previous_count=faces,len(points)
                    validation+=time.monotonic()-phase_started
                    phase='wetted display geometry';phase_started=time.monotonic()
                    masks,_,_,contact=geometry._wetted_display(hull,frames.physical_faces_at(index),points,faces,loops,
                        timings,projected,compute_submerged_edges=False)
                    layers=[('reference waterline',reference[:,:,:2],[119/255]*3+[1.],1.39,True),
                        ('boat mesh',hull[edges,:2],[176/255]*3+[1.],.97,False),
                        ('boat outline',hull[rim,:2],[0,0,0,1],2.08,False),
                        ('physical waterline',contact[:,:,:2],[0,0,1,1],1.81,False)]
                    geometry_seconds+=time.monotonic()-phase_started
                    phase='scene drawing';phase_started=time.monotonic()
                    vertices=np.c_[backend.to_clip(points[:,:2]),water_colors(points[:,2],magnitude)]
                    header=frames.header_at(index)
                    backend.texture('header',hud.header(header),
                        rectangle=(0,hud.header_top,size[0],hud.header_top+hud.header_height))
                    rgba=backend.draw(vertices.astype(np.float32),faces,masks,layers)
                    draw_seconds+=time.monotonic()-phase_started
                    phase='encoder send';phase_started=time.monotonic()
                    writer.send(rgba)
                    send_seconds+=time.monotonic()-phase_started
                    completed+=1
                    if completed==1 or completed%50==0 or completed==count:
                        print(f'render progress: {completed}/{count} frames; elapsed={time.monotonic()-started:.3f}s',flush=True)
            finally:
                if completed==count:phase='encoder finalization'
                phase_started=time.monotonic()
                writer.close();close_seconds=time.monotonic()-phase_started
            phase='full video verification';phase_started=time.monotonic()
            geometry._verify_video(encoded,count,fps,size)
            verify_seconds=time.monotonic()-phase_started
            phase='video publication';phase_started=time.monotonic()
            geometry._publish_video(encoded,output)
        return dict(complete=True,output=str(output),frames=count,fps=fps,size=size,
            source_end_s=float(stamps[indices[-1]]),video_duration_s=count/fps,
            preparation_seconds=preparation,setup_seconds=setup,geometry_seconds=geometry_seconds,
            draw_seconds=draw_seconds,encoder_send_seconds=send_seconds,encoder_start_seconds=encoder_start,
            encoder_close_seconds=close_seconds,verification_seconds=verify_seconds,elapsed_s=time.monotonic()-started,
            topology_changes=changes,ambiguous_euler_labels=frames.ambiguous,validation_seconds=validation,
            mask_path_seconds=0.,**timings,backend=backend.info['backend_type'],adapter=backend.info,
            hud=hud.metadata)
    except BaseException:
        print(f'render failed: phase={phase}; completed={completed}/{count} frames; '
            f'elapsed={time.monotonic()-started:.3f}s active_phase_elapsed={time.monotonic()-phase_started:.3f}s '
            f'mask={timings["mask_seconds"]:.3f}s encoder_close={close_seconds:.3f}s',file=sys.stderr,flush=True)
        raise
    finally:
        if backend is not None:backend.close()

"""Adapters for saved rigid poses and the public explicit-mesh scene API."""
from dataclasses import dataclass
import numpy as np
from ..solve.settings import VIDEO_SIZE_PX
from .hud import header_text


@dataclass
class VideoFrames:
    stamps: np.ndarray
    indices: np.ndarray
    fps: float
    size: tuple
    lower: np.ndarray
    upper: np.ndarray
    hull_faces: np.ndarray
    water_vertices: np.ndarray
    hull_at: object
    surface_at: object
    header_at: object
    ambiguous: bool=False
    hull_faces_at: object=None

    def physical_faces_at(self,index):
        return self.hull_faces if self.hull_faces_at is None else self.hull_faces_at(index)

    @classmethod
    def from_result(cls,result,settings):
        from . import geometry
        positions,rotations,angles,ambiguous=geometry.sampled_poses(result)
        local=result.hull_vertices_m-result.cg_m
        def hull_at(i):return positions[i]+local@rotations[i].T
        if result.view_bounds_xy_m is not None:
            lower,upper=result.view_bounds_xy_m
        else:
            if not np.isfinite(settings.framing) or settings.framing<=1:
                raise ValueError('framing must be finite and exceed one')
            low,high=np.full(2,np.inf),np.full(2,-np.inf)
            for i in range(len(positions)):
                xy=hull_at(i)[:,:2]
                low,high=np.minimum(low,xy.min(axis=0)),np.maximum(high,xy.max(axis=0))
            half=float((high-low).max()/2)*settings.framing
            if half<=0:raise ValueError('projected hull has no finite view extent')
            if np.max(abs(np.r_[low,high]))>half:
                raise ValueError('origin-centred camera must cover the full hull trajectory at the fixed zoom')
            lower,upper=-np.full(2,half),np.full(2,half)
        indices,fps=geometry.render_clock(result)
        size=tuple(int(v) for v in (VIDEO_SIZE_PX if result.video_size_px is None else result.video_size_px))
        def header_at(i):return header_text(result.surface_time_s[i],angles[i,0],angles[i,1],positions[i,2]-positions[0,2])
        return cls(result.surface_time_s,indices,fps,size,lower,upper,result.hull_faces,
            result.surface_vertices_m,hull_at,result.surface_frame,header_at,ambiguous)

    @classmethod
    def from_scene(cls,scene,settings):
        from .scene import convex_envelope
        stamp=scene.time_s
        if len(stamp)<2 or not np.allclose(np.diff(stamp),stamp[1]-stamp[0],atol=1e-9,rtol=1e-7):
            raise ValueError('Movie input needs at least two uniformly spaced times; resampling belongs to the producer adapter')
        required=('roll_deg','pitch_deg','heave_m')
        if any(key not in scene.telemetry for key in required):
            raise ValueError('Top-down header requires producer telemetry: roll_deg, pitch_deg, heave_m')
        if not np.isfinite(settings.framing) or settings.framing<=1:
            raise ValueError('Framing must exceed one')
        oriented=[convex_envelope(hull,scene.hull_faces)[1] for hull in scene.hull_points]
        center=((scene.hull_points.min(axis=1)+scene.hull_points.max(axis=1))/2).mean(axis=0)[:2]
        half=float(abs(scene.hull_points[:,:,:2]-center).max())*settings.framing
        def header_at(i):return header_text(stamp[i],*(scene.telemetry[key][i] for key in required))
        return cls(stamp,np.arange(len(stamp)),float(1/(stamp[1]-stamp[0])),VIDEO_SIZE_PX,
            center-half,center+half,oriented[0],scene.water_points.reshape(-1,3),
            lambda i:scene.hull_points[i],lambda i:(scene.water_points[i],scene.water_faces),header_at,
            hull_faces_at=lambda i:oriented[i])

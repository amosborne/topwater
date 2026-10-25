"""Display-only orthographic, depth-buffered presentation of an ACCEPTED cache.
No wave/profile generation, simulation, trajectory interpolation or encoding.
"""
import os
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[variable]='1'
import argparse
from collections import Counter
from pathlib import Path
import time
import subprocess
import shutil
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import ConvexHull
from scipy.ndimage import binary_dilation
from threadpoolctl import threadpool_limits
from render_calm_buoy_reference import load_inputs, load_reviewed_cache
from render_generated_reference import read, write, sha, resolved, check_deadline


class Raster:
    def __init__(self,width,height,center,elevation,azimuth,span):
        self.width=width;self.height=height
        e,a=np.deg2rad([elevation,azimuth])
        self.toward=np.array([np.cos(e)*np.cos(a),np.cos(e)*np.sin(a),np.sin(e)])
        self.right=np.array([-np.sin(a),np.cos(a),0.])
        self.up=np.cross(self.toward,self.right)
        self.center=np.array([*center,0.])  # no vertical following
        self.scale=width/span
        self.hull_depth=np.full((height,width),-np.inf)
        self.water_depth=np.full((height,width),-np.inf)
        self.hull_inside=np.zeros((height,width),dtype=bool)
        self.hull_color=np.zeros((height,width,3));self.water_color=np.zeros_like(self.hull_color)
    def project(self,p):
        d=p-self.center
        return np.c_[self.width/2+self.scale*(d@self.right),self.height*.54-self.scale*(d@self.up),d@self.toward]
    def triangle(self,triangle,color,water=False,interior=False):
        q=self.project(triangle);x0=max(0,int(np.floor(q[:,0].min())));x1=min(self.width-1,int(np.ceil(q[:,0].max())))
        y0=max(0,int(np.floor(q[:,1].min())));y1=min(self.height-1,int(np.ceil(q[:,1].max())))
        if x1<x0 or y1<y0:return
        x,y=np.meshgrid(np.arange(x0,x1+1)+.5,np.arange(y0,y1+1)+.5)
        den=(q[1,1]-q[2,1])*(q[0,0]-q[2,0])+(q[2,0]-q[1,0])*(q[0,1]-q[2,1])
        if abs(den)<1e-12:return
        u=((q[1,1]-q[2,1])*(x-q[2,0])+(q[2,0]-q[1,0])*(y-q[2,1]))/den
        v=((q[2,1]-q[0,1])*(x-q[2,0])+(q[0,0]-q[2,0])*(y-q[2,1]))/den
        w=1-u-v;z=u*q[0,2]+v*q[1,2]+w*q[2,2]
        depth=self.water_depth if water else self.hull_depth
        colors=self.water_color if water else self.hull_color
        region=depth[y0:y1+1,x0:x1+1];mask=(u>=-1e-10)&(v>=-1e-10)&(w>=-1e-10)&(z>region)
        region[mask]=z[mask]
        if not water:self.hull_inside[y0:y1+1,x0:x1+1][mask]=interior
        pixelcolor=np.asarray(color)
        if pixelcolor.shape==(3,3):pixelcolor=u[...,None]*pixelcolor[0]+v[...,None]*pixelcolor[1]+w[...,None]*pixelcolor[2]
        colors[y0:y1+1,x0:x1+1][mask]=pixelcolor if pixelcolor.shape==(3,) else pixelcolor[mask]
    def line(self,segment,color,width=2,bias=.001):
        q=self.project(segment);n=max(2,int(np.linalg.norm(q[1,:2]-q[0,:2])*2))
        for p in np.linspace(q[0],q[1],n):
            x,y=int(round(p[0])),int(round(p[1]))
            for xx in range(max(0,x-width//2),min(self.width,x-width//2+width)):
                for yy in range(max(0,y-width//2),min(self.height,y-width//2+width)):
                    if p[2]+bias>=max(self.hull_depth[yy,xx],self.water_depth[yy,xx]):
                        self.hull_depth[yy,xx]=p[2]+bias;self.hull_color[yy,xx]=color
    def finish(self,transmission):
        background=np.array([.90,.95,.98]);rgb=np.broadcast_to(background,(self.height,self.width,3)).copy()
        h=np.isfinite(self.hull_depth);rgb[h]=self.hull_color[h]
        front=np.isfinite(self.water_depth)&(self.water_depth>self.hull_depth)
        rgb[front]=self.water_color[front]
        submerged=front&h
        optical_depth=self.water_depth[submerged]-self.hull_depth[submerged]
        transmitted=transmission*np.exp(-.6*optical_depth)
        rgb[submerged]=(1-transmitted[:,None])*self.water_color[submerged]+transmitted[:,None]*self.hull_color[submerged]
        return Image.fromarray(np.uint8(np.clip(rgb,0,1)*255))


def height(eta,x,y,lo,step):
    fx=np.clip((x-lo)/step,0,eta.shape[1]-1);fy=np.clip((y-lo)/step,0,eta.shape[0]-1)
    j=np.minimum(fx.astype(int),eta.shape[1]-2);i=np.minimum(fy.astype(int),eta.shape[0]-2)
    u=fx-j;v=fy-i;first=u+v<=1
    gx=np.where(first,eta[i,j+1]-eta[i,j],eta[i+1,j+1]-eta[i+1,j])/step
    gy=np.where(first,eta[i+1,j]-eta[i,j],eta[i+1,j+1]-eta[i,j+1])/step
    value=np.where(first,eta[i,j]+u*(eta[i,j+1]-eta[i,j])+v*(eta[i+1,j]-eta[i,j]),
        eta[i+1,j+1]+(1-u)*(eta[i+1,j]-eta[i+1,j+1])+(1-v)*(eta[i,j+1]-eta[i+1,j+1]))
    return value,gx,gy


def contact_segments(hull,faces,grid):
    water=np.concatenate([grid[:-1,:-1,None],grid[:-1,1:,None],grid[1:,:-1,None]],axis=2).reshape(-1,3,3)
    water=np.concatenate([water,np.stack([grid[1:,:-1],grid[:-1,1:],grid[1:,1:]],axis=2).reshape(-1,3,3)])
    result=[]
    for face in faces:
        tri=hull[face];normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);normal/=np.linalg.norm(normal)
        d=(water-tri[0])@normal
        candidates=water[(d.min(axis=1)<=0)&(d.max(axis=1)>=0)]
        if not len(candidates):continue
        d=(candidates-tri[0])@normal;b=np.roll(candidates,-1,axis=1);db=np.roll(d,-1,axis=1)
        valid=(d*db<=0)&(abs(d-db)>1e-14)
        ids=np.flatnonzero(valid.sum(axis=1)>=2)
        if not len(ids):continue
        candidates=candidates[ids];b=b[ids];d=d[ids];db=db[ids];valid=valid[ids]
        weights=np.divide(d,d-db,out=np.zeros_like(d),where=abs(d-db)>1e-14)
        points=candidates+weights[:,:,None]*(b-candidates)
        choices=np.argsort(~valid,axis=1)[:,:2];segments=points[np.arange(len(points))[:,None],choices]
        delta=segments[:,1]-segments[:,0];lower=np.zeros(len(segments));upper=np.ones(len(segments));keep=np.ones(len(segments),bool)
        for a,b in zip(tri,np.roll(tri,-1,axis=0)):
            outward=np.cross(b-a,normal);d0=(segments[:,0]-a)@outward;dd=delta@outward
            lower=np.maximum(lower,np.divide(-d0,dd,out=np.full_like(dd,-np.inf),where=dd< -1e-14))
            upper=np.minimum(upper,np.divide(-d0,dd,out=np.full_like(dd,np.inf),where=dd>1e-14))
            keep &= (abs(dd)>1e-14)|(d0<=1e-9)
        keep &= lower<=upper
        result.extend(np.stack([segments[:,0]+lower[:,None]*delta,segments[:,0]+upper[:,None]*delta],axis=1)[keep])
    return result


def render(frames,data,index,config,palette):
    started=time.monotonic();hull=frames['hull'][index];center=frames['center'][index]
    yaw=np.rad2deg(data['heading']['mean_box_yaw_global_rad'])
    raster=Raster(config['width_px'],config['height_px'],center,24,yaw-110,config['viewport_m'])
    projected=raster.project(hull)
    assert np.all((projected[:,0]>15)&(projected[:,0]<config['width_px']-15)&(projected[:,1]>85)&(projected[:,1]<config['height_px']-20)), 'Presentation hull canvas clipping'
    light=np.array([-.35,-.25,.90]);light/=np.linalg.norm(light);centroid=hull.mean(axis=0)
    closed=ConvexHull(hull);planes=closed.equations
    for face in data['faces']:
        triangle=hull[face]
        assert np.min(abs(triangle@planes[:,:3].T+planes[:,3]).max(axis=0))<1e-7
        normal=np.cross(triangle[1]-triangle[0],triangle[2]-triangle[0]);normal/=np.linalg.norm(normal)
        if (triangle.mean(axis=0)-centroid)@normal<0:normal=-normal
        interior=normal@raster.toward<0
        base=np.array([.73,.52,.32]) if interior else np.array([.93,.91,.84])
        raster.triangle(triangle,base*(.68+.32*abs(normal@light)),interior=interior)
    hull_seconds=time.monotonic()-started;water_started=time.monotonic()
    sine=raster.toward[2];horizontal=raster.toward[:2]/sine
    eta=frames['eta'][index];lo=float(frames['dx'][0,0]);step=float(frames['dx'][0,1]-lo);hi=-lo
    gx0=np.diff(eta,axis=1)/step;gy0=np.diff(eta,axis=0)/step
    directional=max(np.max(abs(gx0[:-1]*horizontal[0]+gy0[:,:-1]*horizontal[1])),np.max(abs(gx0[1:]*horizontal[0]+gy0[:,1:]*horizontal[1])))
    assert directional<1,'Selected view does not guarantee a unique heightfield ray hit'
    dark=np.array([.015,.115,.23]);pale=np.array([.77,.94,.99])
    def query(u,v):
        seaX=u*raster.right[0]+v*raster.up[0]-v*raster.up[2]*horizontal[0]
        seaY=u*raster.right[1]+v*raster.up[1]-v*raster.up[2]*horizontal[1]
        z=np.zeros_like(u)
        for iteration in range(12):
            value,gx,gy=height(eta,seaX+z*horizontal[0],seaY+z*horizontal[1],lo,step)
            correction=(z-value)/(1-gx*horizontal[0]-gy*horizontal[1]);z-=correction
            x=seaX+z*horizontal[0];y=seaY+z*horizontal[1]
            active=(x>=lo)&(x<=hi)&(y>=lo)&(y<=hi)
            if not active.any() or np.max(abs(correction[active]))<1e-9:break
        value,gx,gy=height(eta,x,y,lo,step)
        residual=float(np.max(abs(z[active]-value[active]))) if active.any() else 0.
        inside_volume=active.copy()
        for plane in planes:
            inside_volume &= plane[0]*(x+center[0])+plane[1]*(y+center[1])+plane[2]*z+plane[3]<=1e-9
        visible=active&~inside_volume
        depth=-v*raster.up[2]/sine+z/sine
        level=np.clip((z-palette[0])/(palette[1]-palette[0]),0,1)
        shade=.80+.20*np.clip((-gx*light[0]-gy*light[1]+light[2])/np.sqrt(1+gx*gx+gy*gy),0,1)
        color=(dark+level[...,None]*(pale-dark))*shade[...,None]
        return depth,color,visible,residual,iteration+1,int(inside_volume.sum())
    loww,lowh=config['water_pass_px']
    yy,xx=np.indices((lowh,loww),dtype=float)
    u=(xx+.5-loww/2)/(raster.scale*loww/raster.width)
    v=(lowh*.54-yy-.5)/(raster.scale*lowh/raster.height)
    depth,color,visible,lowres_residual,iterations,removed=query(u,v)
    def resize(array,nearest=False):
        return np.asarray(Image.fromarray(array.astype(np.float32)).resize((raster.width,raster.height),Image.Resampling.NEAREST if nearest else Image.Resampling.BILINEAR))
    raster.water_depth[:]=resize(depth)
    raster.water_depth[resize(visible,True)<.5]=-np.inf
    for channel in range(3):raster.water_color[:,:,channel]=resize(color[:,:,channel])
    # Keep all water/hull occlusion and contact at native pixel resolution.
    h=np.isfinite(raster.hull_depth);refine=binary_dilation(h,iterations=2)
    ry,rx=np.nonzero(refine)
    u=(rx+.5-raster.width/2)/raster.scale;v=(raster.height*.54-ry-.5)/raster.scale
    depth,color,visible,fullres_residual,fulliterations,removed_full=query(u,v)
    raster.water_depth[ry,rx]=np.where(visible,depth,-np.inf);raster.water_color[ry,rx]=color
    hyindex,hxindex=np.nonzero(h);u=(hxindex+.5-raster.width/2)/raster.scale;v=(raster.height*.54-hyindex-.5)/raster.scale
    depth=raster.hull_depth[hyindex,hxindex]
    hx=u*raster.right[0]+v*raster.up[0]+depth*raster.toward[0]
    hy=u*raster.right[1]+v*raster.up[1]+depth*raster.toward[1]
    hz=v*raster.up[2]+depth*sine
    hull_level=height(eta,hx,hy,lo,step)[0]
    wet=~raster.hull_inside[hyindex,hxindex]&(hz<hull_level)
    wy,wx=hyindex[wet],hxindex[wet]
    raster.hull_color[wy,wx]=.5*raster.hull_color[wy,wx]+.5*np.array([.12,.42,.57])
    water_seconds=time.monotonic()-water_started;contact_started=time.monotonic()
    xxgrid,yygrid=np.meshgrid(frames['dx'][0]+center[0],frames['dy'][:,0]+center[1])
    grid=np.stack([xxgrid,yygrid,eta],axis=-1);contact=contact_segments(hull,data['faces'],grid)
    for segment in contact:raster.line(segment,np.array([.005,.025,.035]),width=2)
    counts=Counter(tuple(sorted(edge)) for face in data['faces'] for edge in zip(face,np.roll(face,-1)))
    for edge,count in counts.items():
        if count==1:raster.line(hull[list(edge)],np.array([.20,.16,.12]),width=2)
    image=raster.finish(config['underwater_transmission'])
    # Outside the finite accepted patch is a solid background hue, not new water.
    background=(dark+(0-palette[0])/(palette[1]-palette[0])*(pale-dark))*(.80+.20*light[2])
    pix=np.asarray(image).copy();empty=~np.isfinite(raster.water_depth)&~h;pix[empty]=np.uint8(background*255)
    image=Image.fromarray(pix);draw=ImageDraw.Draw(image)
    try:font=ImageFont.truetype('C:/Windows/Fonts/segoeui.ttf',25)
    except OSError:font=ImageFont.load_default()
    q=data['q'][frames['indices'][index]];angles=np.rad2deg(q[3:])
    draw.text((35,25),f't {frames["time"][index]:05.2f} s   Roll {angles[0]:+.2f} deg   Pitch {angles[1]:+.2f} deg   Heave {q[2]:+.2f} m',font=font,fill='#102d3a')
    contact_residual=max((float(np.max(abs(segment[:,2]-height(eta,segment[:,0]-center[0],segment[:,1]-center[1],lo,step)[0]))) for segment in contact),default=0.)
    return image,dict(frame_index=int(index),time_s=float(frames['time'][index]),render_seconds=time.monotonic()-started,contact_segments=len(contact),contact_max_cached_eta_residual_m=contact_residual,ray_height_residual_m=max(lowres_residual,fullres_residual),ray_iterations=max(iterations,fulliterations),water_pass_px=config['water_pass_px'],full_resolution_water_hull_pixels=len(ry),hull_raster_seconds=hull_seconds,water_seconds=water_seconds,contact_and_composite_seconds=time.monotonic()-contact_started,directional_slope_ratio=float(directional),water_pixels_inside_volume_removed=removed_full,camera_elevation_deg=24)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--mode',choices=['preview','encode'],default='preview');parser.add_argument('--preview-review',type=Path);args=parser.parse_args()
    config=read(args.config);assert config['status']=='ready_for_coordinator_presentation_preview'
    assert all(sha(resolved(p))==v for p,v in config['input_sha256'].items())
    accepted=read(resolved(config['accepted_controls']));data=load_inputs(accepted)
    identity=dict(source_sha256=sha(resolved('reference/wecsim/render_calm_buoy_reference.py')),config_sha256=sha(resolved(config['accepted_controls'])),input_sha256=accepted['input_sha256'])
    frames,_=load_reviewed_cache(resolved(config['accepted_preview_review']),identity)
    identity=dict(source_sha256=sha(__file__),config_sha256=sha(args.config),input_sha256=config['input_sha256'])
    if args.mode=='encode':
        assert args.preview_review
        approved=read(args.preview_review);assert approved['decision']=='preview_pass' and approved['reviewed_by']
        report_path=resolved(approved['preflight_json']);assert sha(report_path)==approved['preflight_sha256']
        report=read(report_path);assert report['complete'] and report['mode']=='preview' and report['identity']==identity
    args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    budget=config['preview_maximum_seconds'] if args.mode=='preview' else config['encode_maximum_seconds']
    peak=int(np.argmax(abs(data['q'][:,3])));nearest=int(np.argmin(abs(frames['indices']-peak)))
    indices=sorted(set([0,375,750,nearest]));results=[];palette=[float(frames['eta'].min()),float(frames['eta'].max())]
    angle=np.deg2rad(np.rad2deg(data['heading']['mean_box_yaw_global_rad'])-110)
    horizontal=np.array([np.cos(angle),np.sin(angle)])/np.tan(np.deg2rad(24))
    step=float(frames['dx'][0,1]-frames['dx'][0,0]);gx=np.diff(frames['eta'],axis=2)/step;gy=np.diff(frames['eta'],axis=1)/step
    allframe_ratio=float(max(np.max(abs(gx[:,:-1]*horizontal[0]+gy[:,:,:-1]*horizontal[1])),np.max(abs(gx[:,1:]*horizontal[0]+gy[:,:,1:]*horizontal[1]))))
    assert allframe_ratio<1,'Selected24view has no guaranteed unique ray hit across all751 cachedframes'
    del gx,gy
    proc=None
    try:
        if args.mode=='encode':
            executable=shutil.which('ffmpeg');assert executable
            proc=subprocess.Popen([executable,'-y','-loglevel','error','-f','rawvideo','-pixel_format','rgb24','-video_size',f"{config['width_px']}x{config['height_px']}",'-framerate','12.5','-i','pipe:0','-an','-c:v','libx264','-b:v','3500k','-pix_fmt','yuv420p','-threads','1',str(args.output/'reference-presentation.mp4')],stdin=subprocess.PIPE,stderr=(args.output/'ffmpeg.log').open('wb'))
            print(dict(encoder_pid=proc.pid,parent_pid=os.getpid()),flush=True)
        for index in (indices if args.mode=='preview' else range(751)):
            check_deadline(started,budget);image,proof=render(frames,data,index,config,palette)
            if proc:proc.stdin.write(image.tobytes())
            if index in indices:
                name=f'elev24-frame{index:04d}.png';image.save(args.output/name);proof['file']=name
            results.append(proof)
            if proc and index%100==0:print(dict(frames_completed=index+1,elapsed_s=time.monotonic()-started),flush=True)
        if proc:
            proc.stdin.close();assert proc.wait(timeout=max(1,budget-(time.monotonic()-started)))==0
            print(dict(encoder_pid=proc.pid,exact_process_terminal=True,exit_code=proc.returncode),flush=True)
    finally:
        if proc and proc.poll() is None:proc.terminate();proc.wait(timeout=10)
    assert all(sha(resolved(p))==v for p,v in config['input_sha256'].items())
    write(args.output/'presentation-preview.json',dict(complete=True,mode=args.mode,identity=identity,all751_directional_slope_ratio=allframe_ratio,elapsed_s=time.monotonic()-started,estimated_full751_seconds=751*np.mean([p['render_seconds'] for p in results]),previews=results,palette_fixed_all751_z_range_m=palette,full_native_peak_roll_index=peak,cache_peak_preview_index=nearest,field='immutable accepted diagonal piecewise-linear cache watertriangles only; no wave regeneration',limitations='visual contact with prescribed incident surface, not flooding/wetted-body physics; native trajectory remains provisional'))
    write(args.output/'output-sha256.json',{p.name:sha(p) for p in args.output.iterdir() if p.is_file()})


if __name__=='__main__':
    with threadpool_limits(1):main()

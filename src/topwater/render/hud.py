"""Cached Pillow glyph/graticule textures; drawing is composited by WebGPU."""
from importlib.metadata import distribution
from importlib.resources import files
import numpy as np
from PIL import Image,ImageDraw,ImageFont

PALETTE=np.load(files(__package__).joinpath('coolwarm.npy'),allow_pickle=False)


def header_text(stamp,roll,pitch,heave):
    """Fixed-width columns for the monospace motion header."""
    return (f'TIME {stamp:7.2f} s   ROLL {roll:+7.2f}°   '
            f'PITCH {pitch:+7.2f}°   HEAVE {heave:+7.2f} m')


def water_colors(height,magnitude):
    indices=np.clip(((np.asarray(height)/magnitude+1)*128).astype(int),0,255)
    return PALETTE[indices]


class HUD:
    def __init__(self,size,lower,upper,magnitude,ticks):
        self.size=tuple(size)
        self.lower,self.upper=np.asarray(lower),np.asarray(upper)
        self.magnitude=magnitude
        self.ticks=ticks
        # Reuse the installed, licensed DejaVu font data without importing any
        # Matplotlib plotting or drawing implementation.
        self.font_path=str(distribution('matplotlib').locate_file('matplotlib/mpl-data/fonts/ttf/DejaVuSansMono.ttf'))
        self.header_font=ImageFont.truetype(self.font_path,max(1,round(14.6*self.size[0]/800)))
        self.header_height=max(1,round(25*self.size[0]/800))
        self.header_top=round(.075*self.size[1])
        span=self.upper-self.lower
        self.metadata=dict(common_metres_per_pixel=float(span[1]/self.size[1]),
            gauge_metres_per_pixel=float(span[1]/self.size[1]),world_y_zero_fraction=float(-lower[1]/span[1]),
            gauge_zero_fraction=float(-lower[1]/span[1]),edge_ticks=ticks)

    def graticule(self):
        factor=2
        width,height=np.array(self.size)*factor
        image=Image.new('RGBA',(int(width),int(height)))
        draw=ImageDraw.Draw(image)
        font=ImageFont.truetype(self.font_path,max(1,round(13.2*width/800)))
        def point(x,y):return (x*width,(1-y)*height)
        def line(xy,width_px=1.5):draw.line([point(*p) for p in xy],fill='black',width=max(1,round(width_px*self.size[0]/800*factor)))
        def text(x,y,value,anchor='mm'):draw.text(point(x,y),value,font=font,fill='black',anchor=anchor)
        for axis,name in enumerate(('x','y')):
            lo,hi=self.lower[axis],self.upper[axis]
            tick=self.ticks[name];major=np.asarray(tick['major_world_m']);step=tick['major_step_m']
            for value in tick['minor_world_m']:
                fraction=(value-lo)/(hi-lo)
                if not .08<=fraction<=.92:continue
                is_major=np.any(abs(major-value)<=step*1e-8)
                length=.020 if is_major else .010
                thickness=1.5 if is_major else 1.1
                if axis==0:
                    line([(fraction,.04),(fraction,.04+length)],thickness)
                    line([(fraction,.96),(fraction,.96-length)],thickness)
                else:
                    line([(.04,fraction),(.04+length,fraction)],thickness)
                    line([(.96,fraction),(.96-length,fraction)],thickness)
                if is_major and .12<=fraction<=.88:
                    label=f'{(0. if abs(value)<step*1e-10 else value):g}'
                    if axis==0:text(fraction,.075,label,anchor='mb')
                    else:text(.075,fraction,label,anchor='lm')
        text(.5,.009,'X (m)',anchor='mb')
        # Vertical world-Y label, independent of the data-coordinate viewport.
        label=Image.new('RGBA',(int(width*.1),int(height*.04)))
        ImageDraw.Draw(label).text((label.width/2,label.height/2),'Y (m)',font=font,fill='black',anchor='mm')
        label=label.rotate(90,expand=True)
        image.alpha_composite(label,(round(.009*width),round(height/2-label.height/2)))
        span_y=self.upper[1]-self.lower[1]
        bottom=(-self.magnitude-self.lower[1])/span_y
        top=(self.magnitude-self.lower[1])/span_y
        if not 0<=bottom<top<=1:
            raise ValueError('This preview requires the rounded height gauge inside the fixed world-Y view; no rescaling is applied')
        x0,x1=round(.875*width),round(.893*width)
        y0,y1=round((1-top)*height),round((1-bottom)*height)
        colors=(PALETTE[np.clip(np.linspace(255,0,y1-y0).astype(int),0,255)]*255).astype(np.uint8)
        gauge=Image.fromarray(np.broadcast_to(colors[:,None],(len(colors),max(1,x1-x0),4)).copy())
        image.alpha_composite(gauge,(x0,y0))
        for value in (-self.magnitude,0.,self.magnitude):
            fraction=(value-self.lower[1])/span_y
            line([(.870,fraction),(.875,fraction)])
            text(.865,fraction,('0' if value==0 else f'{value:+g}'),anchor='rm')
        minor=self.ticks['y']['minor_step_m']
        for value in np.arange(np.ceil(-self.magnitude/minor),np.floor(self.magnitude/minor)+1)*minor:
            if abs(value)<minor*1e-8 or abs(value)>=self.magnitude-minor*1e-8:continue
            fraction=(value-self.lower[1])/span_y
            line([(.872,fraction),(.875,fraction)],1.1)
        text(.884,top+.017,'Z (m)',anchor='mb')
        return np.asarray(image.resize(self.size,Image.Resampling.LANCZOS))

    def header(self,value):
        image=Image.new('RGBA',(self.size[0],self.header_height))
        draw=ImageDraw.Draw(image)
        center=(self.size[0]/2,self.header_height/2)
        box=draw.textbbox(center,value,font=self.header_font,anchor='mm')
        draw.rectangle((box[0]-3,box[1]-3,box[2]+3,box[3]+3),fill=(255,255,255,217))
        draw.text(center,value,font=self.header_font,fill='black',anchor='mm')
        return np.asarray(image)

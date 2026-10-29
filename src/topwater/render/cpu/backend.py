"""Native Skia CPU rasterization of the shared scene, without GPU emulation."""
import numpy as np
import skia
from ..backend import DrawingGeometry


class DrawingBackend(DrawingGeometry):
    def __init__(self,size,lower,upper):
        self.width,self.height=map(int,size)
        self.lower,self.upper=np.asarray(lower),np.asarray(upper)
        self.info=dict(device='Skia native CPU rasterizer',adapter_type='CPU',backend_type='Skia')
        self.rgba=np.empty((self.height,self.width,4),np.uint8)
        self.bgra=np.empty_like(self.rgba)
        info=skia.ImageInfo.Make(self.width,self.height,skia.ColorType.kBGRA_8888_ColorType,
            skia.AlphaType.kPremul_AlphaType)
        self.surface=skia.Surface.MakeRasterDirect(info,self.bgra)
        self.canvas=self.surface.getCanvas()
        self.images={}
        self._xy=self._faces=None

    def texture(self,key,rgba,*,rectangle=None):
        # Skia's native premultiplied BGRA blitter avoids converting every
        # mostly transparent HUD pixel through its generic shader per frame.
        rgba=np.asarray(rgba)
        bgra=np.empty_like(rgba)
        bgra[:,:,:3]=((rgba[:,:,[2,1,0]].astype(np.uint16)*rgba[:,:,3:4]+127)//255).astype(np.uint8)
        bgra[:,:,3]=rgba[:,:,3]
        image=skia.Image.fromarray(bgra,
            colorType=skia.ColorType.kBGRA_8888_ColorType,alphaType=skia.AlphaType.kPremul_AlphaType)
        self.images[key]=(image,(0,0) if rectangle is None else rectangle[:2])

    def _pixels(self,xy):
        return (np.asarray(xy)-self.lower)/(self.upper-self.lower)*[self.width,-self.height]+[0,self.height]

    def draw(self,water_vertices,water_faces,mask_loops,line_layers):
        self.canvas.clear(skia.ColorWHITE)
        xy=water_vertices[:,:2]
        if self._xy is None or not np.array_equal(xy,self._xy) or not np.array_equal(water_faces,self._faces):
            self._batches=[]
            max_nodes=np.iinfo(np.uint16).max+1
            # Skia indices are uint16. Larger meshes retain every triangle by
            # remapping bounded face groups, rather than truncating indices.
            groups=[water_faces] if len(xy)<=max_nodes else [
                water_faces[i:i+max_nodes//3] for i in range(0,len(water_faces),max_nodes//3)]
            for faces in groups:
                nodes,indices=np.unique(faces,return_inverse=True)
                pixel=(xy[nodes]+1)*[self.width/2,-self.height/2]+[0,self.height]
                points=[skia.Point(float(x),float(y)) for x,y in pixel]
                self._batches.append((nodes,points,indices.ravel().tolist()))
            self._xy=xy.copy()
            self._faces=np.array(water_faces,copy=True)
        color=np.clip(np.rint(water_vertices[:,2:]*255),0,255).astype(np.uint32)
        packed=(color[:,3]<<24)|(color[:,0]<<16)|(color[:,1]<<8)|color[:,2]
        for nodes,points,indices in self._batches:
            vertices=skia.Vertices.MakeCopy(skia.Vertices.kTriangles_VertexMode,
                points,colors=packed[nodes].tolist(),indices=indices)
            self.canvas.drawVertices(vertices,skia.Paint(Color=skia.ColorWHITE,AntiAlias=True),skia.BlendMode.kModulate)
        mask=skia.Path()
        for loop in mask_loops:
            points=self._pixels(loop)
            if not len(points):continue
            mask.moveTo(*points[0])
            for point in points[1:]:mask.lineTo(*point)
            mask.close()
        self.canvas.drawPath(mask,skia.Paint(Color=skia.ColorSetARGB(255,243,241,232),AntiAlias=True))
        for name,segments,color,width,dashed in line_layers:
            path=skia.Path()
            for start,end in self._pixels(np.asarray(segments).reshape(-1,2,2)):
                path.moveTo(*start);path.lineTo(*end)
            paint=skia.Paint(Color=skia.ColorSetARGB(*(round(v*255) for v in (color[3],*color[:3]))),
                AntiAlias=True,Style=skia.Paint.kStroke_Style,StrokeWidth=width)
            if dashed:paint.setPathEffect(skia.DashPathEffect.Make([5.,2.],0.))
            self.canvas.drawPath(path,paint)
        for key in ('graticule','header'):
            if key in self.images:
                image,(x,y)=self.images[key]
                self.canvas.drawImage(image,x,y)
        self.rgba[:,:,:3]=self.bgra[:,:,[2,1,0]]
        self.rgba[:,:,3]=self.bgra[:,:,3]
        return self.rgba

    def close(self):
        self.images.clear()
        self.canvas=self.surface=None

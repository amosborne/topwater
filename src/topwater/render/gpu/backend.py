"""Offscreen WebGPU drawing on physical GPU hardware."""
import numpy as np
import wgpu
from ..backend import DrawingGeometry


SHADER = """
struct Colored { @builtin(position) position: vec4f, @location(0) color: vec4f };
@vertex fn vertex(@location(0) xy: vec2f, @location(1) color: vec4f) -> Colored {
    var out: Colored; out.position=vec4f(xy,0.,1.); out.color=color; return out;
}
@fragment fn fragment(input: Colored) -> @location(0) vec4f { return input.color; }
struct Textured { @builtin(position) position: vec4f, @location(0) uv: vec2f };
@group(0) @binding(0) var image: texture_2d<f32>;
@group(0) @binding(1) var image_sampler: sampler;
@vertex fn texture_vertex(@location(0) xy: vec2f, @location(1) uv: vec2f) -> Textured {
    var out: Textured; out.position=vec4f(xy,0.,1.); out.uv=uv; return out;
}
@fragment fn texture_fragment(input: Textured) -> @location(0) vec4f {
    return textureSample(image,image_sampler,input.uv);
}
"""


def select_adapter():
    """Select hardware only; CPU GPU-emulation adapters are never used."""
    candidates=[a for a in wgpu.gpu.enumerate_adapters_sync()
                if a.info['adapter_type'] in ('DiscreteGPU','IntegratedGPU')]
    if not candidates:
        raise RuntimeError('No hardware WebGPU adapter is available')
    return max(candidates,key=lambda a:a.info['adapter_type']=='DiscreteGPU')


class DrawingBackend(DrawingGeometry):
    """RGBA triangle drawing, stencil winding fills and textured HUD layers.

    Four-sample coverage antialiasing applies to water, hull and line geometry.
    """
    def __init__(self, size, lower, upper):
        self.width,self.height=map(int,size)
        self.lower,self.upper=np.asarray(lower),np.asarray(upper)
        self.adapter=select_adapter()
        self.info=dict(self.adapter.info)
        self.device=self.adapter.request_device_sync()
        self.queue=self.device.queue
        self.buffers={}
        self._faces=None
        self.textures={}
        self.texture_groups={}
        self.texture_vertices={}
        self.multisample=self.device.create_texture(size=(self.width,self.height,1),
            format='rgba8unorm',sample_count=4,usage=wgpu.TextureUsage.RENDER_ATTACHMENT)
        self.output=self.device.create_texture(size=(self.width,self.height,1),
            format='rgba8unorm',usage=wgpu.TextureUsage.RENDER_ATTACHMENT|wgpu.TextureUsage.COPY_SRC)
        self.stencil=self.device.create_texture(size=(self.width,self.height,1),
            format='depth24plus-stencil8',sample_count=4,usage=wgpu.TextureUsage.RENDER_ATTACHMENT)
        shader=self.device.create_shader_module(code=SHADER)
        vertex=dict(module=shader,entry_point='vertex',buffers=[dict(array_stride=24,
            attributes=[dict(format='float32x2',offset=0,shader_location=0),
                        dict(format='float32x4',offset=8,shader_location=1)])])
        stencil=dict(format='depth24plus-stencil8',depth_write_enabled=False,depth_compare='always')
        def pipeline(**changes):
            fields=dict(layout='auto',vertex=vertex,primitive=dict(topology='triangle-list',cull_mode='none'),
                depth_stencil=stencil,multisample=dict(count=4),
                fragment=dict(module=shader,entry_point='fragment',targets=[dict(format='rgba8unorm')]))
            fields.update(changes)
            return self.device.create_render_pipeline(**fields)
        self.color_pipeline=pipeline()
        self.stencil_pipeline=pipeline(depth_stencil=dict(stencil,
            stencil_front=dict(compare='always',pass_op='increment-wrap'),
            stencil_back=dict(compare='always',pass_op='decrement-wrap')),
            fragment=dict(module=shader,entry_point='fragment',targets=[dict(format='rgba8unorm',write_mask=0)]))
        self.fill_pipeline=pipeline(depth_stencil=dict(stencil,stencil_write_mask=0,
            stencil_front=dict(compare='not-equal'),stencil_back=dict(compare='not-equal')))
        blend=dict(color=dict(src_factor='src-alpha',dst_factor='one-minus-src-alpha',operation='add'),
                   alpha=dict(src_factor='one',dst_factor='one-minus-src-alpha',operation='add'))
        self.texture_pipeline=pipeline(vertex=dict(module=shader,entry_point='texture_vertex',
            buffers=[dict(array_stride=16,attributes=[dict(format='float32x2',offset=0,shader_location=0),
                      dict(format='float32x2',offset=8,shader_location=1)])]),
            fragment=dict(module=shader,entry_point='texture_fragment',targets=[dict(format='rgba8unorm',blend=blend)]))
        self.sampler=self.device.create_sampler(mag_filter='linear',min_filter='linear')
        self.rgba=np.empty((self.height,self.width,4),np.uint8)

    def upload(self,key,data,*,usage=wgpu.BufferUsage.VERTEX):
        data=np.ascontiguousarray(data)
        buffer=self.buffers.get(key)
        if buffer is None or buffer.size<data.nbytes:
            if buffer is not None:buffer.destroy()
            buffer=self.device.create_buffer(size=max(4,data.nbytes),usage=usage|wgpu.BufferUsage.COPY_DST)
            self.buffers[key]=buffer
        if data.nbytes:self.queue.write_buffer(buffer,0,data)
        return buffer

    def texture(self,key,rgba,*,rectangle=None):
        rgba=np.ascontiguousarray(rgba,np.uint8)
        h,w=rgba.shape[:2]
        if key not in self.textures:
            texture=self.device.create_texture(size=(w,h,1),format='rgba8unorm',
                usage=wgpu.TextureUsage.TEXTURE_BINDING|wgpu.TextureUsage.COPY_DST)
            self.textures[key]=texture
            self.texture_groups[key]=self.device.create_bind_group(
                layout=self.texture_pipeline.get_bind_group_layout(0),entries=[
                    dict(binding=0,resource=texture.create_view()),dict(binding=1,resource=self.sampler)])
            # Pixel coordinates use a top-left origin; texture UV follows it.
            x0,y0,x1,y1=(0,0,self.width,self.height) if rectangle is None else rectangle
            xy=np.array([[x0,y1],[x1,y1],[x1,y0],[x0,y1],[x1,y0],[x0,y0]],np.float32)
            xy=xy/[self.width,self.height]*[2,-2]+[-1,1]
            uv=np.array([[0,1],[1,1],[1,0],[0,1],[1,0],[0,0]])
            self.texture_vertices[key]=self.upload('texture-'+key,np.c_[xy,uv].astype(np.float32))
        self.queue.write_texture(dict(texture=self.textures[key]),rgba,
            dict(bytes_per_row=w*4,rows_per_image=h),(w,h,1))

    def line_vertices(self,segments,color,width,*,dashed=False):
        segments=np.asarray(segments).reshape(-1,2,2)
        if not len(segments):return np.empty((0,6),np.float32)
        pixel=(segments-self.lower)/(self.upper-self.lower)*[self.width,self.height]
        if dashed:
            pieces=[]
            for a,b in pixel:
                length=np.linalg.norm(b-a)
                if length:
                    for offset in np.arange(0.,length,7.):
                        pieces.append([a+(b-a)*offset/length,a+(b-a)*min(offset+5.,length)/length])
            pixel=np.asarray(pieces).reshape(-1,2,2)
        delta=pixel[:,1]-pixel[:,0]
        length=np.linalg.norm(delta,axis=1)
        pixel,delta,length=pixel[length>0],delta[length>0],length[length>0]
        perpendicular=np.c_[-delta[:,1],delta[:,0]]/length[:,None]*width/2
        corners=np.stack((pixel[:,0]+perpendicular,pixel[:,1]+perpendicular,
            pixel[:,1]-perpendicular,pixel[:,0]-perpendicular),axis=1)
        xy=corners[:,[0,1,2,0,2,3]].reshape(-1,2)/[self.width,self.height]*2-1
        return self.vertices(xy,color)

    def draw(self,water_vertices,water_faces,mask_loops,line_layers):
        """Render one scene into top-left-origin RGBA pixels."""
        encoder=self.device.create_command_encoder()
        render=encoder.begin_render_pass(color_attachments=[dict(view=self.multisample.create_view(),
            resolve_target=self.output.create_view(),load_op='clear',store_op='discard',clear_value=(1.,1.,1.,1.))],
            depth_stencil_attachment=dict(view=self.stencil.create_view(),depth_load_op='clear',
                depth_store_op='discard',depth_clear_value=1.,stencil_load_op='clear',
                stencil_store_op='discard',stencil_clear_value=0))
        def triangles(key,values,pipeline):
            if len(values):
                render.set_pipeline(pipeline)
                render.set_vertex_buffer(0,self.upload(key,values))
                render.draw(len(values))
        if self._faces is None or not np.array_equal(self._faces,water_faces):
            self.water_indices=self.upload('water-indices',np.asarray(water_faces,dtype=np.uint32),
                usage=wgpu.BufferUsage.INDEX)
            self._faces=np.array(water_faces,copy=True)
        render.set_pipeline(self.color_pipeline)
        render.set_vertex_buffer(0,self.upload('water',water_vertices))
        render.set_index_buffer(self.water_indices,'uint32')
        render.draw_indexed(water_faces.size)
        fans=[]
        for loop in mask_loops:
            if len(loop)>=3:
                fans.append(np.stack((np.repeat(loop[:1],len(loop)-2,axis=0),loop[1:-1],loop[2:]),axis=1).reshape(-1,2))
        if fans:
            triangles('mask',self.vertices(self.to_clip(np.concatenate(fans)),[0,0,0,1]),self.stencil_pipeline)
            triangles('fill',self.vertices(self.fullscreen,[243/255,241/255,232/255,1]),self.fill_pipeline)
        for name,segments,color,width,dashed in line_layers:
            triangles(name,self.line_vertices(segments,color,width,dashed=dashed),self.color_pipeline)
        for key in ('graticule','header'):
            if key in self.textures:
                render.set_pipeline(self.texture_pipeline)
                render.set_bind_group(0,self.texture_groups[key])
                render.set_vertex_buffer(0,self.texture_vertices[key]);render.draw(6)
        render.end();self.queue.submit([encoder.finish()])
        raw=self.queue.read_texture(dict(texture=self.output),
            dict(bytes_per_row=self.width*4,rows_per_image=self.height),(self.width,self.height,1))
        self.rgba[:]=np.frombuffer(raw,np.uint8).reshape(self.height,self.width,4)
        return self.rgba

    def close(self):
        for resource in [*self.buffers.values(),*self.textures.values(),self.multisample,self.output,self.stencil]:
            resource.destroy()
        self.device.destroy()

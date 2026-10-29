# Render performance

Video rendering now shares scene preparation, physical contact, camera, palette,
HUD, cadence and encoding under `topwater/render/`. Mirrored `cpu/backend.py`
and `gpu/backend.py` implementations draw that scene with native Skia CPU
rasterization or hardware WebGPU. Automatic selection uses physical GPU hardware
when available, otherwise native CPU. No GPU-emulation adapter is used.

The current simple-boat-native archive contains 1,301 frames at 10 fps and
800 × 800 pixels. Complete measured renders include setup, geometry, drawing,
encoding, full decode verification and atomic publication. Archive loading is
recorded separately at 1.46 s.

| Same current native example | Complete render and verification |
|---|---:|
| Native CPU before indexed mesh/native image layout | 79.06 s |
| Final native CPU (Skia 144.1) | 53.02 s |
| Hardware GPU before indexed mesh | 29.20 s |
| Final hardware GPU (wgpu 0.32.0, RX 9070 XT/Vulkan) | 29.03 s |

The final CPU path is 33% faster than its first implementation. It caches
premultiplied BGRA HUD images and draws into a native BGRA raster surface,
converting the final frame to RGBA once for the encoder. This uses Skia's
optimized native image blitter instead of converting the mostly transparent
full-frame overlay through a generic shader on every frame.

Both backends receive the indexed water mesh directly: 7,396 vertices and
14,450 triangles in this case. They avoid expanding it to 43,350 vertex/color
records per frame. CPU caches geometry and local index mappings; GPU retains
the index buffer. Skia's uint16 index limit is handled by remapping bounded
triangle groups for larger meshes, preserving every triangle. GPU uses uint32
indices.

The GPU shader performs Gouraud water shading, winding/stencil hull fills,
antialiased wireframe/contact lines and HUD texture compositing. CPU uses native
Skia operations for the same scene. GPU uses four-sample coverage; Skia uses its
native coverage antialiasing. Pixel identity between those rasterizers is not
claimed. Startup, maximum-roll and final decoded frames were compared and
visually inspected. Mean channel differences were 0.22, 0.43 and 0.22 on the
0–255 scale; larger differences were confined to antialiased edges and encoding.
Pose, waterline, camera and meter-scaled HUD agree.

The focused 99-test verification covered CPU/GPU scene interiors, concave masks
with nested holes, non-aligned 642-pixel readback, indices above 65,535, varying
surface topology, exposed holes, saved resolution/cadence, failed-publication
preservation, and explicit non-rigid scenes with mixed hull winding. The unique
OS-temporary pytest base was removed after the run.

Physical hull/water intersection remains CPU work and H.264 encoding remains
CPU work. Those now account for much of the hardware render's remaining time.
Complete measurements, input/source hashes and sampled-frame differences are in
[the final benchmark report](../tests/data/native_render_backends.json).

## Earlier conservative contact optimization

The shared waterline kernel first discards water triangles outside the whole
hull's bounding box, retaining the existing tolerance, face order and exact
triangle-intersection algorithm. On the older 1,626-frame, 12.5 fps native
archive, that change improved the previous renderer from 64.16 s to a median
54.55 s. Baseline and optimized MP4s were byte-identical. These are a different
archive and rendering implementation, so they should not be read as a controlled
comparison with the current table. Evidence remains in
[the historical report](../tests/data/render_contact_culling.json).

An intermediate WebGPU software-emulation experiment took 228.50 s on the
current example. It was rejected and is not part of the final CPU path.

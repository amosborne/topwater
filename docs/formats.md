# Numeric file formats

NPZ files contain only named numeric arrays and load with `allow_pickle=False`.
Required members are exact: missing, extra, duplicate, object-valued or malformed
members are rejected. Float arrays must be finite. Connectivity/offsets are integers.
There are no embedded solver markers, provenance objects or restart metadata.

Units are metres, kilograms, seconds and radians. World axes are east/north/up;
body axes are bow/port/up, both right handed. Quaternions `(w,x,y,z)` actively map
body coordinates to world coordinates. Positions represent the world CG.

## boat.npz

| Array | Shape | Meaning |
|---|---|---|
| vertices_m | N x 3 | Complete physical hull in body coordinates |
| faces | F x 3 | Outward triangles; zero-based indices |
| mass_kg | scalar | Total dry mass including loading |
| cg_m | 3 | Body-frame centre of gravity |
| inertia_cg_kg_m2 | 3 x 3 | Full symmetric inertia tensor about CG |

Supply a convex, uncapped shell with exactly one opening, inferred from
single-incidence edges. Internal closures are geometric helpers, not wetted
pressure surfaces. The hull, CG and inertia must use the same body frame.

## sea.npz

Four equal-length vectors: `amplitudes_m`, `angular_frequencies_rad_s`,
`directions_rad`, `phases_rad`. Amplitudes are nonnegative and angular frequencies
positive. Empty arrays mean calm water. Direction is propagation CCW from world
east. Each carrier uses:

```text
amplitude * cos(k*(x*cos(direction) + y*sin(direction))
                - omega*(t - dwell_midpoint) + phase)
```

Rest/ramp/dwell are solve inputs; the sea archive stores original carriers and
does not duplicate the derived finite packet spectrum.

## result.npz

| Array | Shape | Meaning |
|---|---|---|
| time_s | T | Strictly increasing motion times starting at zero |
| position_m | T x 3 | Absolute world CG positions |
| orientation_wxyz | T x 4 | Unit body-to-world quaternions |
| hull_vertices_m | N x 3 | Exported full hull, body coordinates |
| hull_faces | F x 3 | Exported hull connectivity |
| cg_m | 3 | CG in exported hull coordinates |
| surface_time_s | S | Surface times, same endpoints as motion |
| surface_vertices_m | sum(Ns) x 3 | Packed world total water-air interface |
| surface_faces | sum(Fs) x 3 | Packed triangles, indices LOCAL to each frame |
| surface_vertex_offsets | S+1 | Frame boundaries in packed vertices |
| surface_face_offsets | S+1 | Frame boundaries in packed faces |

Offsets start at zero, increase strictly and end at the corresponding packed
array length. Pose a hull vertex with `position + R @ (vertex - cg)`.
The result's hull is authoritative for that result.

Surface geometry is mandatory and represents total water elevation. It preserves
the sampling grid without cutting/remeshing at the hull. Portions inside the
posed hull are invalid fluid and implicitly excluded by that hull; finite stored
values there are not physical predictions. Consumers compute intersections when
a displayed waterline is needed. No waterline, spectrum or validity array is stored.

`Boat`, `Sea` and `Result` in `topwater.workflow_data` read/write these archives.
The top-level physical input dataclasses remain available. The retired native
solver submodules are no longer public APIs. Format checks establish data
integrity, not physical accuracy.

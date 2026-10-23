"""Guarded GPU paired source potential with bounded CPU precision recovery."""
from dataclasses import dataclass
from decimal import DecimalException
from pathlib import Path
import numpy as np

from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence, QuadratureError
from topwater._rankine_extended import triangle_values
from .cuda_operator import cp

_source = (Path(__file__).with_name('cuda_pair.cuh').read_text()+'\n'
           +Path(__file__).with_suffix('.cu').read_text())
_kernel = cp.RawKernel(_source, 'paired_source_potential', options=('--std=c++11', '--fmad=false'))


@dataclass(frozen=True)
class SourcePotential:
    values: object
    diagnostics: dict


def _recover_pair(point, panels, column, image, controls):
    reflected = point.copy(); reflected[2] = 2*image-point[2]
    # Rejection guards have already run through public panel_influence.
    if np.array_equal(point, reflected):
        return 0.
    previous = None
    for refinement in range(min(controls.max_refinement+1, 5)):
        try:
            real, er = triangle_values(point, panels.triangles[column], panels.normals[column], 32+16*refinement)
            mirror, em = triangle_values(reflected, panels.triangles[column], panels.normals[column], 32+16*refinement)
        except (DecimalException, ValueError, ZeroDivisionError):
            previous = None
            continue
        if not all(np.isfinite(v) and v > 0 for v in [real[0], mirror[0]]):
            previous = None
            continue
        value = real[0]-mirror[0]
        error = er[0]+em[0]+4*np.finfo(float).eps*(abs(real[0])+abs(mirror[0])+abs(value))
        tolerance = controls.quadrature_atol*panels.diameters[column]+controls.quadrature_rtol*abs(value)
        if previous is not None and np.isfinite(value) and abs(value-previous)+error <= tolerance:
            return value
        previous = value
    raise QuadratureError('paired source cancellation exceeds requested binary64 accuracy')


def paired_potential(points, panels, image, controls=None):
    controls = LaplaceControls() if controls is None else controls
    points = np.asarray(points)
    image = np.asarray(image)
    if (not isinstance(panels, TrianglePanels) or not isinstance(controls, LaplaceControls)
            or np.iscomplexobj(points) or points.ndim != 2 or points.shape[1:] != (3,)
            or not np.isfinite(points).all() or image.shape != () or image.dtype.kind not in 'fiu'
            or not np.isfinite(image)):
        raise ValueError('finite bulk points, triangle panels, controls and image height required')
    points = np.ascontiguousarray(points, dtype=np.float64)
    image = float(image)
    count = len(points)*len(panels)
    if count*9 > controls.max_operator_bytes:
        raise MemoryError('paired source matrix and guard flags exceed the operator memory budget')
    values = cp.empty((len(points), len(panels)), dtype=cp.float64)
    flags = cp.empty(count, dtype=cp.uint8)
    if count:
        geometry = [cp.asarray(np.ascontiguousarray(a)) for a in
            [panels.triangles, panels.geometric_normals, panels.edge_outward, panels.tangents,
             panels.edge_lengths, panels.areas, panels.centroids, panels.normal_sign, panels.diameters]]
        _kernel(((count+127)//128,), (128,), (*geometry, cp.asarray(points), np.int32(len(points)),
            np.int32(len(panels)), np.float64(image), np.float64(controls.quadrature_atol),
            np.float64(controls.quadrature_rtol), values, flags))
    indices = cp.asnumpy(cp.flatnonzero(flags))
    status = cp.asnumpy(flags[cp.asarray(indices)])
    rows = np.unique(indices//len(panels))
    # Recompute flagged rows with unchanged public geometry/accuracy guards.
    # No geometry rejection is converted into a guessed value.
    for first in range(0, len(rows), controls.target_block_size):
        selected = rows[first:first+controls.target_block_size]
        reflected = points[selected].copy(); reflected[:, 2] = 2*image-reflected[:, 2]
        panel_influence(points[selected], panels, controls=controls)
        panel_influence(reflected, panels, controls=controls)
    for index, bits in zip(indices, status, strict=True):
        if bits & 6:
            row, column = divmod(int(index), len(panels))
            values[row, column] = _recover_pair(points[row], panels, column, image, controls)
    return SourcePotential(values, dict(gpu_flagged_pairs=len(indices), cpu_guarded_rows=len(rows),
        individual_accuracy_flags=int(np.count_nonzero(status & 2)),
        paired_accuracy_flags=int(np.count_nonzero(status & 4)), geometry_guard_flags=int(np.count_nonzero(status & 1)),
        matrix_bytes=int(values.nbytes), guard_bytes=int(flags.nbytes),
        scope='Same analytic sample primitive; explicit nonself points; CPU guards/recovery for flagged pairs'))

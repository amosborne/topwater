"""Blocked Laplace operator assembly; boundary solves are a separate layer."""

from dataclasses import asdict, dataclass
import json
from time import perf_counter
from zipfile import ZipFile

import numpy as np

from .numerics import DenseOperatorPlan, LaplaceControls
from .rankine import TrianglePanels, _evaluate_block


@dataclass(frozen=True, eq=False)
class LaplaceOperators:
    """Physical S and principal-value K (no half jump), plus measured metadata."""

    single_layer: np.ndarray
    double_layer: np.ndarray
    diagnostics: dict

    def save(self, path):
        """Save actual matrices and metadata without pickle or compression cost."""
        np.savez(path, single_layer=self.single_layer, double_layer=self.double_layer,
                 metadata=json.dumps(self.diagnostics, allow_nan=False))

    @classmethod
    def load(cls, path, *, controls=None):
        """Validate storage/version/trace ownership before accepting saved operators."""
        controls = LaplaceControls() if controls is None else controls
        # Inspect member headers before NumPy allocates an array. Metadata alone
        # cannot enforce a memory budget on a corrupt or mislabeled archive.
        with ZipFile(path) as archive:
            names = archive.namelist()
            expected = {"metadata.npy", "single_layer.npy", "double_layer.npy"}
            if len(names) != 3 or set(names) != expected:
                raise ValueError("invalid Laplace operator archive members")
            if archive.getinfo("metadata.npy").file_size > 65536:
                raise ValueError("operator metadata exceeds 64 KiB")
            shapes = []
            for name in ("single_layer.npy", "double_layer.npy"):
                with archive.open(name) as member:
                    version = np.lib.format.read_magic(member)
                    if version == (1, 0):
                        shape, fortran, dtype = np.lib.format.read_array_header_1_0(member)
                    elif version == (2, 0):
                        shape, fortran, dtype = np.lib.format.read_array_header_2_0(member)
                    else:
                        raise ValueError("unsupported operator array format")
                    if len(shape) != 2 or shape[0] != shape[1] or dtype != np.dtype("float64") or fortran:
                        raise ValueError("invalid saved operator shape, dtype or layout")
                    plan = DenseOperatorPlan(shape[0], controls)
                    if archive.getinfo(name).file_size != member.tell()+plan.operator_bytes//2:
                        raise ValueError("operator member size does not match header")
                    shapes.append(shape)
            if shapes[0] != shapes[1]:
                raise ValueError("operator shapes differ")
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata"]))
            if type(metadata.get("format_version")) is not int or metadata["format_version"] != 1:
                raise ValueError("unsupported Laplace operator format_version")
            if metadata.get("double_layer_trace") != "principal_value":
                raise ValueError("double-layer trace convention must be principal_value")
            plan = DenseOperatorPlan(metadata["panel_count"], controls)
            if shapes[0] != (plan.panel_count, plan.panel_count):
                raise ValueError("operator shape does not match metadata")
            arrays = []
            for name in ("single_layer", "double_layer"):
                array = data[name]
                if (array.shape != (plan.panel_count, plan.panel_count)
                        or array.dtype != np.dtype("float64") or not np.all(np.isfinite(array))):
                    raise ValueError(f"invalid saved {name} shape, dtype or values")
                array.setflags(write=False)
                arrays.append(array)
            if metadata.get("operator_bytes") != plan.operator_bytes:
                raise ValueError("operator_bytes does not match saved matrix sizes")
            return cls(*arrays, metadata)


def assemble_laplace(panels, *, controls=None):
    """Assemble independent analytic entries without modifying diagonal row sums.

    Geometry preparation is outside assembly wall time. Peak process RSS and
    actual thread/library metadata are collected by the fresh-process benchmark,
    not inferred from allocated operator sizes here.
    """
    if not isinstance(panels, TrianglePanels):
        raise ValueError("panels must be cached TrianglePanels")
    controls = LaplaceControls() if controls is None else controls
    plan = DenseOperatorPlan(len(panels), controls)
    start = perf_counter()
    single = np.empty((len(panels), len(panels)))
    double = np.empty_like(single)
    fallback_pairs = evaluations = 0
    for first in range(0, len(panels), plan.block_rows):
        last = min(first+plan.block_rows, len(panels))
        result = _evaluate_block(panels.centroids[first:last], panels, controls=controls,
                                 self_indices=np.arange(first, last),
                                 outputs=(single[first:last], double[first:last]))
        fallback_pairs += result.fallback_pairs
        evaluations += result.quadrature_evaluations
    elapsed = perf_counter()-start
    single.setflags(write=False)
    double.setflags(write=False)
    metadata = dict(
        format_version=1, double_layer_trace="principal_value", backend="numba-analytic",
        panel_count=len(panels), operator_bytes=plan.operator_bytes,
        target_block_size=plan.block_rows, assembly_seconds=elapsed,
        fallback_pairs=fallback_pairs, quadrature_evaluations=evaluations,
        controls=asdict(controls), estimated_error_max=None, peak_rss_bytes=None,
    )
    return LaplaceOperators(single, double, metadata)

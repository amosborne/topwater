from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from topwater.numerics import DenseOperatorPlan, LaplaceControls


def test_dense_storage_and_portable_preflight_record():
    # Storage arithmetic is independent of any application mesh.
    plan = DenseOperatorPlan(2048)
    assert plan.operator_bytes == 64 * 1024**2
    fixture = (Path(__file__).parent / "fixtures/laplace-preflight.json").read_text()
    assert json.loads(plan.to_json()) == json.loads(fixture)
    assert DenseOperatorPlan.from_json(fixture) == plan
    assert LaplaceControls.from_json(plan.controls.to_json()) == plan.controls


def test_budget_boundary_and_large_mesh_rejection_before_allocation():
    size = 2 * 8 * 2048**2
    controls = LaplaceControls(max_operator_bytes=size)
    DenseOperatorPlan(2048, controls)
    with pytest.raises(MemoryError, match="excludes solve/scratch"):
        DenseOperatorPlan(2048, replace(controls, max_operator_bytes=size - 1))
    with pytest.raises(MemoryError):
        DenseOperatorPlan(32768)


def test_small_problem_uses_only_required_block_rows():
    assert DenseOperatorPlan(np.int64(16)).block_rows == 16
    assert DenseOperatorPlan(512, LaplaceControls(target_block_size=64)).block_rows == 64


@pytest.mark.parametrize(("field", "value"), [
    ("quadrature_rtol", 0), ("quadrature_atol", float("nan")),
    ("quadrature_rtol", True), ("max_refinement", -1),
    ("target_block_size", 1.5), ("target_block_size", False),
    ("max_operator_bytes", 0),
])
def test_invalid_controls_cannot_reach_numerical_work(field, value):
    with pytest.raises(ValueError, match=field):
        LaplaceControls(**{field: value})


def test_corrupted_plan_cannot_misreport_storage():
    data = json.loads(DenseOperatorPlan(512).to_json())
    data["operator_bytes"] = 1
    with pytest.raises(ValueError, match="operator_bytes"):
        DenseOperatorPlan.from_json(json.dumps(data))
    data = json.loads(LaplaceControls().to_json())
    data["format_version"] = True
    with pytest.raises(ValueError, match="format_version"):
        LaplaceControls.from_json(json.dumps(data))

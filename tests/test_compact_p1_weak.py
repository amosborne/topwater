"""Resource failures must occur before the expensive allocation/kernel launch."""
import pytest
pytest.importorskip('cupy')
from topwater.rankine import QuadratureError
from tests.experiments.free_surface import compact_p1_weak as weak


def unexpected(*args,**kwargs):
    raise AssertionError('resource preflight allowed expensive numerical work')


def test_oversized_patch_rejected_before_device_allocation(monkeypatch):
    monkeypatch.setattr(weak.cp,'asarray',unexpected)
    with pytest.raises(MemoryError,match='storage budget'):
        weak.WeakPatch((10000,10000),.2,.2)


def test_sample_budget_rejected_before_raw_kernel_launch(monkeypatch):
    monkeypatch.setattr(weak,'_lags',unexpected)
    with pytest.raises(QuadratureError,match='sample/term budget'):
        weak.WeakPatch((2,3),.31,.23,weak.WeakControls(maximum_evaluations=1))

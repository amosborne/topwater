import numpy as np
import pytest
from .finest_cost_probe import sample_starts,projection


def test_samples_cover_both_halves_and_rare_facets_without_crossing_split():
    labels=np.zeros(2000,int);labels[444]=1;labels[1555]=2
    starts=sample_starts(labels)
    assert {0,450,900,1000,1450,1900}<=set(starts)
    for base in (0,1000):
        rows=np.concatenate([np.arange(s,s+100) for s in starts if base<=s<base+1000])
        assert set(labels[rows])==set(labels[base:base+1000])
        assert rows.min()>=base and rows.max()<base+1000
    with pytest.raises(ValueError):sample_starts(labels,maximum=1)


def test_probe_projection_uses_slowest_block_and_complete_field_pairs():
    value=projection([.1,.3,.2],3456,2388787200,12.)
    assert value['projected_field_seconds']==pytest.approx(2073.6)
    assert value['projected_upper_seconds']==pytest.approx((2073.6+12)*1.25+60)
    with pytest.raises(ValueError):projection([np.nan],3456,2388787200,12.)

import json
import numpy as np
import pytest
from .common import read_owned,digest
from .fixed_field import source_interval,continuation_plan


def test_source_interval_charges_both_uncertainties():
    coarse=np.ones((6,6));fine=coarse.copy();u=np.zeros((6,6));v=u.copy()
    u[0,5]=.0006;v[0,5]=.0006
    d=source_interval(coarse,fine,u,v,np.ones((6,6)))
    assert d['source_interval_failed'].sum()==1
    assert d['source_interval_ratio'][0,5]==pytest.approx(.0012/.001001)
    with pytest.raises(ValueError):source_interval(coarse,fine,-u,v,np.ones((6,6)))
    with pytest.raises(ValueError):source_interval(coarse[:3],fine,u,v,np.ones((6,6)))


def test_continuation_plan_counts_native_and_all_retained_buffers():
    old=dict(orders=[dict(seconds=2.,points=100),dict(seconds=1.,points=200)],native_seconds=5.,latest_memory=dict(peak_committed_bytes=1000))
    p=continuation_plan(old)
    assert p['new_field_pairs']==1161135360
    assert p['projected_upper_seconds']==pytest.approx((2./277200*1161135360+5.)*1.25+60.)
    assert p['retained_field_bytes']==542080*6*16
    assert p['target_bytes']==542080*64
    assert p['estimated_peak_bytes']==1000+542080*160+160*1024**2


def test_selected_owned_arrays_keep_full_identity_checks(tmp_path):
    path=tmp_path/'input.json';archive=path.with_suffix('.source');archive.mkdir()
    (archive/'contract.json').write_text('{}');(archive/'producer.py').write_text('producer')
    np.savez(path.with_suffix('.npz'),wanted=np.ones((6,6)),unused=np.arange(20))
    r=dict(kind='example',complete=True,passed=False,arrays_sha256=digest(path.with_suffix('.npz')),
        contract_sha256=digest(archive/'contract.json'),producer_sha256={'producer.py':digest(archive/'producer.py')})
    path.write_text(json.dumps(r))
    _,a=read_owned(path,r['contract_sha256'],'example',selected=['wanted'])
    assert set(a)=={'wanted'}
    with pytest.raises(ValueError):read_owned(path,r['contract_sha256'],'example',selected=['missing'])
    (archive/'producer.py').write_text('changed')
    with pytest.raises(ValueError):read_owned(path,r['contract_sha256'],'example',selected=['wanted'])

import json
import numpy as np
import pytest
from .common import read_owned,digest
from .fixed_field import source_interval


def test_source_interval_charges_both_uncertainties():
    coarse=np.ones((6,6));fine=coarse.copy();u=np.zeros((6,6));v=u.copy()
    u[0,5]=.0006;v[0,5]=.0006
    d=source_interval(coarse,fine,u,v,np.ones((6,6)))
    assert d['source_interval_failed'].sum()==1
    assert d['source_interval_ratio'][0,5]==pytest.approx(.0012/.001001)
    with pytest.raises(ValueError):source_interval(coarse,fine,-u,v,np.ones((6,6)))
    with pytest.raises(ValueError):source_interval(coarse[:3],fine,u,v,np.ones((6,6)))


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

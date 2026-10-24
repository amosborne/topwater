from io import BytesIO
import json
from hashlib import sha256
from zipfile import ZipFile
import numpy as np
import pytest
from .finest_consumer import consume_orders,verify_arrays,checkpoint_arrays
from .finest_cost_probe import owned_zip
from .frozen_source_quadrature import evaluate_trace


def test_consumer_retains_same_caller_owned_buffers_through_all_checkpoints():
    data={};a={'world_cg':np.zeros(3)};seen=[]
    for order,n in [(12,2),(16,3),(20,4)]:
        key=f'target_{order}_';data.update({key+'points':np.zeros((n,3)),key+'owners':np.zeros(n,int),key+'weights':np.ones(n),key+'normal':np.tile([0,0,1.],(n,1))})
    def fill(order,points,out):out[:]=order+1j;seen.append(out)
    def checkpoint(order):
        assert all(a[f'order_{o}_potential'] is value for o,value in zip([12,16,20],seen))
        assert a[f'target_{order}_points'] is data[f'target_{order}_points']
    values=consume_orders(data,[12,16,20],a,fill,lambda:None,checkpoint)
    assert values.shape==(3,6,6) and sum(x.nbytes for x in seen)==9*6*16


def test_native_evaluator_fills_exact_supplied_buffer():
    class Engine:
        def build_S_matrix(self,points,mesh,**kwargs):return np.ones((len(points),2))
    out=np.empty((5,6),complex)
    assert evaluate_trace(Engine(),None,np.zeros((5,3)),np.ones((2,6)),2,lambda:None,out=out) is out
    np.testing.assert_array_equal(out,2)
    with pytest.raises(ValueError):evaluate_trace(Engine(),None,np.zeros((5,3)),np.ones((2,6)),2,lambda:None,out=np.empty((5,6)))


def test_checkpoint_roundtrip_rejects_changed_bytes(tmp_path):
    arrays={'field':np.arange(12).reshape(2,6)+1j};output=tmp_path/'arrays.json';path=output.with_suffix('.npz');report={}
    checkpoint_arrays(output,report,arrays)
    expected=sha256(path.read_bytes()).hexdigest();assert report['arrays_sha256']==expected;verify_arrays(path,arrays,expected)
    changed=dict(arrays);changed['field']=arrays['field']+1
    with pytest.raises(AssertionError):verify_arrays(path,changed,expected)
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError):verify_arrays(path,arrays,expected)


def test_selected_historical_loading_still_checks_complete_payload_identity(tmp_path):
    raw=BytesIO();np.savez(raw,needed=np.ones(2),unneeded=np.ones(20));payload=raw.getvalue();contract=b'owned'
    report=dict(kind='test',complete=True,contract_sha256=sha256(contract).hexdigest(),arrays_sha256=sha256(payload).hexdigest(),producer_sha256={})
    bundle=tmp_path/'input.zip'
    with ZipFile(bundle,'x') as z:
        z.writestr('validation/test.json',json.dumps(report));z.writestr('validation/test.npz',payload);z.writestr('validation/test.source/contract.json',contract)
    _,arrays=owned_zip(bundle,'test','test',report['contract_sha256'],['needed'])
    assert set(arrays)=={'needed'}
    with pytest.raises(ValueError):owned_zip(bundle,'test','test',report['contract_sha256'],['missing'])

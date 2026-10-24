import pytest
from .common_output_preflight import work_model,resource_model


def test_finest_work_includes_all_orders_native_calibration_and_both_matrices():
    w=work_model(3456,[124416,221184,345600])
    assert w['field_pairs']==2388787200 and w['native_calibration_pairs']==11943936
    assert w['native_matrix_bytes']==2*100*3456*16
    assert w['new_field_bytes']==(691200+3456)*6*16
    assert w['block_calls']==6948


def test_resource_gate_rejects_complete_job_even_when_single_order_would_fit():
    c=dict(time_safety_factor=1.25,preparation_allowance_seconds=60,native_scratch_multiplier=4,
        runtime_overhead_bytes=0,future_maximum_elapsed_seconds=1800,maximum_peak_bytes=10**10,maximum_additional_disk_bytes=10**10)
    model=resource_model(work_model(3456,[124416,221184,345600]),[1e-6],0,0,0,0,c)
    assert model['calibrated_native_seconds'][0]>2400 and not model['ceiling_checks']['time']
    assert 3456*345600*1e-6*1.25+60<1800
    with pytest.raises(ValueError):resource_model(work_model(3456,[1,2,3]),[],0,0,0,0,c)

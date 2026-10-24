"""Structural validation gates only; no Capytaine/fluid/equilibrium solve."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from .common import entrywise,mesh_quality,study_preflight,read_owned,digest
from .reference import require_qualified
from . import __main__ as cli


ACCEPTANCE=dict(reference_mesh_absolute_normalized=1e-6,reference_mesh_relative=.001)


def test_each_physical_entry_has_its_own_allowance():
    reference=np.eye(6)*1e9;candidate=reference.copy();candidate[0,4]=1e-3
    assert not entrywise(candidate,reference,np.zeros((6,6)),1e-5,.01)['passed']


def test_complex_phase_and_reference_uncertainty_are_not_hidden():
    assert not entrywise(np.array([1j]),np.array([1+0j]),np.zeros(1),1e-5,.01)['passed']
    assert not entrywise(np.ones(1),np.ones(1),np.array([.02]),1e-5,.01)['passed']


def test_both_finest_increments_and_trend_are_required():
    assert not mesh_quality(np.array([0.,1.,1.]),ACCEPTANCE)[1]['passed']
    assert not mesh_quality(np.array([1.,1.00002,1.00005]),ACCEPTANCE)[1]['passed']
    assert mesh_quality(np.array([1.,1.00002,1.00003]),ACCEPTANCE)[1]['passed']
    with pytest.raises(ValueError):
        mesh_quality(np.ones(2),ACCEPTANCE)


def test_actual_wet_count_cannot_silently_drop_finest_mesh():
    plan=study_preflight(17,[1,2,3,4],8*1024**3)
    assert plan[-1]['panels']==17*4**4
    with pytest.raises(MemoryError):
        study_preflight(68,[1,2,3,4],8*1024**3)


def test_missing_finite_cases_cannot_qualify_reference():
    contract=dict(reference_numerics=dict(body_subdivision_levels=[1,2,3],capytaine_version='3.0.0'),
                  cases=dict(frequencies_hz=[.1],propagation_headings_rad=[0.]))
    report=dict(levels=[dict(level=n,complete=True) for n in [1,2,3]],capytaine_version='3.0.0')
    with pytest.raises((ValueError,KeyError)):
        require_qualified(report,dict(Ainf=np.ones((3,6,6))),contract,{})


def test_corrupted_or_wrong_kind_artifact_fails(tmp_path):
    path=tmp_path/'reference.json';archive=path.with_suffix('.source');archive.mkdir()
    (archive/'contract.json').write_text('{}');(archive/'owner.py').write_text('original')
    np.savez(path.with_suffix('.npz'),A=np.eye(6))
    report=dict(kind='reference',contract_sha256=digest(archive/'contract.json'),complete=True,passed=True,
        arrays_sha256=digest(path.with_suffix('.npz')),producer_sha256={'owner.py':digest(archive/'owner.py')})
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        read_owned(path,report['contract_sha256'],'compare-inertia',passed=True)
    (archive/'owner.py').write_text('changed')
    with pytest.raises(ValueError):
        read_owned(path,report['contract_sha256'],'reference',passed=True)


def test_unsupported_wave_gate_cannot_be_enabled_by_contract_flag():
    report={}
    cli.compare_waves(None,dict(time_domain_wave_controls={'status':'supported'}),report,{})
    assert report['complete'] and not report['passed'] and report['status']=='unsupported'


def test_qualify_cannot_accept_fabricated_passed_wave_artifact(monkeypatch):
    args=SimpleNamespace(frozen=Path('frozen'),reference=Path('reference'),inertia=Path('inertia'),waves=Path('waves'))
    monkeypatch.setattr(cli,'frozen_input',lambda *a:({'input_sha256':'same'},{}))
    monkeypatch.setattr(cli,'digest',lambda path:str(path))
    monkeypatch.setattr(cli,'hashes',lambda paths:{'candidate':'owned'})
    monkeypatch.setattr(cli,'qualified_reference',lambda *a:({},{}))
    monkeypatch.setattr(cli,'read_owned',lambda *a,**k:(dict(input_sha256='same',frozen_sha256='frozen',
        candidate_source_sha256={'candidate':'owned'},reference_sha256='reference'),{}))
    report={'contract_sha256':'contract'}
    cli.qualify(args,{},report,{})
    assert not report['passed']
    assert not report['gates']['small_amplitude_time_domain_wave_responses']['passed']

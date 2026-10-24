"""Structural validation gates only; no Capytaine/fluid/equilibrium solve."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from .common import (entrywise,mesh_quality,study_preflight,read_owned,digest,
                     load_contract,case_envelope,physical_identity,hashes,CANDIDATE,VALIDATORS,REFERENCE)
from .reference import require_qualified
from . import __main__ as cli


ACCEPTANCE=dict(reference_mesh_absolute_normalized=1e-6,reference_mesh_relative=.001)


def contract_data():
    return dict(format_version=1,
        implementation=dict(entry_point='topwater.inertia.instantaneous_added_mass'),
        physical_problem=dict(fixture='tests.reference_boat',body_function='reference_body',pose_function='equilibrium',
            water_depth='infinite',forward_speed_m_s=0.,initial_body_velocity='zero',
            time_convention='Re(amplitude*exp(-i*omega*t))',wave_direction='propagation, radians CCW from world+x',
            incident_phase_origin_world_xy_m=[0.,0.],incident_elevation_amplitude_m=.01,
            density_kg_m3=1025.,gravity_m_s2=9.81,mode_order=['Surge','Sway','Heave','Roll','Pitch','Yaw']),
        reference_numerics=dict(body_subdivision_levels=[1,2,3],cpu_threads=1,richardson_extrapolation=False,
            maximum_elapsed_seconds=1800,maximum_estimated_workspace_bytes=8*1024**3),
        candidate_numerics=dict(body_subdivision_levels=[1,2],cpu_threads=1),
        cases=dict(frequencies_hz=[.1,.5],propagation_headings_rad=[0.,1.]),
        acceptance=dict(reference_required_finest_increments=2,
            required_gates=['reference','all_six_added_mass','small_amplitude_time_domain_wave_responses'],
            candidate_relative=.01,candidate_absolute_normalized=1e-5,**ACCEPTANCE))


@pytest.mark.parametrize('key,value',[
    ('time_convention','Re(amplitude*exp(+i*omega*t))'),('wave_direction','from north'),
    ('incident_phase_origin_world_xy_m',[1.,0.]),('initial_body_velocity','moving'),
    ('water_depth',4.),('forward_speed_m_s',1.),('incident_elevation_amplitude_m',0.),
    ('incident_elevation_amplitude_m',float('nan'))])
def test_changed_physical_conventions_are_rejected(tmp_path,key,value):
    contract=contract_data();contract['physical_problem'][key]=value
    path=tmp_path/'contract.json';path.write_text(json.dumps(contract))
    with pytest.raises(ValueError):
        load_contract(path)


@pytest.mark.parametrize('section,key,value',[
    ('reference_numerics','maximum_elapsed_seconds',float('nan')),
    ('reference_numerics','maximum_estimated_workspace_bytes',float('inf')),
    ('cases','frequencies_hz',[.1,.1]),('cases','propagation_headings_rad',[0.,0.]),
    ('cases','propagation_headings_rad',[0.,2*np.pi]),
    ('cases','frequencies_hz',[.7])])
def test_nonfinite_limits_duplicate_cases_or_outofband_frequency_fail(tmp_path,section,key,value):
    contract=contract_data();contract[section][key]=value
    path=tmp_path/'contract.json';path.write_text(json.dumps(contract))
    with pytest.raises(ValueError):
        load_contract(path)


def test_actual_beam_wavelength_bound_allows_equality():
    contract=contract_data()
    with pytest.raises(ValueError):
        case_envelope(contract,10.)
    minimum=9.81/(2*np.pi*.5**2)
    case_envelope(contract,minimum/2)


def test_stable_physical_identity_ignores_execution_metadata_but_owns_data():
    report=dict(contract_sha256='contract',input_sha256='input',elapsed_seconds=1.,current_commit='first')
    arrays=dict(wetted_triangles=np.arange(9,dtype='<f8').reshape(1,3,3),dry_mass=np.eye(6))
    identity=physical_identity(report,arrays)
    changed=dict(report,elapsed_seconds=100.,current_commit='next',passed=True)
    assert physical_identity(changed,arrays)==identity
    assert physical_identity(dict(report,input_sha256='other'),arrays)!=identity
    moved={name:value.copy() for name,value in arrays.items()};moved['wetted_triangles'][0,0,0]=.01
    assert physical_identity(report,moved)!=identity
    assert physical_identity(report,dict(arrays,dry_mass=arrays['dry_mass'].astype('<f4')))!=identity


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
    monkeypatch.setattr(cli,'frozen_input',lambda *a:({'input_sha256':'same','physical_identity_sha256':'physical'},{}))
    monkeypatch.setattr(cli,'digest',lambda path:str(path))
    monkeypatch.setattr(cli,'hashes',lambda paths:{'candidate':'owned'})
    monkeypatch.setattr(cli,'qualified_reference',lambda *a:({},{}))
    monkeypatch.setattr(cli,'current_sources',lambda *a:None)
    monkeypatch.setattr(cli,'read_owned',lambda *a,**k:(dict(input_sha256='same',frozen_physical_identity_sha256='physical',
        candidate_source_sha256={'candidate':'owned'},reference_sha256='reference'),{}))
    report={'contract_sha256':'contract'}
    cli.qualify(args,{},report,{})
    assert not report['passed']
    assert not report['gates']['small_amplitude_time_domain_wave_responses']['passed']


def test_stale_validator_source_rejects_previously_passed_inertia(monkeypatch):
    args=SimpleNamespace(frozen=Path('frozen'),reference=Path('reference'),inertia=Path('inertia'),waves=Path('waves'))
    monkeypatch.setattr(cli,'frozen_input',lambda *a:({'input_sha256':'same','physical_identity_sha256':'physical'},{}))
    monkeypatch.setattr(cli,'digest',lambda path:str(path))
    monkeypatch.setattr(cli,'qualified_reference',lambda *a:({},{}))
    producer=hashes(CANDIDATE+VALIDATORS)
    producer['tests/validation/__main__.py']='stale'
    monkeypatch.setattr(cli,'read_owned',lambda *a,**k:(dict(input_sha256='same',frozen_physical_identity_sha256='physical',
        candidate_source_sha256=hashes(CANDIDATE),reference_sha256='reference',producer_sha256=producer),{}))
    report={'contract_sha256':'contract'}
    cli.qualify(args,{},report,{})
    assert not report['gates']['all_six_added_mass']['passed']
    assert 'source' in report['gates']['all_six_added_mass']['reason']
    assert 'tests/validation/__main__.py' in REFERENCE

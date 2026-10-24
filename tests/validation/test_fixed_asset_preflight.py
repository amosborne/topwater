import numpy as np
import pytest
import json
from pathlib import Path
from topwater import InitialState
from .fixed_asset_preflight import source_inventory,require_calm_state,recorded_green_identity


def test_finite_sources_do_not_imply_ainf_or_complete_diffraction_coverage():
    fields=('sources','potential','pressure','normal_velocity','source_potential','source_normal_velocity')
    arrays={f'level{level}_{name}':np.zeros((3,n,7),complex) for level,n in enumerate([2,3,4,5],1) for name in fields}
    arrays['raw_Ainf']=np.zeros((4,6,6))
    result=source_inventory(arrays,[2,3,4,5])
    assert not result['has_ainf_densities'] and result['missing_diffraction_frequency_heading_cases']==6
    assert result['radiation_frequency_cases']==18 and result['diffraction_frequency_heading_cases']==3
    arrays['level4_sources'][0,0,0]=np.nan
    with pytest.raises(ValueError):source_inventory(arrays,[2,3,4,5])


def test_asset_adoption_rejects_moving_fixture_even_if_both_states_match():
    calm=InitialState();moving=InitialState(linear_velocity=[.01,0,0])
    require_calm_state(calm,calm)
    with pytest.raises(AssertionError):require_calm_state(moving,calm)
    with pytest.raises(ValueError):require_calm_state(moving,moving)


def test_named_green_contract_uses_recorded_settings_and_rejects_mismatch():
    spec=json.loads(Path('validation/canonical-boat.json').read_text())['reference_numerics']
    settings=dict(finite={'tabulation_nr':676,'tabulation_nz':372},infinite={'tabulation_nr':0,'tabulation_nz':0})
    result=recorded_green_identity(spec,settings,settings['finite'])
    assert result['recorded']==settings and isinstance(result['finite_profile'],str)
    with pytest.raises(ValueError):recorded_green_identity(spec,settings,{'tabulation_nr':0,'tabulation_nz':0})

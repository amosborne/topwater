"""Isolated verification process: execute an immutable archived radiation RHS."""
import json
from pathlib import Path
import sys

archive=Path(sys.argv[1]).resolve()
sys.path[:0]=[str(archive/'src'),str(archive)]

import numpy as np
from topwater import SimulationInputs
from tests.experiments.free_surface.boat_radiation_time import Radiation

parent=json.loads(Path(sys.argv[2]).read_text())
inputs=SimulationInputs.from_json((archive/'inputs.json').read_text())
system=Radiation(inputs,parent['grid_size'],parent['box_length_m'],
                 parent.get('body_refinement',1),parent['sponge_strength_per_second'])
with np.load(Path(sys.argv[2]).with_suffix('.npz')) as saved:
    state=system.cp.asarray(saved['final_state'])
    time=float(saved['times'][-1])
rate,force,_,_=system.rhs(time,state,2*np.pi*parent['frequency_hz'],parent['ramp_seconds'])
loaded={name:str(Path(module.__file__).resolve()) for name,module in list(sys.modules.items())
        if (name in ('topwater','tests') or name.startswith('topwater.') or name.startswith('tests.'))
        and getattr(module,'__file__',None)}
if any(not Path(path).is_relative_to(archive) for path in loaded.values()):
    raise ValueError('archived reference imported a current producer')
np.savez(sys.argv[3],rate=system.cp.asnumpy(rate),force=system.cp.asnumpy(force),
         panels=system.panels.triangles,exterior=system.exterior,hidden=system.hidden)
Path(sys.argv[3]).with_suffix('.json').write_text(json.dumps(loaded,indent=2))

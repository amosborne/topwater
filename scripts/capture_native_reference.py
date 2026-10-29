"""Capture the unchanged WEC-Sim example for MATLAB-free regression tests.

Run with ``uv run --group wecsim python scripts/capture_native_reference.py``.
Requires the existing shared MATLAB engine. No public example is overwritten.
"""
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy.io import loadmat

from topwater.solve.wecsim import engine as backend
from topwater.solve.preparation import prepare_case, waves_from_sea
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.workflow_data import Boat, Sea


def main():
    root = Path(__file__).resolve().parents[1]
    case = root/'examples/simple-boat-wecsim'
    destination = root/'tests/data/simple_boat_wecsim.npz'
    if destination.exists():
        raise FileExistsError(f'Reference already exists: {destination}')
    fork = backend.preflight()
    from topwater.solve.wecsim import exchange as producer
    setup = prepare_case(Boat.load(case/'boat.npz').to_body(),
                         waves_from_sea(Sea.load(case/'sea.npz')), 100., rest_s=10., ramp_s=5.)
    setup.train.certify_incident_startup(setup)
    with tempfile.TemporaryDirectory(prefix='topwater-reference-') as scratch:
        workspace = Path(scratch)
        last = max(34, int(np.ceil(setup.train.omega[-1]*10)))
        h = setup.hydrostatic
        hydro = producer.build_hydrodynamics(
            setup.wet_triangles_m, setup.mean_cg_m, setup.mass_matrix, workspace/'hydro',
            density=WATER_DENSITY, gravity=GRAVITY, depth=np.inf,
            hydrostatic_stiffness=h.stiffness, cb=h.center_of_buoyancy,
            volume=h.displaced_volume, excitation_omega=np.arange(1, last+1)/10,
            radiation_omega=np.arange(1, max(120, last)+1)/10,
            green_cache_dir=Path.home()/'.cache/topwater/capytaine-3.0.0')
        motion = backend.integrate_motion(setup, hydro, workspace/'body', fork)
        # The guard restores engine settings; numeric run variables remain available.
        import matlab.engine
        engine = matlab.engine.connect_matlab(backend.ENGINE_NAME)
        try:
            path = (workspace/'audit.mat').as_posix().replace("'", "''")
            engine.eval(f"save('{path}','x','solid','hf','-v7')", nargout=0)
        finally:
            del engine
        audit = loadmat(workspace/'audit.mat', simplify_cells=True)
        metadata = dict(fork_commit=backend.SUPPORTED_FORK_COMMIT, dwell_s=100.,
                        rest_s=10., ramp_s=5., dt_s=.01,
                        inputs={name: hashlib.sha256((case/name).read_bytes()).hexdigest()
                                for name in ('boat.npz', 'sea.npz')})
        data = {name: getattr(motion, name) for name in motion.__dataclass_fields__}
        data.update({name: getattr(hydro, name) for name in
                     ('radiation_omega', 'radiation_damping', 'excitation_omega',
                      'excitation_headings', 'excitation', 'Ainf', 'M', 'K', 'density', 'gravity')})
        data.update(force_components=audit['x'][:, 18:], solid=audit['solid'],
                    residual_added_mass=audit['hf']['fAddedMass'],
                    radiation_kernel=audit['hf']['irkb'],
                    metadata_json=json.dumps(metadata))
        np.savez_compressed(destination, **data)
        print(f'Reference captured: {destination} ({destination.stat().st_size} bytes)', flush=True)


if __name__ == '__main__':
    main()

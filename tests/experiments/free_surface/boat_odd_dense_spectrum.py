"""Full finite-dimensional odd physical-trace spectrum, verification only.

Assemble D=M^-1(Kphysical-F*JD^-1*N) with diagonal-mass similarity. No
symmetrization, modal filtering, new transient, or runtime dense-operator path.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from .regularized_snapshot import RegularizedSnapshot, owned
from .odd_spectrum_engine import run_spectrum




def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    odd, body_probe = owned(args.odd), owned(args.body_probe)
    parents = {name: sha256(getattr(args, name).read_bytes()).hexdigest()
               for name in ('water', 'phase', 'mode', 'body_probe')}
    if odd['parents'] != parents or body_probe['parents'] != {k: parents[k] for k in ('water', 'phase', 'mode')}:
        raise ValueError('odd/body-probe/snapshot ownership differs')
    for name in ('regularized_snapshot.py', 'physical_p1_patch.py', 'compact_p1_trace.py',
                 'compact_p1_maps.py', 'boat_spatial_work.py', 'panel_gradient_reference.py'):
        source = Path('tests/experiments/free_surface')/name
        if sha256(source.read_bytes()).hexdigest() != odd['producer_sha256'].get(str(source)):
            raise ValueError('odd physical action source changed: '+name)
    with np.load(args.odd.with_suffix('.npz')) as z:
        F4, F8 = z['q4_F'].copy(), z['q8_F'].copy()
        jd, psi = z['jd'].copy(), z['physical_psi'].copy()
        saved_rate = z['q8_eta_rate'].copy()
    mass, nmap = snapshot.patch.mass, snapshot.coarse['n']
    count = len(psi)
    if (F4.shape != F8.shape or F8.shape != (count, len(jd))
            or jd.shape != (len(jd), len(jd)) or not np.isfinite(F8).all()):
        raise ValueError('saved odd operator dimensions/values invalid')
    if sha256((mass@psi).tobytes()).hexdigest() != odd['physical_psi_moment_sha256']:
        raise ValueError('odd physical state ownership differs')
    run_spectrum(snapshot, nmap, jd, F4, F8, psi, args, report, reference_rate=saved_rate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'odd', 'body-probe', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--batch-columns', type=int, default=64)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.batch_columns, args.maximum_bytes) <= 0:
        parser.error('fresh output/source archive and positive bounds required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True)
    hashes = {}
    for source in [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'), *Path(__file__).parent.glob('*.py')]:
        source_relative = source.resolve().relative_to(Path.cwd())
        target = archive/source_relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(source_relative)] = sha256(source.read_bytes()).hexdigest()
    report = dict(complete=False, producer_sha256=hashes, parents={name: sha256(getattr(args, name).read_bytes()).hexdigest()
        for name in ('water', 'phase', 'mode', 'odd', 'body_probe')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))

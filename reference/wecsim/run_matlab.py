"""Run one assigned script in the existing reference engine; never start an engine."""
import argparse
import io
import json
from pathlib import Path
import time
import matlab.engine

parser = argparse.ArgumentParser()
parser.add_argument('script', type=Path)
parser.add_argument('--log', type=Path, required=True)
parser.add_argument('--wait-for-engine', type=float, default=0)
args = parser.parse_args()
assert args.script.is_file()
assert not args.log.exists()
deadline = time.time() + args.wait_for_engine
while 'CAPYTAINE_REFERENCE' not in matlab.engine.find_matlab():
    if time.time() >= deadline:
        raise RuntimeError('Existing CAPYTAINE_REFERENCE engine is not visible')
    time.sleep(0.5)
Path('artifacts/reference/.engine-last-used').touch()
engine = matlab.engine.connect_matlab('CAPYTAINE_REFERENCE')
started = time.time()
out, err = io.StringIO(), io.StringIO()
try:
    print(json.dumps({'engine_pid':engine.eval("feature('getpid')",nargout=1),
                      'script':str(args.script.resolve()),'started_unix':started}), flush=True)
    script = str(args.script.resolve()).replace('\\','/').replace("'", "''")
    if args.script.read_text(encoding='utf-8').lstrip().startswith('function '):
        folder = str(args.script.resolve().parent).replace('\\','/').replace("'", "''")
        engine.eval(f"addpath('{folder}'); {args.script.stem}()", nargout=0, stdout=out, stderr=err)
    else:
        engine.eval(f"run('{script}')", nargout=0, stdout=out, stderr=err)
except BaseException:
    args.log.write_text(out.getvalue() + '\nSTDERR\n' + err.getvalue(),encoding='utf-8')
    print(out.getvalue(), flush=True)
    raise
else:
    args.log.write_text(out.getvalue() + '\nSTDERR\n' + err.getvalue(),encoding='utf-8')
    print(out.getvalue(), flush=True)
    print(json.dumps({'complete':True,'elapsed_s':time.time()-started}),flush=True)
finally:
    engine.exit()

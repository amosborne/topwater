"""Allocated sequential H5 readback then harmonic-only diagnostics; no dynamics."""
import subprocess
import sys

commands = [
    [sys.executable, 'reference/wecsim/run_matlab.py', 'reference/wecsim/check_generated_h5.m',
     '--log', 'artifacts/reference/generated-bank-consumer-20261024/h5-readback-attempt1.log'],
    [sys.executable, 'reference/wecsim/screen_generated_banks.py', '--output',
     'artifacts/reference/generated-bank-consumer-20261024/harmonic-attempt1']]
for command in commands:
    process = subprocess.Popen(command)
    print(f'Owned sequential child {process.pid}: {command}', flush=True)
    code = process.wait()
    print(f'Exact sequential child {process.pid} terminal exit{code}', flush=True)
    if code != 0:
        raise SystemExit(code)

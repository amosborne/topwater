"""Fail-closed execution bookkeeping for the allocated geometry job only.

Claim, verify and run are separate invocations. Only run can spawn the numeric
producer, and it requires the receipt from the successful separate verification.
Completion/idle readback/lease release remain the execution owner's responsibility.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import subprocess
import sys
from importlib.metadata import distribution, version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT/'docs/experts/execution-state.json'
LEASE = ROOT/'artifacts/numerical.lease'


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def frozen_identity(proposal):
    for path, digest in proposal['input_sha256'].items():
        if sha(ROOT/path) != digest:
            raise ValueError('Frozen input identity changed: '+path)
    for name, expected in proposal['dependency_identity'].items():
        package = distribution(name)
        if version(name) != expected['version']:
            raise ValueError('Frozen dependency version changed: '+name)
        if not all(sha(package.locate_file(path)) == digest for path, digest in expected['file_sha256'].items()):
            raise ValueError('Frozen dependency implementation changed: '+name)


def claimed_ledger(owner, manifest):
    state = read_json(LEDGER)
    if not LEASE.is_file():
        raise ValueError('Numerical lease is absent')
    if state['state'] != 'claimed' or state['execution_owner'] != owner or state['manifest'] != manifest:
        raise ValueError('Claimed ledger owner/state/manifest does not match this allocation')
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['record-claim', 'verify', 'run'])
    parser.add_argument('--proposal', type=Path, required=True)
    parser.add_argument('--owner', required=True)
    args = parser.parse_args()
    proposal_path = (ROOT/args.proposal).resolve()
    manifest = proposal_path.relative_to(ROOT).as_posix()
    proposal = read_json(proposal_path)
    if proposal['owner'] != args.owner or not proposal['status'].startswith('ready_for_coordinator_execution_allocation'):
        raise ValueError('Proposal is not ready for this allocated owner')
    frozen_identity(proposal)
    receipt_path = proposal_path.parent/proposal['verification_receipt_name']
    if args.stage == 'record-claim':
        state = read_json(LEDGER)
        if not LEASE.is_file() or state['state'] != 'idle':
            raise ValueError('New exclusive lease and idle ledger are required')
        state.update(state='claimed', execution_owner=args.owner, assignment=proposal['assignment'],
            claimed_at=utc(), manifest=manifest,
            active_numerical_job=dict(state='not_yet_launched', host=platform.node(),
                output=proposal['output'], recording_process_pid=os.getpid()),
            note='Claim recorded only; separate successful verify invocation is required before run.')
        write_json(LEDGER, state)
        print('CLAIM_RECORDED: no numeric process spawned')
        return 0
    state = claimed_ledger(args.owner, manifest)
    if args.stage == 'verify':
        if receipt_path.exists():
            raise ValueError('Use a unique fresh verification receipt')
        write_json(receipt_path, dict(owner=args.owner, manifest=manifest,
            proposal_sha256=sha(proposal_path), claimed_at=state['claimed_at'], verified_at=utc(),
            verification_process_pid=os.getpid(), recording_process_pid=state['active_numerical_job']['recording_process_pid'],
            ledger_state='claimed', lease_exists=True, all_frozen_identities_match=True))
        print('CLAIM_READBACK_VERIFIED: owner/manifest/state/lease/frozen identities match; no numeric process spawned')
        return 0
    receipt = read_json(receipt_path)
    if not (receipt['owner'] == args.owner and receipt['manifest'] == manifest and
            receipt['claimed_at'] == state['claimed_at'] and receipt['ledger_state'] == 'claimed' and
            receipt['lease_exists'] and receipt['all_frozen_identities_match'] and
            receipt['proposal_sha256'] == sha(proposal_path)):
        raise ValueError('Separate verification receipt does not match the current claimed allocation')
    if receipt['verification_process_pid'] == os.getpid():
        raise ValueError('Verification must be a separate process invocation')
    command = [sys.executable, *proposal['producer_argv']]
    # This is the ONLY spawn site; every failure above exits before reaching it.
    process = subprocess.Popen(command, cwd=ROOT)
    ledger_error = None
    try:
        state['state'] = 'running'
        state['active_numerical_job'].update(state='running', producer_pid=process.pid,
            owner_launcher_pid=os.getpid(), launched_at=utc(), command=command,
            verification_receipt=str(receipt_path.relative_to(ROOT)))
        state['note'] = 'Numeric producer launched after separate successful claimed-ledger readback. Sources stay frozen until exact terminal evidence.'
        write_json(LEDGER, state)
    except BaseException as error:
        ledger_error = error
        print(json.dumps(dict(producer_pid=process.pid, ledger_write_failed=True,
            action='Waiting for the exact child terminal before reporting this failure')), flush=True)
    finally:
        # A bookkeeping failure after spawn must never leave an unmonitored
        # producer behind a failed parent command. Preserve the claim and wait.
        code = process.wait()
        print(json.dumps(dict(producer_pid=process.pid, exact_process_terminal=True, exit_code=code)), flush=True)
    if ledger_error is not None:
        raise ledger_error
    return code


if __name__ == '__main__':
    raise SystemExit(main())

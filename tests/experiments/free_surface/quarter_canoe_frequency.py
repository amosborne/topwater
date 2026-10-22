"""First actual uncompressed PML frequency control: omega2 EE/heave only.

This worker has a mandatory Windows9GiB per-process commit limit. It never
materializes factor.L/U merely to measure their arrays. A single computed
heave response is not a six-mode matrix or a radiation-family qualification.
No mesh, boundary, physical tolerance or PML parameter is fitted here.
"""
import argparse
import ctypes
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from . import p4_tetra
from .p4_canoe import loads
from .quarter_canoe_control import scalar_operator
from .quarter_canoe_preflight import digest, WIDTH, BUFFER
from .quarter_canoe import POCKET


COMMIT_CAP = 9*1024**3


class ProcessMemoryLimit:
    """One worker's self-job; assignment failure stops before numerical work.

    Microsoft documents this flag as limiting process COMMITTED virtual memory,
    not working set. Nested jobs may fail; an unbounded fallback is forbidden.
    https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-jobobject_extended_limit_information
    https://learn.microsoft.com/en-us/windows/win32/procthread/nested-jobs
    """
    def __init__(self):
        if sys.platform != 'win32':
            raise RuntimeError('this actual-candidate worker requires its Windows process-memory cap')
        from ctypes import wintypes

        class Basic(ctypes.Structure):
            _fields_ = [('PerProcessUserTimeLimit', ctypes.c_longlong),
                        ('PerJobUserTimeLimit', ctypes.c_longlong), ('LimitFlags', wintypes.DWORD),
                        ('MinimumWorkingSetSize', ctypes.c_size_t), ('MaximumWorkingSetSize', ctypes.c_size_t),
                        ('ActiveProcessLimit', wintypes.DWORD), ('Affinity', ctypes.c_size_t),
                        ('PriorityClass', wintypes.DWORD), ('SchedulingClass', wintypes.DWORD)]

        class IO(ctypes.Structure):
            _fields_ = [(name, ctypes.c_ulonglong) for name in
                ('ReadOperationCount', 'WriteOperationCount', 'OtherOperationCount',
                 'ReadTransferCount', 'WriteTransferCount', 'OtherTransferCount')]

        class Extended(ctypes.Structure):
            _fields_ = [('BasicLimitInformation', Basic), ('IoInfo', IO)]+[(name, ctypes.c_size_t)
                for name in ('ProcessMemoryLimit', 'JobMemoryLimit', 'PeakProcessMemoryUsed', 'PeakJobMemoryUsed')]

        class Counters(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD)]+[
                (name, ctypes.c_size_t) for name in ('PeakWorkingSetSize', 'WorkingSetSize',
                'QuotaPeakPagedPoolUsage', 'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage',
                'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage', 'PrivateUsage')]

        self.Extended, self.Counters = Extended, Counters
        self.kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        self.kernel.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
        self.kernel.CreateJobObjectW.restype = wintypes.HANDLE
        self.kernel.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
        self.kernel.SetInformationJobObject.restype = wintypes.BOOL
        self.kernel.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p,
                                                         wintypes.DWORD, ctypes.c_void_p]
        self.kernel.QueryInformationJobObject.restype = wintypes.BOOL
        self.kernel.GetCurrentProcess.restype = wintypes.HANDLE
        self.kernel.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        self.kernel.AssignProcessToJobObject.restype = wintypes.BOOL
        self.kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        self.kernel.CloseHandle.restype = wintypes.BOOL
        self.psapi = ctypes.WinDLL('psapi', use_last_error=True)
        self.psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        self.psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
        self.handle = self.kernel.CreateJobObjectW(None, None)
        if not self.handle:
            raise ctypes.WinError(ctypes.get_last_error())
        info = Extended(); info.BasicLimitInformation.LimitFlags = 0x100  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
        info.ProcessMemoryLimit = COMMIT_CAP
        if not self.kernel.SetInformationJobObject(self.handle, 9, ctypes.byref(info), ctypes.sizeof(info)):
            error = ctypes.get_last_error(); self.close(); raise ctypes.WinError(error)
        if not self.kernel.AssignProcessToJobObject(self.handle, self.kernel.GetCurrentProcess()):
            error = ctypes.get_last_error(); self.close(); raise ctypes.WinError(error)
        # Query back the enforced policy; never rely solely on requested values.
        self.metrics()

    def metrics(self):
        info = self.Extended()
        if not self.kernel.QueryInformationJobObject(self.handle, 9, ctypes.byref(info), ctypes.sizeof(info), None):
            raise ctypes.WinError(ctypes.get_last_error())
        if not info.BasicLimitInformation.LimitFlags & 0x100 or info.ProcessMemoryLimit != COMMIT_CAP:
            raise RuntimeError('OS did not retain the required process-memory limit')
        counters = self.Counters(); counters.cb = ctypes.sizeof(counters)
        if not self.psapi.GetProcessMemoryInfo(self.kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        return dict(commit_limit_bytes=COMMIT_CAP, job_peak_process_commit_bytes=int(info.PeakProcessMemoryUsed),
            private_commit_bytes=int(counters.PrivateUsage), peak_working_set_bytes=int(counters.PeakWorkingSetSize),
            working_set_bytes=int(counters.WorkingSetSize),
            policy='Windows JOB_OBJECT_LIMIT_PROCESS_MEMORY; no fallback if assignment fails')

    def close(self):
        if self.handle:
            self.kernel.CloseHandle(self.handle); self.handle = None


def storage(matrix):
    return int(matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)


def _input_contract(report):
    """Exact fixed-case and source identity; no approximate reuse criterion."""
    contract = {key: report[key] for key in ('mesh_report_sha256', 'mesh_npz_sha256',
        'capytaine_report_sha256', 'omega', 'density', 'gravity', 'parity', 'computed_modes',
        'pml_width', 'pml_strength', 'pml_onset', 'volume_quadrature_order', 'surface_quadrature_order')}
    # Writer/logging/ordering changes do not alter an already assembled matrix.
    # Numerical helpers and every argument supplied to scalar_operator remain
    # exact identities. The original worker hash is retained as provenance.
    contract['assembly_source_sha256'] = {name: value for name, value in report['source_sha256'].items()
                                         if name != 'quarter_canoe_frequency.py'}
    return contract


def _array_hash(array):
    array = np.ascontiguousarray(array)
    value = sha256(str(array.shape).encode()+array.dtype.str.encode())
    value.update(array.tobytes())
    return value.hexdigest()


def save_input(path, system, load, rhs, report):
    """One plain failed-factor retry artifact, not a general operator cache."""
    if path.exists() or path.with_suffix('.json').exists():
        raise FileExistsError('persisted input already exists; choose a new output or explicitly reuse it')
    arrays = dict(data=system.data, indices=system.indices, indptr=system.indptr,
                  shape=np.asarray(system.shape, dtype=np.int64), body_load=load, rhs=rhs)
    np.savez(path, **arrays)
    manifest = dict(contract=_input_contract(report),
        writer_source_sha256=report['source_sha256']['quarter_canoe_frequency.py'],
        array_sha256={name: _array_hash(value) for name, value in arrays.items()},
        npz_sha256=sha256(path.read_bytes()).hexdigest())
    path.with_suffix('.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


def reuse_input(path, report, expected_load):
    manifest = json.loads(path.with_suffix('.json').read_text())
    if manifest['contract'] != _input_contract(report):
        raise ValueError('retry input source, mesh, reference or fixed physical case differs')
    if sha256(path.read_bytes()).hexdigest() != manifest['npz_sha256']:
        raise ValueError('retry scalar input file hash differs')
    with np.load(path) as saved:
        arrays = {name: saved[name].copy() for name in
                  ('data', 'indices', 'indptr', 'shape', 'body_load', 'rhs')}
    for name, array in arrays.items():
        if _array_hash(array) != manifest['array_sha256'][name]:
            raise ValueError(f'retry scalar input {name} hash differs')
    if not np.array_equal(arrays['body_load'], expected_load):
        raise ValueError('retry input body load differs from fresh physical quadrature')
    if not np.array_equal(arrays['rhs'], expected_load[:, 2].astype(complex)):
        raise ValueError('retry input heave RHS differs from the fixed physical load')
    system = csc_matrix((arrays['data'], arrays['indices'], arrays['indptr']),
                        shape=tuple(arrays['shape']))
    system.check_format(full_check=True)
    if system.shape != (len(expected_load),)*2 or not np.all(np.isfinite(system.data)):
        raise ValueError('invalid retry scalar matrix')
    return system, arrays['rhs'], manifest


def run(mesh_path, reference_path, output, budget, *, input_path=None):
    from .canoe_wet_reference import pinned_reference, verify_coverage
    from .quarter_canoe_preflight import reflected_body
    from .quarter_canoe import make_plc

    if input_path is None:
        proposed = output.with_suffix('.input.npz')
        if proposed.exists() or proposed.with_suffix('.json').exists():
            raise FileExistsError('persisted input already exists; choose a new output or explicitly reuse it')
    meta = json.loads(mesh_path.with_suffix('.json').read_text())
    if (not meta.get('mesh_only') or meta['spacing_envelope'] != 16 or meta['wet_subdivisions'] != 1
            or meta['pml_width'] != [WIDTH, WIDTH] or meta['pml_strength'] != 8. or meta['real_terminal'] != 1.
            or meta['pml_onset'] != (POCKET+BUFFER).tolist()):
        raise ValueError('worker requires the fixed reviewed uncompressed n16/wet1 mesh preflight')
    with np.load(mesh_path) as saved:
        points, tets, body, parents, center = [saved[name].copy() for name in
            ('points', 'tets', 'body_faces', 'body_parents', 'center')]
    for name, values in [('linear_points', points), ('linear_tets', tets), ('body_faces', body), ('body_parent_ids', parents)]:
        if digest(values) != meta['array_sha256'][name]:
            raise ValueError(f'saved {name} differs from its actual preflight hash')
    mesh = p4_tetra.make_mesh(points, tets)
    if (len(mesh.points) != meta['quarter_p4_nodes'] or digest(mesh.points) != meta['array_sha256']['p4_master_points']
            or digest(mesh.elements) != meta['array_sha256']['p4_elements']):
        raise ValueError('rebuilt P4 masters differ from the actual mesh preflight')
    domain, reference, capytaine, provenance = pinned_reference(reference_path)
    if not np.array_equal(center, domain.rest.position):
        raise ValueError('saved torque center differs from pinned reference')
    # Regenerate only the cheap physical PLC parent-reflection permutations,
    # not another tetrahedral mesh or fluid operator.
    plc = make_plc(meta['x'], meta['y'], meta['z'], subdivisions=1)
    candidate, candidate_parents = reflected_body(mesh.points[body][:, ::-1], parents, plc.parent_reflections)
    coverage = verify_coverage(candidate, candidate_parents, reference, domain.contact.body_source_faces)
    report = dict(qualified=False, stage='mesh and physical-reference hashes verified',
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('quarter_canoe_frequency.py', 'quarter_canoe_control.py', 'quarter_canoe.py',
                         'quarter_canoe_preflight.py', 'canoe_wet_reference.py', 'pml_corner.py',
                         'pml_strip.py', 'p4_tetra.py', 'p4_canoe.py')},
        mesh_report_sha256=sha256(mesh_path.with_suffix('.json').read_bytes()).hexdigest(),
        mesh_npz_sha256=sha256(mesh_path.read_bytes()).hexdigest(),
        capytaine_report_path=str(reference_path.resolve()),
        capytaine_report_sha256=sha256(reference_path.read_bytes()).hexdigest(),
        reference=provenance, physical_coverage=coverage, memory=budget.metrics(),
        omega=2., density=1000., gravity=9.81, parity='EE', computed_modes=['heave'],
        pml_width=[WIDTH, WIDTH], pml_strength=8., pml_onset=(POCKET+BUFFER).tolist(),
        volume_quadrature_order=12, surface_quadrature_order=12, lu_ordering='MMD_AT_PLUS_A',
        scope='One omega2 EE/heave coefficient from the fixed uncompressed PML candidate. '
              'No uncomputed matrix entries, full-family gate, convergence or nonlinear qualification.')

    def checkpoint(stage):
        report['stage'] = stage; report['memory'] = budget.metrics()
        output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(dict(stage=stage, memory=report['memory'])), flush=True)

    load, _ = loads(mesh, body, center)
    if input_path is None:
        checkpoint('before scalar operator assembly')
        if report['volume_quadrature_order'] != report['surface_quadrature_order']:
            raise ValueError('the shared scalar assembler uses the same surface/volume quadrature order')
        matrix = scalar_operator(mesh, report['omega'], report['pml_onset'], report['pml_width'],
            report['pml_strength'], order=report['volume_quadrature_order'], max_nodes=120000)
        # EE needs no plane Dirichlet constraints. The complex FS Robin term
        # fixes the potential gauge. No physical matrix values change below.
        system = matrix.tocsc(); del matrix
        rhs = np.asarray(load[:, 2], complex)
    else:
        checkpoint('validating persisted scalar input; no operator assembly')
        system, rhs, _ = reuse_input(input_path, report, load)
    report['input'] = dict(rows=system.shape[0], nonzeros=system.nnz, sparse_array_bytes=storage(system),
                          rhs_array_bytes=rhs.nbytes, storage_is_not_lu_fill_bound=True)
    if system.shape[0] > 120000 or system.nnz > 16000000 or storage(system) > 400000000:
        checkpoint('input preflight rejected before LU')
        raise RuntimeError('actual sparse input exceeds the declared pre-LU cap')
    if input_path is None:
        input_path = output.with_suffix('.input.npz')
        manifest = save_input(input_path, system, load, rhs, report)
    else:
        manifest = json.loads(input_path.with_suffix('.json').read_text())
    report['persisted_input'] = dict(path=str(input_path.resolve()), npz_sha256=manifest['npz_sha256'],
        manifest_sha256=sha256(input_path.with_suffix('.json').read_bytes()).hexdigest())
    checkpoint('actual sparse input measured; before LU')
    factor = None
    try:
        factor = splu(system, permc_spec='MMD_AT_PLUS_A')
        phi = factor.solve(rhs)
    except MemoryError as error:
        factor = None
        report['factor_failure'] = dict(type=type(error).__name__, message=str(error))
        checkpoint('LU MemoryError under the unchanged process commit cap')
        raise
    response = 4*load[:, 2]@phi
    added = report['density']*response.real
    damping = report['density']*report['omega']*response.imag
    exact_added = float(capytaine[0]['added_mass'][2][2])
    exact_damping = float(capytaine[0]['radiation_damping'][2][2])
    report['computed_heave'] = dict(raw_response_real=float(response.real), raw_response_imag=float(response.imag),
        added_mass=float(added), radiation_damping=float(damping),
        capytaine_added_mass=exact_added, capytaine_radiation_damping=exact_damping,
        added_mass_relative_difference=float((added-exact_added)/exact_added),
        radiation_damping_relative_difference=float((damping-exact_damping)/exact_damping),
        relative_linear_residual=float(np.linalg.norm(system@phi-rhs)/np.linalg.norm(rhs)),
        factor_nonzeros=int(factor.nnz), factor_array_bytes_not_materialized=True,
        exp_time_convention='exp(-i omega t)',
        fourfold_work_factor=4, full_family_gate_evaluated=False)
    # Retain the physical nodal solution for independent pressure/work checks.
    # Release LU before serialization; measuring its L/U arrays is unnecessary.
    del factor
    potential_path = output.with_suffix('.npz')
    np.savez(potential_path, potential=phi, mode_index=np.array(2), omega=np.array(2.))
    report['potential'] = dict(path=str(potential_path.resolve()),
        array_sha256=sha256(np.asarray(phi, dtype='<c16').tobytes()).hexdigest(),
        npz_sha256=sha256(potential_path.read_bytes()).hexdigest(),
        dtype='<c16', ordering='Exact quarter P4 master ordering pinned by the mesh report hashes')
    checkpoint('computed omega2 EE/heave only')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True, help='actual n16 preflight NPZ, with sibling JSON')
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('artifacts/quarter-canoe-omega2-heave-mmd.json'))
    parser.add_argument('--reuse-input', type=Path, help='exact .input.npz plus sibling .input.json from this source/case')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    budget = ProcessMemoryLimit()  # Mandatory and queried back before rebuilding the mesh.
    try:
        with threadpool_limits(1):
            run(args.mesh, args.reference, args.output, budget, input_path=args.reuse_input)
    finally:
        budget.close()

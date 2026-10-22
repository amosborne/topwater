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


def run(mesh_path, reference_path, output, budget):
    from .canoe_wet_reference import pinned_reference, verify_coverage
    from .quarter_canoe_preflight import reflected_body
    from .quarter_canoe import make_plc

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
        volume_quadrature_order=12, surface_quadrature_order=12,
        scope='One omega2 EE/heave coefficient from the fixed uncompressed PML candidate. '
              'No uncomputed matrix entries, full-family gate, convergence or nonlinear qualification.')

    def checkpoint(stage):
        report['stage'] = stage; report['memory'] = budget.metrics()
        output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(dict(stage=stage, memory=report['memory'])), flush=True)

    checkpoint('before scalar operator assembly')
    matrix = scalar_operator(mesh, 2., POCKET+BUFFER, [WIDTH, WIDTH], 8., order=12, max_nodes=120000)
    load, _ = loads(mesh, body, center)
    # EE needs no plane Dirichlet constraints. This is exactly the all-unknown
    # scalar PML problem; the complex FS Robin term fixes the potential gauge.
    system = matrix.tocsc(); del matrix
    rhs = np.asarray(load[:, 2], complex)
    report['input'] = dict(rows=system.shape[0], nonzeros=system.nnz, sparse_array_bytes=storage(system),
                          rhs_array_bytes=rhs.nbytes, storage_is_not_lu_fill_bound=True)
    if system.shape[0] > 120000 or system.nnz > 16000000 or storage(system) > 400000000:
        checkpoint('input preflight rejected before LU')
        raise RuntimeError('actual sparse input exceeds the declared pre-LU cap')
    checkpoint('actual sparse input measured; before LU')
    factor = splu(system)
    phi = factor.solve(rhs)
    response = 4*load[:, 2]@phi
    added = 1000*response.real; damping = 2000*response.imag
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
    checkpoint('computed omega2 EE/heave only')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True, help='actual n16 preflight NPZ, with sibling JSON')
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('artifacts/quarter-canoe-omega2-heave.json'))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    budget = ProcessMemoryLimit()  # Mandatory and queried back before rebuilding the mesh.
    try:
        with threadpool_limits(1):
            run(args.mesh, args.reference, args.output, budget)
    finally:
        budget.close()

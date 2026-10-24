"""Owned validation artifacts and the predeclared entrywise criteria."""
from hashlib import sha256
import json
from pathlib import Path
import platform
import shutil
import subprocess
from importlib.metadata import version
import numpy as np

PHYSICAL = tuple('src/topwater/'+name+'.py' for name in
    ('__init__','inputs','_validation','frames','hydrostatics','_clipping','constants'))+('tests/reference_boat.py',)
CANDIDATE = PHYSICAL+tuple('src/topwater/'+name+'.py' for name in
    ('inertia','bem','boundary','rankine','_rankine_kernel','_rankine_extended','numerics'))
VALIDATORS = tuple('tests/validation/'+name+'.py' for name in ('__init__','__main__','common','reference'))
REFERENCE = PHYSICAL+('tests/capytaine_reference.py',)+VALIDATORS+('pyproject.toml','uv.lock')


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def hashes(paths):
    return {name:digest(name) for name in paths}


def physical_identity(report,arrays):
    """Stable physical binding; execution metadata is provenance only."""
    payload=dict(contract_sha256=report['contract_sha256'],input_sha256=report['input_sha256'],
        arrays={name:dict(dtype=np.asarray(value).dtype.str,shape=list(np.asarray(value).shape),
            sha256=sha256(np.ascontiguousarray(value).tobytes()).hexdigest())
            for name,value in sorted(arrays.items())})
    return sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def load_contract(path):
    value=json.loads(Path(path).read_text())
    p=value['physical_problem'];r=value['reference_numerics'];c=value['candidate_numerics']
    if (value['format_version']!=1 or p['fixture']!='tests.reference_boat'
            or p['body_function']!='reference_body' or p['pose_function']!='equilibrium'
            or p['water_depth']!='infinite' or p['forward_speed_m_s']!=0
            or p['initial_body_velocity']!='zero'
            or p['time_convention']!='Re(amplitude*exp(-i*omega*t))'
            or p['wave_direction']!='propagation, radians CCW from world+x'
            or p['incident_phase_origin_world_xy_m']!=[0.,0.]
            or value['implementation']['entry_point']!='topwater.inertia.instantaneous_added_mass'
            or p['mode_order']!=['Surge','Sway','Heave','Roll','Pitch','Yaw']):
        raise ValueError('unsupported canonical physical/implementation contract')
    for levels,minimum in ((r['body_subdivision_levels'],3),(c['body_subdivision_levels'],2)):
        if len(levels)<minimum or levels!=sorted(set(levels)) or any(type(n)!=int or n<0 for n in levels):
            raise ValueError('declared ordered nested reference/candidate levels required')
    if r['cpu_threads']!=1 or c['cpu_threads']!=1 or r['richardson_extrapolation'] is not False:
        raise ValueError('unsupported thread/extrapolation policy')
    cases=value['cases'];a=value['acceptance']
    for name in ('density_kg_m3','gravity_m_s2','incident_elevation_amplitude_m'):
        if not np.isfinite(p[name]) or p[name]<=0:
            raise ValueError('finite positive physical constants and incident amplitude required')
    for name in ('maximum_elapsed_seconds','maximum_estimated_workspace_bytes'):
        if isinstance(r[name],bool) or not np.isfinite(r[name]) or r[name]<=0:
            raise ValueError('finite positive reference resource ceilings required')
    if (not cases['frequencies_hz'] or not cases['propagation_headings_rad']
            or not np.isfinite(cases['frequencies_hz']).all() or min(cases['frequencies_hz'])<=0
            or max(cases['frequencies_hz'])>.6
            or len(set(cases['frequencies_hz']))!=len(cases['frequencies_hz'])
            or not np.isfinite(cases['propagation_headings_rad']).all()
            or len(set(cases['propagation_headings_rad']))!=len(cases['propagation_headings_rad'])
            or len(set(np.remainder(cases['propagation_headings_rad'],2*np.pi)))!=len(cases['propagation_headings_rad'])
            or a['reference_required_finest_increments']!=2
            or set(a['required_gates'])!={'reference','all_six_added_mass','small_amplitude_time_domain_wave_responses'}):
        raise ValueError('missing finite cases or required gate/uncertainty policy')
    for name in ('reference_mesh_relative','reference_mesh_absolute_normalized','candidate_relative','candidate_absolute_normalized'):
        if not np.isfinite(a[name]) or a[name]<=0:
            raise ValueError('positive finite predeclared acceptance criteria required')
    return value


def case_envelope(contract,beam):
    if not np.isfinite(beam) or beam<=0:
        raise ValueError('finite positive actual full physical hull beam required')
    f=np.asarray(contract['cases']['frequencies_hz'])
    wavelength=contract['physical_problem']['gravity_m_s2']/(2*np.pi*f*f)
    if np.any(f>.6) or np.any(wavelength<2*beam):
        raise ValueError('declared waves violate f<=0.6 Hz or deepwater wavelength>=2*actual boat beam')
    return dict(beam_m=float(beam),wavelengths_m=wavelength.tolist(),minimum_wavelength_m=float(2*beam))


def prepare(output,contract_path,kind):
    output=Path(output)
    if any(p.exists() for p in (output,output.with_suffix('.npz'),output.with_suffix('.source'))):
        raise FileExistsError('preserve checkpoints: output and archive must be fresh')
    output.parent.mkdir(parents=True,exist_ok=True)
    archive=output.with_suffix('.source');archive.mkdir()
    paths=[*Path('src/topwater').rglob('*.py'),*Path('tests/validation').rglob('*.py'),
        Path('tests/reference_boat.py'),Path('tests/capytaine_reference.py'),Path('pyproject.toml'),Path('uv.lock')]
    manifest={}
    for path in paths:
        relative=path.resolve().relative_to(Path.cwd()).as_posix()
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target);manifest[relative]=digest(path)
    shutil.copyfile(contract_path,archive/'contract.json')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    return dict(kind=kind,complete=False,passed=False,contract_sha256=digest(contract_path),
        current_commit=commit,producer_sha256=manifest,python=platform.python_version(),
        numpy=np.__version__,scipy=version('scipy'),
        implementation_source_sha256=hashes(CANDIDATE))


def save(output,report,arrays):
    output=Path(output);np.savez(output.with_suffix('.npz'),**arrays)
    report['arrays_sha256']=digest(output.with_suffix('.npz'))
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


def read_owned(path,contract_hash,kind,*,passed=False):
    path=Path(path);report=json.loads(path.read_text())
    if (report['kind']!=kind or report['contract_sha256']!=contract_hash
            or report['complete'] is not True or (passed and report['passed'] is not True)
            or report['arrays_sha256']!=digest(path.with_suffix('.npz'))):
        raise ValueError('missing, unqualified, stale or wrong-kind validation artifact: '+str(path))
    archive=path.with_suffix('.source')
    if digest(archive/'contract.json')!=contract_hash:
        raise ValueError('archived contract identity differs')
    for name,expected in report['producer_sha256'].items():
        if digest(archive/name)!=expected:
            raise ValueError('corrupted archived producer: '+name)
    with np.load(path.with_suffix('.npz')) as z:
        arrays={name:z[name].copy() for name in z.files}
    return report,arrays


def current_sources(report,paths):
    if any(report['producer_sha256'].get(name)!=digest(name) for name in paths):
        raise ValueError('current source identity differs from required producer graph')


def midpoint_children(triangles):
    a,b,c=triangles.transpose(1,0,2);ab,bc,ca=(a+b)/2,(b+c)/2,(c+a)/2
    return np.concatenate([np.stack(t,axis=1) for t in ((a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca))])


def triangles_at(base,level):
    for _ in range(level):
        base=midpoint_children(base)
    return base


def study_preflight(base_count,levels,limit):
    # Twelve dense complex128 arrays plus geometry/solver overhead. An estimate,
    # not an OS memory limit. NO declared finest mesh can be silently omitted.
    rows=[dict(level=n,panels=base_count*4**n,
        estimated_workspace_bytes=12*16*(base_count*4**n)**2+32768*(base_count*4**n)+64*1024**2) for n in levels]
    if any(row['estimated_workspace_bytes']>limit for row in rows):
        raise MemoryError('declared study exceeds workspace budget; revise contract, do not omit finest mesh')
    return rows


def entrywise(candidate,reference,uncertainty,atol,rtol):
    candidate,reference,uncertainty=map(np.asarray,(candidate,reference,uncertainty))
    if (candidate.shape!=reference.shape or uncertainty.shape!=reference.shape
            or not all(np.isfinite(a).all() for a in (candidate,reference,uncertainty))
            or np.any(uncertainty<0)):
        raise ValueError('complete finite identity-matched entry arrays required')
    error=abs(candidate-reference)+uncertainty;allowance=atol+rtol*abs(reference)
    return dict(passed=bool(np.all(error<=allowance)),maximum_tolerance_ratio=float(np.max(error/allowance)))


def mesh_quality(values,acceptance):
    values=np.asarray(values)
    if len(values)<3 or not np.isfinite(values).all():
        raise ValueError('three or more complete finite nested levels required')
    previous,last=abs(values[-2]-values[-3]),abs(values[-1]-values[-2])
    atol=acceptance['reference_mesh_absolute_normalized'];rtol=acceptance['reference_mesh_relative']
    allowance=atol+rtol*abs(values[-1])
    passed=bool(np.all(previous<=allowance)&np.all(last<=allowance)&np.all(last<=previous+atol))
    return np.maximum(previous,last),dict(passed=passed,
        previous_increment_max_ratio=float(np.max(previous/allowance)),
        latest_increment_max_ratio=float(np.max(last/allowance)),
        trend_passed=bool(np.all(last<=previous+atol)))

"""Six-motion acceleration screen from OWNED physical-psi phase operators.

No trajectory or work certificate. Matrices remain raw; symmetric parts are
used only to report diagnostic eigenvalues. No matrix repair or extra force.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.linalg import block_diag,lu_factor,lu_solve
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

from topwater import SimulationInputs
from topwater.constants import WATER_DENSITY
from topwater.frames import rotation_matrix
from topwater.inertia import instantaneous_added_mass
from topwater.rankine import TrianglePanels
from topwater.bem import assemble_dirichlet_halfspace
from .boat_spatial_work import source_matrix,subdivide


def run(args,report):
    phase=json.loads(args.phase.read_text())
    if not phase['complete'] or sha256(args.phase.with_suffix('.npz').read_bytes()).hexdigest()!=phase['arrays_sha256']:
        raise ValueError('complete owned phase arrays required')
    contact_path=args.phase.with_suffix('.source')/'contact.json'
    if sha256(contact_path.read_bytes()).hexdigest()!=phase['parent_contact_sha256']:
        raise ValueError('phase contact ownership changed')
    contact=json.loads(contact_path.read_text());audit_path=Path(contact['parents']['audit']['path'])
    if sha256(audit_path.read_bytes()).hexdigest()!=contact['parents']['audit']['sha256']:
        raise ValueError('reference audit ownership changed')
    audit=json.loads(audit_path.read_text())
    if sha256(audit_path.with_suffix('.npz').read_bytes()).hexdigest()!=audit['arrays_sha256']:
        raise ValueError('reference body arrays changed')
    inputs=SimulationInputs.from_json((audit_path.with_suffix('.source')/'inputs.json').read_text())
    if sha256(inputs.to_json().encode()).hexdigest()!=phase['input_sha256']:
        raise ValueError('physical input ownership changed')
    with np.load(audit_path.with_suffix('.npz')) as z: triangles=z['triangles'].copy()
    if sha256(triangles.tobytes()).hexdigest()!=phase['panel_sha256']:
        raise ValueError('physical body panel ownership changed')
    panels=TrianglePanels(triangles,-1);cg=inputs.initial_state.position
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    load=(panels.areas[:,None]*modes).T
    rotation=rotation_matrix(inputs.initial_state.orientation)
    dry=block_diag(np.eye(3)*inputs.body.mass,rotation@inputs.body.inertia_cg@rotation.T)
    values,vectors=np.linalg.eigh(dry)
    inverse_root=(vectors/np.sqrt(values))@vectors.T
    independent=instantaneous_added_mass(panels,cg,density=WATER_DENSITY)
    source_limit=WATER_DENSITY*load@assemble_dirichlet_halfspace(panels).single_layer@lu_solve(lu_factor(source_matrix(panels)),modes)
    refinement_rows=[dict(level=0,panels=len(panels),direct_diagonal=np.diag(independent.added_mass).tolist(),
        source_diagonal=np.diag(source_limit).tolist())]
    refined=triangles
    for level in (1,2):
        print(json.dumps(dict(body_refinement=level,start=True)),flush=True)
        refined=subdivide(refined);rp=TrianglePanels(refined,-1)
        rmodes=np.column_stack((rp.normals,np.cross(rp.centroids-cg,rp.normals)))
        direct=instantaneous_added_mass(rp,cg,density=WATER_DENSITY).added_mass
        indirect=WATER_DENSITY*(rp.areas[:,None]*rmodes).T@assemble_dirichlet_halfspace(rp).single_layer@lu_solve(
            lu_factor(source_matrix(rp)),rmodes)
        refinement_rows.append(dict(level=level,panels=len(rp),direct_diagonal=np.diag(direct).tolist(),
            source_diagonal=np.diag(indirect).tolist()))
        print(json.dumps(refinement_rows[-1]),flush=True)
    rows=[];saved=dict(dry_mass=dry,halfspace_reference=independent.added_mass,source_halfspace_limit=source_limit)
    with np.load(args.phase.with_suffix('.npz')) as z:
        p=z['pmap'];n=z['nmap'];s=z['even'];j=z['neumann'];count=len(z['active_nodes'])
        mass=csr_matrix((z['mass_data'],z['mass_indices'],z['mass_indptr']),shape=(count,count))
        mass_factor=splu(mass.tocsc())
        for record in phase['orders']:
            order=record['order'];cross=z[f'q{order}_cross']
            inverse_cross=mass_factor.solve(cross)
            schur=j-n@inverse_cross
            response=lu_solve(lu_factor(schur),modes)
            added=WATER_DENSITY*load@(s-p@inverse_cross)@response
            normalized=inverse_root@added@inverse_root
            norm=np.linalg.norm(normalized)
            row=dict(order=order,added_mass_raw=added.tolist(),
                dry_mass_normalized_reciprocity_defect=float(np.linalg.norm(normalized-normalized.T)/norm),
                dry_mass_normalized_symmetric_part_eigenvalues=np.linalg.eigvalsh((normalized+normalized.T)/2).tolist(),
                symmetric_part_eigenvalues_raw_units=np.linalg.eigvalsh((added+added.T)/2).tolist(),
                actual_dry_plus_added_condition_number=float(np.linalg.cond(dry+added)),
                mass_normalized_dry_plus_added_condition_number=float(np.linalg.cond(np.eye(6)+normalized)),
                independent_halfspace_principal_relative_difference=((np.diag(added)-np.diag(independent.added_mass))
                    /np.diag(independent.added_mass)).tolist(),
                source_halfspace_principal_relative_difference=((np.diag(added)-np.diag(source_limit))/np.diag(source_limit)).tolist())
            rows.append(row);saved[f'q{order}_added_mass']=added
    np.savez(args.output.with_suffix('.npz'),**saved)
    report.update(complete=True,qualified=False,input_sha256=phase['input_sha256'],panel_sha256=phase['panel_sha256'],
        rows=rows,halfspace_reference_diagnostics=independent.diagnostics,
        body_only_halfspace_refinement=refinement_rows,
        source_vs_direct_halfspace_principal_relative_difference=((np.diag(source_limit)-np.diag(independent.added_mass))
            /np.diag(independent.added_mass)).tolist(),
        arrays_sha256=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest(),scope=__doc__,
        limits=['Mixed physical-psi local closure with prescribed far chi is not the full-plane zero-potential boundary problem',
            'Independent halfspace result is a useful limit with different boundary closure and P0 formulation, not an exact oracle for this discretization',
            'Reciprocity defects can include collocation, finite spaces and far closure; no single cause is established',
            'No global state energy balance, evolving eta, nonlinear6DOF trajectory or speed qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('phase','output'): parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists(): parser.error('fresh output/source archive required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parent_phase_sha256=sha256(args.phase.read_bytes()).hexdigest())
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))

"""ONE geometry-only near/far count on the archived 472-face 0/10/20 ms chart.

No CUDA, fluid assembly, Taylor cache, factorization, threshold sweep or pilot.
Separation is bounded for supplied endpoint geometries and their piecewise
linear NODAL chords. The actual nonlinear ALE path between samples is not
enclosed. Contact/slender and opening-ratio choices are screening heuristics.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from topwater import InitialState
from .boat_geometry_template import evaluate_template
from .boat_world_field import refine_boundary

OPENING_FACTOR=2.0
SLENDER_QUALITY=1e-3
OBSERVED_RHS=64
OBSERVED_COLD_BUILDS=29
USES_PER_COLD_BUILD=OBSERVED_RHS/OBSERVED_COLD_BUILDS  # Includes startup/older-chart reuse, not anchor lifetime.


def _norm_bounds(a,b):
    """Outward float64 bounds on distances of the supplied coordinate arrays.

    Basic operations are widened individually; centers are arbitrary float
    points, not claims of exact geometric centroids. Reject nonfinite inputs.
    These do not enclose the exact CSG geometry behind the chart evaluator.
    """
    delta=a-b;lo=np.nextafter(delta,-np.inf);hi=np.nextafter(delta,np.inf)
    upper=np.maximum(abs(lo),abs(hi))
    lower=np.where((lo<=0)&(hi>=0),0.,np.minimum(abs(lo),abs(hi)))
    low2=np.maximum(0.,np.nextafter(lower*lower,-np.inf))
    high2=np.nextafter(upper*upper,np.inf)
    low_sum=np.zeros(low2.shape[:-1]);high_sum=np.zeros(high2.shape[:-1])
    for axis in range(3):
        low_sum=np.maximum(0.,np.nextafter(low_sum+low2[...,axis],-np.inf))
        high_sum=np.nextafter(high_sum+high2[...,axis],np.inf)
    return (np.maximum(0.,np.nextafter(np.sqrt(low_sum),-np.inf)),
            np.nextafter(np.sqrt(high_sum),np.inf))


def _sum_upper(a,b):
    return np.nextafter(a+b,np.inf)


def _gap_lower(distance_lower,radius):
    return np.nextafter(np.nextafter(distance_lower-radius[:,None],-np.inf)-radius[None,:],-np.inf)


def _summary(values,mask):
    selected=values[mask]
    return dict(minimum=float(selected.min()),maximum=float(selected.max())) if len(selected) else None


def _verified_report(path,states,geometry):
    data=json.loads(path.read_text())
    if not data['complete']:raise ValueError('complete archived cost/control report required')
    for actual,expected in ((states,data['states_sha256']),
            (geometry.with_suffix('.npz'),data['geometry_sha256']),
            (geometry,data['geometry_json_sha256'])):
        if sha256(actual.read_bytes()).hexdigest()!=expected:
            raise ValueError('control input hash mismatch: '+str(actual))
    archive=path.with_name(path.stem+'-source')
    for name,digest in data['source_sha256'].items():
        if sha256((archive/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('archived producer hash mismatch: '+name)
    return data


def run(states,geometry,reverse_report,factor_report,output):
    states,geometry,reverse_report,factor_report,output=map(Path,
        (states,geometry,reverse_report,factor_report,output))
    if output.exists():raise FileExistsError(output)
    folder=Path(__file__).parent
    names=('boat_hybrid_mask_probe.py','boat_geometry_template.py','boat_world_field.py')
    report=dict(complete=False,qualified=False,scope=__doc__,
        source_sha256={name:sha256((folder/name).read_bytes()).hexdigest() for name in names},
        policy=dict(opening_factor=OPENING_FACTOR,slender_quality=SLENDER_QUALITY,
            contact='All rows/columns incident to a body/free shared node exact; conservative heuristic, not mathematically necessary.',
            slender='All rows/columns of a face with endpoint quality below threshold exact; heuristic, no interval aspect-ratio proof.',
            touching='Any shared node ID exact, including self.',
            margin='Maximum held-out vertex displacement from anchor; retrospective geometric budget, not online lifetime evidence.'))
    try:
        reverse=_verified_report(reverse_report,states,geometry)
        factor=_verified_report(factor_report,states,geometry)
        # Retain actual geometry owner, not a port of its clipping/refinement.
        for name in names[1:]:
            if report['source_sha256'][name]!=factor['source_sha256'][name]:
                raise ValueError('current geometry producer differs from archived control: '+name)
        metadata=json.loads(geometry.read_text())
        if not metadata['complete'] or not metadata['geometry_roundoff_checks_passed']:
            raise ValueError('accepted complete original boundary required')
        reference=dict(np.load(geometry.with_suffix('.npz'),allow_pickle=False))
        saved=np.load(states,allow_pickle=False)
        times=saved['accepted_times'][:3];sequence=saved['accepted_states'][:3]
        if (times.shape!=(3,) or not np.allclose(times,[0.,.01,.02],rtol=0,atol=1e-14)
                or sequence.shape!=(3,43) or not np.isfinite(sequence).all()):
            raise ValueError('this bounded screen requires the actual 0/10/20 ms 43-state fixture')
        frames=[];guards=[];n=len(reference['parent_nodes'])
        for state in sequence:
            arrays,guard=evaluate_template(reference,pose=InitialState(state[:3],state[3:7]/np.linalg.norm(state[3:7])),
                parent_eta=state[13:13+n],parent_psi=state[13+n:])
            frames.append(refine_boundary(arrays,region='body'));guards.append(guard)
        anchor=frames[0];faces=anchor['faces'];nf=len(faces)
        if nf!=472:raise ValueError('this screen owns only the archived 472-face case')
        basis=('faces','body_faces','free_faces','free_nodes','wall_faces','bottom_faces',
               'body_source_group','free_source_parent','free_parent_indices','node_chart_key')
        for arrays in frames[1:]:
            if (arrays['points'].shape!=anchor['points'].shape
                    or any(not np.array_equal(arrays[key],anchor[key]) for key in basis)):
                raise ValueError('held-out frame left the exact semantic chart/basis')
        triangles=np.array([arrays['points'][faces] for arrays in frames])
        if not np.isfinite(triangles).all():raise ValueError('nonfinite owned geometry')
        centers=triangles[0].mean(axis=1)
        radius=_norm_bounds(triangles[0],centers[:,None,:])[1].max(axis=1)
        moves=_norm_bounds(triangles,triangles[:1])[1].max(axis=2)
        motion=moves.max(axis=0)
        swept_radius=_sum_upper(radius,motion)
        distance=_norm_bounds(centers[:,None,:],centers[None,:,:])[0]
        gap=_gap_lower(distance,swept_radius)
        size_sum=_sum_upper(swept_radius[:,None],swept_radius[None,:])
        threshold=np.nextafter(OPENING_FACTOR*size_sum,np.inf)
        geometric_far=(distance>threshold)&(gap>0)
        touching=(faces[:,None,:,None]==faces[None,:,None,:]).any(axis=(2,3))
        contact_nodes=np.intersect1d(faces[anchor['body_faces']].ravel(),faces[anchor['free_faces']].ravel())
        contact_face=np.isin(faces,contact_nodes).any(axis=1)
        edges=np.roll(triangles,-1,axis=2)-triangles
        area=np.linalg.norm(np.cross(triangles[:,:,1]-triangles[:,:,0],triangles[:,:,2]-triangles[:,:,0]),axis=2)/2
        longest2=np.max(np.sum(edges*edges,axis=3),axis=2)
        if np.any(area<=0) or np.any(longest2<=0):raise ValueError('degenerate held-out face')
        quality=area/longest2
        slender=quality.min(axis=0)<SLENDER_QUALITY
        contact_pairs=contact_face[:,None]|contact_face[None,:]
        slender_pairs=slender[:,None]|slender[None,:]
        protected=touching|contact_pairs|slender_pairs
        far=geometric_far&~protected
        if not np.array_equal(far,far.T) or np.diag(far).any():raise ValueError('mask symmetry/self failure')
        endpoints=[];intervals=[]
        for time,tri in zip(times,triangles,strict=True):
            # Keep the chosen anchor centers: recentered enclosing spheres
            # could overlap even when the original certified spheres separate.
            rad=_norm_bounds(tri,centers[:,None,:])[1].max(axis=1)
            current_gap=_gap_lower(distance,rad)
            if far.any() and np.min(current_gap[far])<=0:raise ValueError('held-out endpoint separation failed')
            endpoints.append(dict(time=float(time),far_ball_gap_lower_m=_summary(current_gap,far)))
        for k in range(2):
            segment_radius=_sum_upper(radius,np.maximum(moves[k],moves[k+1]))
            segment_gap=_gap_lower(distance,segment_radius)
            if far.any() and np.min(segment_gap[far])<=0:raise ValueError('nodal-chord separation bound failed')
            intervals.append(dict(times=times[k:k+2].tolist(),far_ball_gap_lower_m=_summary(segment_gap,far)))
        count=int(far.sum());total=nf*nf;near=total-count;fraction=count/total
        # Cost figures are archived measurements, NOT timings taken by this script.
        reverse_rows=[row for row in reverse['rows'] if row['label']=='reverse']
        if reverse['order']!=3 or reverse['faces']!=472 or not reverse_rows:
            raise ValueError('matching q3 reverse cost report required')
        gradient=float(reverse_rows[-1]['seconds'])
        mixed_factor=float(factor['samples'][0]['candidate_diagnostics']['mixed_factor_and_solve_seconds'])
        gradient_samples=[float(row['seconds']) for row in reverse_rows]
        mixed_samples=[float(row['candidate_diagnostics']['mixed_factor_and_solve_seconds']) for row in factor['samples']]
        report.update(complete=True,nodes=len(anchor['points']),faces=nf,times=times.tolist(),
            input_sha256={str(path):sha256(path.read_bytes()).hexdigest() for path in
                (states,geometry,geometry.with_suffix('.npz'),reverse_report,factor_report)},
            node_chart_key=str(anchor['node_chart_key']),guard_diagnostics=guards,
            separation=dict(scope='Supplied floating endpoint triangles and all piecewise-linear nodal chords joining them; NOT nonlinear ALE paths.',
                arithmetic='Individually outward-widened float64 basic operations on supplied coordinates; not an enclosure of exact CSG geometry.',
                maximum_vertex_motion_m=float(motion.max()),global_far_gap_lower_m=_summary(gap,far),
                endpoints=endpoints,intervals=intervals,actual_nonlinear_inter_sample_path_certified=False),
            counts=dict(directed_total=total,directed_far=count,directed_exact_near=near,far_fraction=fraction,
                unordered_far=count//2,geometric_far_before_heuristics=int(geometric_far.sum()),
                touching_directed=int(touching.sum()),contact_nodes=len(contact_nodes),contact_faces=int(contact_face.sum()),
                slender_faces=int(slender.sum()),minimum_endpoint_quality=float(quality.min()),
                geometric_far_removed_by_heuristics=int((geometric_far&protected).sum())),
            mask_sha256=sha256(np.packbits(far.ravel()).tobytes()).hexdigest(),
            retained_payload_estimate=dict(label='Arithmetic estimate only; no allocation/cache/peak-memory measurement.',
                packed_far_six_outputs_bytes=count*6*19*8,packed_far_seven_outputs_bytes=count*7*19*8,
                packed_far_pair_index_bytes=count*8,fixed_boolean_mask_bytes=total,
                full_1088_six_outputs_bytes=1088**2*6*19*8,full_1088_seven_outputs_bytes=1088**2*7*19*8,
                excludes='Factors/E/H/current geometry/temporary work/CuPy pool/local kernel memory; 256MiB is not a user cap.'),
            cost_screen=dict(label='Optimistic uniform-pair count proxy, not a measured runtime or load-error estimate.',
                observed_rhs_count=OBSERVED_RHS,observed_cold_build_count=OBSERVED_COLD_BUILDS,
                rhs_per_cold_build_screening_ratio=USES_PER_COLD_BUILD,anchor_lifetime_measured=False,
                archived_whole_reverse_gradient_seconds=gradient,archived_mixed_factor_and_solve_seconds=mixed_factor,
                archived_timing_selection='Last reverse-labeled row and first factor sample; screening arithmetic uses these explicitly selected archived observations.',
                archived_reverse_seconds_samples=gradient_samples,archived_mixed_factor_seconds_samples=mixed_samples,
                V_factor_cost_separately_measured=False,factorization_eliminated_by_mask=False,
                maximum_pair_speedup_if_far_free=total/near,
                optimistic_gradient_seconds_if_far_free=gradient*(1-fraction),
                optimistic_gradient_savings_seconds_per_rhs=gradient*fraction,
                optimistic_anchor_build_break_even_seconds=USES_PER_COLD_BUILD*gradient*fraction,
                caveats='64/29 includes initial cold build and reuse of older charts, not actual anchor lifetime. Whole-gradient packing/solves are not all pair-scaled; near/far pair costs vary, Taylor hits cost bandwidth, and masks do not remove current matrix factorization. This one conservative policy cannot reject every less-conservative hybrid.'))
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('states','geometry','reverse_report','factor_report','output'):
        parser.add_argument('--'+name.replace('_','-'),required=True)
    print(json.dumps(run(**vars(parser.parse_args())),indent=2))

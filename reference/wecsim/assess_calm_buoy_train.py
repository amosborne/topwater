"""Read-only reassessment of immutable train2 under explicit REF002 calm policy."""
import argparse
from collections import deque
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits
from calm_buoy_train import CalmBuoyTrain


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


class CombinedBoxBounds:
    def __init__(self,train,order,deadline):
        C = train.C.reshape(-1)
        theta = np.broadcast_to(train.theta[:,None,None],train.C.shape).reshape(-1)
        sign = np.broadcast_to(np.array([1.,-1.])[None,:,None],train.C.shape).reshape(-1)
        k = np.broadcast_to(train.k,train.C.shape).reshape(-1)
        omega = np.broadcast_to(train.omega,train.C.shape).reshape(-1)
        active = abs(C) > 0
        self.C = C[active]
        ka = k[active]
        assert np.all(ka > 0),'DC coefficient must remain naturally zero'
        kx = sign[active]*ka*np.cos(theta[active])
        ky = sign[active]*ka*np.sin(theta[active])
        self.k = np.c_[kx,ky,ka]
        self.V = np.c_[omega[active]*kx/ka,omega[active]*ky/ka,-1j*omega[active]]
        self.powers = [np.array([v**i for i in range(order+1)]) for v in [1j*kx,1j*ky,ka]]
        self.weight_eta = abs(self.C)
        self.weight_u = np.sqrt(2)*omega[active]*abs(self.C)
        self.order = order
        self.deadline = deadline

    def evaluate(self,boxes,dimension,batch):
        assert dimension in [2,3]
        bounds = []
        centers_value = []
        for start in range(0,len(boxes),batch):
            if time.monotonic() > self.deadline:
                raise TimeoutError('Declared combined-field assessment budget exceeded')
            b = np.asarray(boxes[start:start+batch])
            center = (b[:,0]+b[:,1])/2
            h = (b[:,1]-b[:,0])/2
            if dimension == 2:
                phase = self.C*np.exp(1j*(center@self.k[:,:2].T))
                weights = self.weight_eta
            else:
                assert np.all(b[:,1,2] <= 1e-12)
                phase = self.C*np.exp(1j*(center[:,:2]@self.k[:,:2].T)+center[:,2,None]*self.k[:,2])
                weights = self.weight_u
            total = np.zeros(len(b))
            center_magnitude = None
            for alpha in itertools.product(range(self.order+1),repeat=dimension):
                if sum(alpha) > self.order:
                    continue
                if time.monotonic() > self.deadline:
                    raise TimeoutError('Declared combined-field assessment budget exceeded')
                factor = np.ones(len(self.C),dtype=complex)
                halfwidth_factor = np.ones(len(b))
                for axis,power in enumerate(alpha):
                    factor *= self.powers[axis][power]
                    halfwidth_factor *= h[:,axis]**power/math.factorial(power)
                if dimension == 2:
                    value = abs((phase@factor).real)
                else:
                    value = np.linalg.norm((phase@(factor[:,None]*self.V)).real,axis=1)
                if sum(alpha) == 0:
                    center_magnitude = value.copy()
                total += value*halfwidth_factor
            r = h@abs(self.k[:,:dimension]).T
            remainder = r**(self.order+1)@weights/math.factorial(self.order+1)
            bounds.extend((total+remainder).tolist())
            centers_value.extend(center_magnitude.tolist())
        return np.asarray(bounds),np.asarray(centers_value)

    def certify(self,boxes,dimension,batch,limit,maxlevel):
        pending = deque((box,0) for box in boxes)
        leaves = []; upper = []; centers = []; levels = []
        budget_exceeded = False
        violation = False
        while pending:
            group = [pending.popleft() for _ in range(min(batch,len(pending)))]
            try:
                bound,value = self.evaluate([b for b,_ in group],dimension,batch)
            except TimeoutError:
                pending.extendleft(reversed(group))
                budget_exceeded = True
                break
            for (box,level),b,v in zip(group,bound,value):
                if b > limit and v <= limit and level < maxlevel:
                    middle = box.mean(axis=0)
                    for side in itertools.product([0,1],repeat=dimension):
                        child = np.array([[box[0,d] if s == 0 else middle[d] for d,s in enumerate(side)],
                            [middle[d] if s == 0 else box[1,d] for d,s in enumerate(side)]])
                        pending.append((child,level+1))
                else:
                    leaves.append(box); upper.append(float(b)); centers.append(float(v)); levels.append(level)
                    violation = violation or bool(v > limit)
            if violation:
                break
        complete = not pending
        proved_pass = complete and bool(upper) and max(upper) <= limit
        return dict(proved_pass=proved_pass,actual_center_violation=violation,
            coverage_complete=complete,budget_exceeded=budget_exceeded,
            unresolved_leaf_boxes=sum(b > limit for b,v in zip(upper,centers) if v <= limit),
            unvisited_boxes=len(pending),maximum_upper_bound=max(upper,default=None),
            maximum_center_value=max(centers,default=None),leaf_boxes=np.asarray(leaves),
            leaf_upper_bounds=np.asarray(upper),leaf_center_values=np.asarray(centers),levels=np.asarray(levels),
            pending_boxes=np.asarray([b for b,_ in pending]))


def boxes_for(lo,hi,width):
    axes = [np.linspace(l,r,int(np.ceil((r-l)/w))+1) for l,r,w in zip(lo,hi,width)]
    return [np.array([[axes[d][i] for d,i in enumerate(index)],
                       [axes[d][i+1] for d,i in enumerate(index)]])
            for index in itertools.product(*(range(len(a)-1) for a in axes))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    started = time.monotonic()
    control_path = Path('reference/wecsim/calm-buoy-assessment-controls.json')
    controls = json.loads(control_path.read_text())
    source = Path(controls['immutable_train_directory'])
    old = json.loads((source/'diagnostic.json').read_text())
    before = {p.as_posix():sha(p) for p in source.iterdir() if p.is_file()}
    train = CalmBuoyTrain(source/'train.npz')
    forcing = loadmat(source/'body-force.mat',simplify_cells=True)
    F0 = forcing['force'][0]
    assert np.array_equal(F0,np.asarray(old['initial']['initial_force_N_Nm']))
    S = np.asarray(old['initial']['carrier_force_scale'])
    physical_path = Path('reference/capytaine/canonical-boat.json')
    physical = json.loads(physical_path.read_text())
    eps = controls['epsilon']; assert eps == .001
    mass = physical['body_mass_kg']; gravity = physical['water']['gravity_m_s2']
    hull = physical['minimum_hull']
    material = np.asarray(hull['vertices_material_m'])
    Lchar = np.ptp(material,axis=0)[:2].max()
    force_limit = eps*mass*gravity
    moment_limit = force_limit*Lchar
    force_norm = np.linalg.norm(F0[:3]); moment_norm = np.linalg.norm(F0[3:])
    force_policy = dict(epsilon=eps,mass_kg=mass,gravity_m_s2=gravity,Lchar_m=float(Lchar),
        native_force_N_Nm=F0.tolist(),carrier_S_same_native_basis=S.tolist(),
        peraxis_relative_limits=(eps*S).tolist(),
        norm_force_N=float(force_norm),norm_force_limit_N=force_limit,
        norm_moment_Nm=float(moment_norm),norm_moment_limit_Nm=float(moment_limit),
        force_over_weight=float(force_norm/(mass*gravity)),moment_over_weight_length=float(moment_norm/(mass*gravity*Lchar)),
        policy_status='coordinator-declared engineering numericalcalm assumption, not physicalaccuracy guarantee',
        old_REF001_absolute_cap_N_Nm=.1,old_policy_rejected=True)
    gates = dict(force_norm=bool(force_norm <= force_limit),moment_norm=bool(moment_norm <= moment_limit),
        peraxis_carrier_relative=bool(np.all(abs(F0) <= eps*S)),arrival=old['gates']['arrival'])
    args.output.mkdir(parents=True,exist_ok=False)
    write(args.output/'force-policy.json',force_policy)
    bounds = None
    if all(gates.values()):
        q = hull['equilibrium_quaternion_wxyz']
        R = Rotation.from_quat([*q[1:],q[0]]).as_matrix()
        v = (material-hull['cg_material_m'])@R.T+physical['world_cg_m']
        wet = np.asarray(physical['minimum_wetted_triangles_world_m']).reshape(-1,3)
        lo = np.r_[v[:,:2].min(axis=0),wet[:,2].min()]
        hi = np.r_[v[:,:2].max(axis=0),0.]
        widths = np.asarray(controls['maximum_box_width_xyz_m'])
        xyboxes = boxes_for(lo[:2],hi[:2],widths[:2])
        xyzboxes = boxes_for(lo,hi,widths)
        proof = CombinedBoxBounds(train,controls['Taylor_order'],started+controls['maximum_elapsed_seconds'])
        eta_proof = proof.certify(xyboxes,2,controls['center_batch_size'],old['initial_eta_limit_m'],controls['maximum_refinement_level'])
        speed_proof = proof.certify(xyzboxes,3,controls['center_batch_size'],old['initial_speed_limit_m_s'],controls['maximum_refinement_level'])
        raw = {}
        for name,result in [('eta',eta_proof),('speed',speed_proof)]:
            for key in ['leaf_boxes','leaf_upper_bounds','leaf_center_values','levels','pending_boxes']:
                raw[name+'_'+key] = result.pop(key)
        np.savez(args.output/'box-bounds.npz',**raw)
        bounds = dict(eta=eta_proof,speed=speed_proof,region_min_xyz_m=lo.tolist(),region_max_xyz_m=hi.tolist(),
            initial_eta_boxes=len(xyboxes),initial_velocity_boxes=len(xyzboxes),Taylor_order=controls['Taylor_order'],
            method='combined ALLmode derivatives before abs/norm; exact exponential remainder; adaptive COMPUTATIONALboxes only')
        gates['eta'] = eta_proof['proved_pass']
        gates['velocity'] = speed_proof['proved_pass']
    assert all(sha(Path(p)) == digest for p,digest in before.items())
    actual_violation = bool(bounds and (bounds['eta']['actual_center_violation'] or bounds['speed']['actual_center_violation']))
    field_status = 'proved_calm' if all(gates.values()) else ('actual_center_violation' if actual_violation else 'bound_unresolved_or_budget_limited')
    report = dict(complete=True,accepted=bool(all(gates.values())),gates=gates,field_status=field_status,
        force_policy=force_policy,bounds=bounds,immutable_train=str(source),immutable_input_sha256=before,
        old_rejection_preserved=True,initial_eta_limit_m=old['initial_eta_limit_m'],
        initial_speed_limit_m_s=old['initial_speed_limit_m_s'],arrival=old['arrival'],
        elapsed_s=time.monotonic()-started,body_simulated=False,
        source_sha256=sha(__file__),controls_sha256=sha(control_path),
        scope='New explicitly declared assessment ONLY of same train2, no physicalinput change/historyseed/retiming')
    write(args.output/'diagnostic.json',report)
    print(json.dumps(report),flush=True)
    return 0 if report['accepted'] else 2


if __name__ == '__main__':
    with threadpool_limits(1):
        raise SystemExit(main())

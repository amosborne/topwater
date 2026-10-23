"""Owned CPU replay and body-only rebuild for regularized-water diagnostics."""
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.fft import rfft2, irfft2
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import splu
from topwater import SimulationInputs
from topwater.rankine import TrianglePanels, panel_influence
from tests.laplace_reference import triangle_samples
from .physical_p1_patch import PhysicalP1Patch
from .compact_p1_maps import body_maps
from .boat_spatial_work import source_matrix
from .panel_gradient_reference import checked_influence


def owned(path):
    path = Path(path)
    report = json.loads(path.read_text())
    if not report['complete'] or sha256(path.with_suffix('.npz').read_bytes()).hexdigest() != report['arrays_sha256']:
        raise ValueError('complete SHA-owned snapshot required: '+str(path))
    return report


class SavedPhysicalFlux:
    """Exact recorded kernel/contact-row action replayed with CPU FFTs."""
    def __init__(self, patch, arrays):
        self.patch = patch
        self.shape = tuple(int(v) for v in arrays['weak_full_shape'])
        self.fft_shape = tuple(int(v) for v in arrays['weak_fft_shape'])
        kernel = arrays['weak_full_kernel']
        self.rows = arrays['physical_contact_active_rows'].copy()
        self.contact = arrays['physical_contact_rows'].copy()
        if self.shape != (len(patch.y)-2, len(patch.x)-2) or self.contact.shape != (len(self.rows), len(patch.active_nodes)):
            raise ValueError('recorded operator dimensions differ')
        position = np.full(len(patch.nodes), -1, dtype=int)
        position[patch.interior_nodes] = np.arange(len(patch.interior_nodes))
        self.active = position[patch.active_nodes]
        embed = np.zeros(self.fft_shape)
        dy, dx = np.indices(kernel.shape)
        embed[(dy-self.shape[0]+1)%self.fft_shape[0], (dx-self.shape[1]+1)%self.fft_shape[1]] = kernel
        self.spectrum = rfft2(embed)
        replay = self(arrays['operator_audit_coefficients'])
        self.audit_error = float(abs(replay-arrays['operator_audit_loads']).max())
        np.testing.assert_allclose(replay, arrays['operator_audit_loads'], rtol=1e-11, atol=1e-12)

    def __call__(self, values):
        if np.iscomplexobj(values):
            return self(values.real)+1j*self(values.imag)
        vector = values.ndim == 1
        if vector:
            values = values[:, None]
        expanded = self.patch.expand(values)[self.patch.interior_nodes].reshape(*self.shape, -1)
        padded = np.zeros((*self.fft_shape, values.shape[1]))
        padded[:self.shape[0], :self.shape[1]] = expanded
        result = irfft2(rfft2(padded, axes=(0, 1))*self.spectrum[:, :, None],
            s=self.fft_shape, axes=(0, 1))[:self.shape[0], :self.shape[1]].reshape(-1, values.shape[1])[self.active].copy()
        result[self.rows] = self.contact @ values
        return result[:, 0] if vector else result


class RegularizedSnapshot:
    def __init__(self, water, phase, mode):
        self.water, self.phase, self.mode = owned(water), owned(phase), owned(mode)
        hashes = {name: sha256(Path(path).read_bytes()).hexdigest() for name, path in
                  [('water', water), ('phase', phase)]}
        if self.mode['parents'] != hashes or self.water['parent_phase_sha256'] != hashes['phase']:
            raise ValueError('mode/water/phase ownership differs')
        if any(r['input_sha256'] != self.mode['input_sha256'] for r in (self.phase, self.water)):
            raise ValueError('snapshot physical inputs differ')
        for name in ('src/topwater/rankine.py', 'src/topwater/_rankine_kernel.py',
                     'src/topwater/_rankine_extended.py', 'src/topwater/numerics.py',
                     'src/topwater/constants.py',
                     'tests/experiments/free_surface/panel_gradient_reference.py',
                     'tests/experiments/free_surface/boat_spatial_work.py',
                     'tests/experiments/free_surface/compact_p1_maps.py'):
            path = Path(name)
            expected = self.phase['producer_sha256'].get(str(path))
            if sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('recorded phase kernel source changed: '+name)
        for name in ('physical_p1_patch.py', 'compact_p1_trace.py', 'compact_p1_maps.py'):
            relative = Path('tests/experiments/free_surface')/name
            if sha256(relative.read_bytes()).hexdigest() != self.mode['producer_sha256'].get(str(relative)):
                raise ValueError('recorded field/extension source changed: '+name)
        contact_path = Path(phase).with_suffix('.source')/'contact.json'
        if sha256(contact_path.read_bytes()).hexdigest() != self.phase['parent_contact_sha256']:
            raise ValueError('contact ownership differs')
        contact = json.loads(contact_path.read_text())
        audit_path = Path(contact['parents']['audit']['path'])
        # Recorded paths were relative to the archived original repository.
        if not audit_path.is_absolute():
            audit_path = Path(water).resolve().parent.parent/audit_path
        audit = owned(audit_path)
        if sha256(audit_path.read_bytes()).hexdigest() != contact['parents']['audit']['sha256']:
            raise ValueError('body/input audit changed')
        self.inputs = SimulationInputs.from_json((audit_path.with_suffix('.source')/'inputs.json').read_text())
        if sha256(self.inputs.to_json().encode()).hexdigest() != self.mode['input_sha256']:
            raise ValueError('physical input ownership differs')
        with np.load(audit_path.with_suffix('.npz')) as z:
            self.triangles = z['triangles'].copy()
        if sha256(self.triangles.tobytes()).hexdigest() != self.mode['panel_sha256']:
            raise ValueError('body panel ownership differs')
        with np.load(Path(water).with_suffix('.npz')) as z:
            data = {name: z[name].copy() for name in z.files}
        with np.load(Path(phase).with_suffix('.npz')) as z:
            self.patch = PhysicalP1Patch(data['x'], data['y'], z['polygon'])
            np.testing.assert_array_equal(self.patch.active_nodes, z['active_nodes'])
            for current, key in zip((self.patch.mass.data, self.patch.mass.indices, self.patch.mass.indptr),
                                    ('mass_data', 'mass_indices', 'mass_indptr'), strict=True):
                np.testing.assert_array_equal(current, z[key])
            order = self.phase['orders'][-1]['order']
            self.coarse = dict(n=z['nmap'].copy(), j=z['neumann'].copy(), cross=z[f'q{order}_cross'].copy())
        self.flux = SavedPhysicalFlux(self.patch, data)
        np.testing.assert_array_equal(self.patch.active_nodes, data['active_nodes'])
        self.mass_factor = splu(self.patch.mass.tocsc())
        with np.load(Path(mode).with_suffix('.npz')) as z:
            self.psi = z['mode_0_psi'].copy()
            self.eta = z['mode_0_eta'].copy()
            self.growth = complex(z['eigenvalues'][0])

    def mass_solve(self, values):
        if np.iscomplexobj(values):
            return self.mass_factor.solve(values.real)+1j*self.mass_factor.solve(values.imag)
        return self.mass_factor.solve(values)

    def body_system(self, triangles, order, *, maximum_pairs=500_000_000, maximum_bytes=512*1024**2):
        patch = self.patch
        panels = TrianglePanels(triangles, -1)
        nb = len(panels)
        points = len(patch.cut_triangles)*order*order
        pairs = nb*len(patch.panels)+4*nb*nb+nb*points
        memory = 8*(4*nb*len(patch.nodes)+4*nb*nb+128*nb+1024*len(patch.panels))
        if pairs > maximum_pairs or memory > maximum_bytes:
            raise RuntimeError('body-only construction preflight exceeds work/storage budget')
        p, n = body_maps(patch, panels.centroids, panels.normals)
        n = patch.condense_maps(n)
        del p
        j = source_matrix(panels)
        for first in range(0, nb, 128):
            last = min(first+128, nb)
            image = checked_influence(panels.centroids[first:last]*[1., 1., -1.], panels)
            j[first:last] += 2*np.einsum('mni,mi->mn', image.gradient_single_layer,
                                       panels.normals[first:last]*[1., 1., -1.])
        bary, weights = triangle_samples(np.eye(3)[None], order)
        bary = bary[0]
        cross = np.zeros((len(patch.nodes), nb))
        for first in range(0, len(patch.cut_triangles), 16):
            last = min(first+16, len(patch.cut_triangles))
            xy = np.einsum('qj,tji->tqi', bary, patch.cut_triangles[first:last])
            basis = np.einsum('qj,tji->tqi', bary, patch.basis_values[first:last])
            w = patch.cut_areas[first:last, None, None]*weights[None, :, None]*basis
            targets = np.column_stack((xy.reshape(-1, 2), np.zeros(len(xy)*len(weights))))
            nodes = np.repeat(patch.faces[patch.parent_id[first:last]], len(weights), axis=0)
            w = w.reshape(-1, 3)
            for start in range(0, len(targets), 128):
                stop = min(start+128, len(targets))
                source = 2*panel_influence(targets[start:stop], panels).single_layer
                for local in range(3):
                    np.add.at(cross, nodes[start:stop, local], w[start:stop, local, None]*source)
        return dict(n=n, j=j, cross=cross[patch.active_nodes], source_pairs_upper_bound=pairs,
                    constructor_array_bytes_bound=memory)

    def solve(self, system, psi):
        inverse_cross = self.mass_factor.solve(system['cross'])
        schur = system['j']-system['n']@inverse_cross
        factor = lu_factor(schur)
        sigma = -lu_solve(factor, system['n']@psi)
        chi = psi-inverse_cross@sigma
        return chi, sigma, inverse_cross, factor

    def body_field(self, triangles, chi, sigma, points, normals):
        full = self.patch.expand(chi)
        p, n = body_maps(self.patch, points, normals)
        value, normal = p@full, n@full
        del p, n
        panels = TrianglePanels(triangles, -1)
        for first in range(0, len(points), 128):
            last = min(first+128, len(points))
            real = checked_influence(points[first:last], panels)
            image = checked_influence(points[first:last]*[1., 1., -1.], panels)
            value[first:last] += (real.single_layer+image.single_layer)@sigma
            gradient = real.gradient_single_layer+image.gradient_single_layer*[1., 1., -1.]
            normal[first:last] += np.einsum('mni,mi->mn', gradient, normals[first:last])@sigma
        return value, normal

"""Independent Rankine comparison, not wave-force/hydrodynamic correlation.

Run: uv run --locked --group correlation python -m tests.test_capytaine_rankine
Capytaine's far/small-distance approximations are reported separately from its
analytic branch; they are not a reason to relax Topwater's quadrature gates.
"""

import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest

cpt = pytest.importorskip("capytaine")

from topwater import HullMesh
from topwater.bem import assemble_laplace
from topwater.rankine import TrianglePanels
from .rankine_reference import panel_reference
from .analytic_geometry import sphere_mesh


def compare_rankine(refinement):
    assert cpt.__version__ == '3.0.0', 'review the reference branch rules before upgrading'
    vertices, faces = sphere_mesh(refinement)
    panels = TrianglePanels.from_mesh(HullMesh(vertices, faces))
    # Preserve the exact triangles/order. Do not let mesh healing change the case.
    mesh = cpt.Mesh(vertices=vertices, faces=faces, auto_clean=False, auto_check=False)
    assert_allclose(mesh.faces_centers, panels.centroids, atol=0, rtol=0)
    assert_allclose(mesh.faces_areas, panels.areas, atol=1e-14, rtol=1e-12)
    assert_allclose(mesh.faces_normals, panels.normals, atol=1e-14, rtol=1e-12)
    assert np.all(panels.edge_lengths/mesh.faces_radiuses[:, None] >= float(np.float32(1e-3)))
    green = cpt.Delhommeau(tabulation_nr=0, tabulation_nz=0, tabulation_cache_dir=None)
    reference_s, reference_k = green.evaluate_rankine_only(
        mesh, mesh, adjoint_double_layer=False, diagonal_term_in_double_layer=False)
    operators = assemble_laplace(panels)
    # Capytaine uses negative G. Its Fortran PI initializer is rounded to default
    # real precision before conversion to float64 in 3.0.0. Record the actual
    # constant and explicitly normalize to the same mathematical 1/(4*pi*r).
    reference_pi = float(green.fortran_core.constants.pi)
    conversion = -reference_pi/np.pi
    converted_s, converted_k = reference_s*conversion, reference_k*conversion
    total = stable_count = 0
    coverage = np.zeros(len(panels), dtype=int)
    stable_s_error = stable_k_error = 0.
    for first in range(0, len(panels), 128):
        last = min(first+128, len(panels))
        delta = panels.centroids[first:last, None, None, :]-panels.triangles[None]
        height = np.einsum('mnj,nj->mn', delta[:, :, 0], panels.geometric_normals)
        edge_distance = np.einsum('mnej,nej->mne', delta, panels.edge_outward)
        separation = np.linalg.norm(panels.centroids[first:last, None]-panels.centroids[None], axis=2)
        radius = mesh.faces_radiuses[None]
        # These are the reference implementation's explicit branches, not fitted
        # masks based on observed agreement. See the tagged Fortran source.
        stable = ((separation <= 7*radius)
                  & ((np.abs(height) >= float(np.float32(1e-4))*radius)
                     | (np.abs(height) <= 64*np.finfo(float).eps*panels.diameters))
                  & np.all(np.abs(edge_distance) >= float(np.float32(1e-5)), axis=2))
        diameter = np.broadcast_to(panels.diameters, stable.shape)
        ours_s, ours_k = operators.single_layer[first:last], operators.double_layer[first:last]
        theirs_s, theirs_k = converted_s[first:last], converted_k[first:last]
        assert_allclose(theirs_s[stable]/diameter[stable], ours_s[stable]/diameter[stable],
                        atol=1e-10, rtol=1e-8)
        assert_allclose(theirs_k[stable], ours_k[stable], atol=1e-10, rtol=1e-8)
        stable_s_error = max(stable_s_error, float(np.max(np.abs(theirs_s-ours_s)[stable])))
        stable_k_error = max(stable_k_error, float(np.max(np.abs(theirs_k-ours_k)[stable])))
        stable_count += int(stable.sum())
        total += stable.size
        stable[np.arange(last-first), np.arange(first, last)] = False
        coverage += stable.sum(axis=0)
    assert np.all(coverage > 0), 'each source panel needs independent nonself coverage'
    # Independently check the worst disagreements too, including the reference's
    # approximated/cutoff branches. Do not turn them into a second loose gate.
    worst = []
    for values, ours, component in ((converted_s, operators.single_layer, 0),
                                     (converted_k, operators.double_layer, 1)):
        scale = panels.diameters[None, :] if component == 0 else 1.
        row, col = np.unravel_index(np.argmax(np.abs(values-ours)/scale), ours.shape)
        integral = panel_reference(panels.centroids[row], panels.triangles[col])[component]
        scale = panels.diameters[col] if component == 0 else 1.
        assert_allclose(ours[row, col]/scale, integral/scale, atol=1e-10, rtol=1e-8)
        worst.append(dict(component='S' if component == 0 else 'K', target=int(row),
                          source=int(col), topwater=float(ours[row, col]),
                          capytaine_converted=float(values[row, col]), independent=float(integral)))
    return dict(capytaine=cpt.__version__, panel_count=len(panels), pair_count=total,
                stable_pair_count=stable_count, reference_pi=reference_pi,
                stable_nonself_pair_count=int(coverage.sum()),
                minimum_nonself_comparisons_per_source=int(coverage.min()),
                stable_max_s_error=stable_s_error, stable_max_k_error=stable_k_error,
                raw_max_s_sign_converted_error=float(np.max(np.abs(reference_s+operators.single_layer))),
                raw_max_k_sign_converted_error=float(np.max(np.abs(reference_k+operators.double_layer))),
                worst_disagreements=worst)


@pytest.mark.parametrize('refinement', [2, 3, 4])
def test_capytaine_rankine_analytic_branch_and_independent_worst_pairs(refinement):
    report = compare_rankine(refinement)
    path = Path(f'artifacts/capytaine-rankine-{refinement}.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+"\n", encoding='utf-8')


if __name__ == '__main__':
    reports = [compare_rankine(refinement) for refinement in (2, 3, 4)]
    path = Path('artifacts/capytaine-rankine.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(reports, indent=2)+"\n", encoding='utf-8')
    print(json.dumps(reports, indent=2))

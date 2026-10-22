"""Render the simple boat's actual mesh and solved calm-water pose.

Run: uv run --locked --with matplotlib python -m tests.plot_reference_boat
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
import numpy as np

from topwater.frames import rotation_matrix
from .reference_boat import (
    BEAM, LENGTH, MOTOR_MASS, MOTOR_POSITION, PAYLOAD_MASS, PAYLOAD_POSITION,
    SHELL_MASS, equilibrium, reference_body,
)


def main():
    body = reference_body()
    rest, result = equilibrium()
    rotation = rotation_matrix(rest.orientation)
    vertices = (body.mesh.vertices-body.cg)@rotation.T+rest.position
    triangles = vertices[body.mesh.faces[body.mesh.physical_faces]]
    view = np.array([.35, -.94, .75])
    view /= np.linalg.norm(view)
    right = np.cross([0, 0, 1], view)
    right /= np.linalg.norm(right)
    project = np.column_stack((right, np.cross(view, right)))
    normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    light = np.array([-.3, -.5, .8])
    light /= np.linalg.norm(light)
    colors = (.65+.35*np.maximum(normals@light, 0))[:, None]*[.70, .84, .88]
    order = np.argsort(triangles.mean(axis=1)@view)
    fig, ax = plt.subplots(figsize=(11, 4), facecolor='white')
    fig.subplots_adjust(left=.03, right=.97, bottom=.05, top=.73)
    water = np.array([[-LENGTH/2-.5, -BEAM/2-.5, 0], [LENGTH/2+.5, -BEAM/2-.5, 0],
                      [LENGTH/2+.5, BEAM/2+.5, 0], [-LENGTH/2-.5, BEAM/2+.5, 0]])
    ax.add_collection(PolyCollection([water@project], facecolors='#edf6fa',
                                    edgecolors='#b6d5e3', linewidths=.7))
    ax.add_collection(PolyCollection(triangles[order]@project, facecolors=colors[order],
                                    edgecolors='#35576a', linewidths=.8))
    # Only draw intersections on viewer-facing hull panels.
    midpoint = result.waterline_segments.mean(axis=1)
    visible = midpoint[:, 1] < rest.position[1]
    ax.add_collection(LineCollection(result.waterline_segments[visible]@project,
                                    colors='#168bb9', linewidths=1.5))
    rim = vertices[np.r_[body.mesh.opening_rim, body.mesh.opening_rim[0]]]@project
    ax.plot(rim[:, 0], rim[:, 1], color='#173b50', linewidth=1.3)
    for point, label, offset in (
        (PAYLOAD_POSITION, f'{PAYLOAD_MASS:g} kg console / gear / people', (5, 35)),
        (MOTOR_POSITION, f'{MOTOR_MASS:g} kg motor', (-15, 35)),
    ):
        location = ((point-body.cg)@rotation.T+rest.position)@project
        ax.plot(*location, 'o', color='#d47924', markeredgecolor='white', markersize=7)
        ax.annotate(label, location, xytext=offset, textcoords='offset points',
                    ha='center', fontsize=10, color='#87450d',
                    arrowprops=dict(arrowstyle='-', color='#87450d'))
    ax.autoscale_view()
    ax.margins(x=.06, y=.12)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.suptitle('Simple 24 ft reference boat', y=.97, fontsize=18, color='#173b50')
    fig.text(.5, .88, f'24 ft x 9 ft | {SHELL_MASS:g} kg uniform shell + two point masses',
             ha='center', color='#526573', fontsize=11)
    fig.text(.5, .81, f'{body.mass:g} kg total | solved heave and pitch | artificial cap omitted',
             ha='center', color='#526573', fontsize=10)
    destination = Path(__file__).resolve().parents[1]/'docs/images/reference-boat.png'
    fig.savefig(destination, dpi=160, facecolor='white', bbox_inches='tight', pad_inches=.15)
    plt.close(fig)
    print(destination)


if __name__ == '__main__':
    main()

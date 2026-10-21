"""Render the canonical canoe at its solved calm-water heave equilibrium.

Run: uv run --with matplotlib python -m tests.plot_reference_hull
Matplotlib is needed only to regenerate the image, not to use Topwater.
"""

from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
import numpy as np

from topwater import SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from .reference_hull import PAYLOAD_MASS, PAYLOAD_POSITION, SEMI_AXES, input_data


def main():
    data = SimulationInputs.from_json(json.dumps(input_data(refinement=4)))
    rest, result = equilibrate_heave(data.body, data.initial_state)
    vertices = data.body.mesh.vertices-data.body.cg+rest.position
    rim = data.body.mesh.opening_rim
    # Omit the explicitly artificial closure so the physical opening is visible.
    faces = data.body.mesh.faces[data.body.mesh.physical_faces]
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0],
                       triangles[:, 2] - triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    light = np.array([-0.3, -0.5, 0.8])
    light /= np.linalg.norm(light)
    intensity = 0.65 + 0.35 * np.maximum(normals @ light, 0)
    colors = intensity[:, None] * np.array([0.70, 0.84, 0.88])
    colors[triangles.mean(axis=1)[:, 2] < 0] *= np.array([.75, .9, 1.])
    # An orthographic view drawn in 2-D avoids distorted axes and preserves
    # physical proportions. Sort triangles from back to front for occlusion.
    view = np.array([0.34, -0.94, 0.70])
    view /= np.linalg.norm(view)
    right = np.cross([0, 0, 1], view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)
    projected = triangles @ np.column_stack([right, up])
    order = np.argsort(triangles.mean(axis=1) @ view)
    fig = plt.figure(figsize=(10, 3.8), facecolor="white")
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.72])
    water = np.array([[-2, -.7, 0], [2, -.7, 0], [2, .7, 0], [-2, .7, 0]])
    ax.add_collection(PolyCollection([water @ np.column_stack([right, up])],
        facecolors='#edf6fa', edgecolors='#b6d5e3', linewidths=.7))
    ax.add_collection(PolyCollection(
        projected[order], facecolors=colors[order], edgecolors="#35576a", linewidths=0.4,
    ))
    # Only the near-side waterline is visible on this convex reference hull.
    normals_xy = result.waterline_segments.mean(axis=1)[:, :2]/SEMI_AXES[:2]**2
    visible = normals_xy @ view[:2] > 0
    ax.add_collection(LineCollection(
        result.waterline_segments[visible] @ np.column_stack([right, up]),
        colors='#168bb9', linewidths=1.5))
    outline = vertices[np.r_[rim, rim[0]]] @ np.column_stack([right, up])
    ax.plot(outline[:, 0], outline[:, 1], color="#173b50", linewidth=1.2)
    payload = (PAYLOAD_POSITION-data.body.cg+rest.position) @ np.column_stack([right, up])
    ax.plot(*payload, marker="o", markersize=7, color="#d47924", markeredgecolor="white")
    ax.annotate(f"{PAYLOAD_MASS:g} kg payload CG", payload, xytext=(20, 25),
                textcoords="offset points", color="#87450d", fontsize=10,
                arrowprops=dict(arrowstyle="-", color="#87450d"))
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")
    fig.suptitle("Reference canoe in calm water", fontsize=17, color="#173b50", y=0.98)
    a, b, c = SEMI_AXES
    fig.text(0.5, 0.84, f"Single-person reference · {2*a:g} m × {2*b:g} m × {c:g} m · artificial cap omitted",
             ha="center", color="#526573", fontsize=11)
    fig.text(.5, .76, f"Solved buoyancy = weight · {data.body.mass:.1f} kg · "
             f"{1000*(-vertices[:, 2].min()):.0f} mm draft · flat water at z = 0",
             ha='center', color='#526573', fontsize=10)
    destination = Path(__file__).resolve().parents[1] / "docs/images/reference-hull.png"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=160, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(destination)


if __name__ == "__main__":
    main()

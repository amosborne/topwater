"""Render the canonical geometry, without water or a simulated floating pose.

Run: uv run --with matplotlib python -m tests.plot_reference_hull
Matplotlib is needed only to regenerate the image, not to use Topwater.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from .reference_hull import PAYLOAD_MASS, PAYLOAD_POSITION, SEMI_AXES, mesh_arrays, opening_boundaries


def main():
    vertices, faces = mesh_arrays(refinement=3)
    rim, cap = opening_boundaries(vertices, faces)
    # Omit the explicitly artificial closure so the physical opening is visible.
    faces = np.delete(faces, cap, axis=0)
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0],
                       triangles[:, 2] - triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    light = np.array([-0.3, -0.5, 0.8])
    light /= np.linalg.norm(light)
    intensity = 0.65 + 0.35 * np.maximum(normals @ light, 0)
    colors = intensity[:, None] * np.array([0.70, 0.84, 0.88])
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
    ax.add_collection(PolyCollection(
        projected[order], facecolors=colors[order], edgecolors="#35576a", linewidths=0.4,
    ))
    outline = vertices[np.r_[rim, rim[0]]] @ np.column_stack([right, up])
    ax.plot(outline[:, 0], outline[:, 1], color="#173b50", linewidth=1.2)
    payload = PAYLOAD_POSITION @ np.column_stack([right, up])
    ax.plot(*payload, marker="o", markersize=7, color="#d47924", markeredgecolor="white")
    ax.annotate(f"{PAYLOAD_MASS:g} kg payload CG", payload, xytext=(20, 25),
                textcoords="offset points", color="#87450d", fontsize=10,
                arrowprops=dict(arrowstyle="-", color="#87450d"))
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")
    fig.suptitle("Reference hull mesh", fontsize=17, color="#173b50", y=0.98)
    a, b, c = SEMI_AXES
    fig.text(0.5, 0.84, f"Single-person reference · {2*a:g} m × {2*b:g} m × {c:g} m · artificial cap omitted",
             ha="center", color="#526573", fontsize=11)
    destination = Path(__file__).resolve().parents[1] / "docs/images/reference-hull.png"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=160, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(destination)


if __name__ == "__main__":
    main()

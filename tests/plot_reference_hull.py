"""Render the canonical geometry, without water or a simulated floating pose.

Run: uv run --with matplotlib python -m tests.plot_reference_hull
Matplotlib is needed only to regenerate the image, not to use Topwater.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from .reference_hull import mesh_arrays


def main():
    vertices, faces = mesh_arrays(refinement=3)
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
    view = np.array([0.34, -0.94, 0.36])
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
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")
    fig.suptitle("Reference hull mesh", fontsize=17, color="#173b50", y=0.98)
    fig.text(0.5, 0.84, "Triaxial ellipsoid · 6 m × 2 m × 1 m · 512 triangles",
             ha="center", color="#526573", fontsize=11)
    destination = Path(__file__).resolve().parents[1] / "docs/images/reference-hull.png"
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=160, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    print(destination)


if __name__ == "__main__":
    main()

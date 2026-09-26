#!/usr/bin/env python3
"""
Read data.json and write geometry.stl (in mm) for 3D printing.

Marching cubes on the exact signed distance field of the spheres, clipped to
the 28 x 28 x 54 mm box. Not mirrored.
"""

import json
import math
import numpy as np
import pyvista as pv
from scipy import ndimage

SPHERES_FILE = "data.json"
SCALE = 28.0                    # mm, data.json is a unit cube in x/y
LX, LY, LZ = 28.0, 28.0, 54.0   # mm, box size
H = 0.2                         # mm, marching-cubes grid spacing
PAD = 2                         # grid points outside the box on each side
STL_FILE = "geometry.stl"

# --- load spheres (mm) ---
data = json.load(open(SPHERES_FILE))
centers = np.array([s["center"] for s in data["spheres"]]) * SCALE
radii = np.array([s["radius"] for s in data["spheres"]]) * SCALE
print(f"Loaded {len(radii)} spheres from {SPHERES_FILE}")

# --- signed distance field (negative inside the solid) ---
x = (np.arange(math.ceil(LX / H) + 2 * PAD + 1) - PAD) * H
y = (np.arange(math.ceil(LY / H) + 2 * PAD + 1) - PAD) * H
z = (np.arange(math.ceil(LZ / H) + 2 * PAD + 1) - PAD) * H
f = np.full((len(x), len(y), len(z)), 1e3, dtype=np.float32)

for (cx, cy, cz), r in zip(centers, radii):
    # only compute the distance in a block around the sphere
    m = r + PAD * H
    i0, i1 = np.searchsorted(x, cx - m), np.searchsorted(x, cx + m, side="right")
    j0, j1 = np.searchsorted(y, cy - m), np.searchsorted(y, cy + m, side="right")
    k0, k1 = np.searchsorted(z, cz - m), np.searchsorted(z, cz + m, side="right")
    dist = np.sqrt(
        (x[i0:i1, None, None] - cx) ** 2
        + (y[None, j0:j1, None] - cy) ** 2
        + (z[None, None, k0:k1] - cz) ** 2
    ) - r
    block = f[i0:i1, j0:j1, k0:k1]
    np.minimum(block, dist, out=block)

# clip to the box
f = np.maximum(f, np.maximum(-x, x - LX)[:, None, None])
f = np.maximum(f, np.maximum(-y, y - LY)[None, :, None])
f = np.maximum(f, np.maximum(-z, z - LZ)[None, None, :]).astype(np.float32)
# grid points exactly on the surface give non-manifold marching-cubes edges
f[f == 0] = 1e-6

# --- marching cubes ---
grid = pv.ImageData(dimensions=f.shape, spacing=(H, H, H), origin=(x[0], y[0], z[0]))
grid.point_data["sdf"] = f.ravel(order="F")
mesh = grid.contour([0.0], scalars="sdf")

# normals must point outward (positive signed volume), slicers expect that
tri = mesh.faces.reshape(-1, 4)[:, 1:]
p = mesh.points.astype(np.float64)
signed_volume = np.einsum("ij,ij->i", p[tri[:, 0]], np.cross(p[tri[:, 1]], p[tri[:, 2]])).sum() / 6
if signed_volume < 0:
    mesh = mesh.flip_faces()

mesh.save(STL_FILE)

# --- checks ---
_, n_bodies = ndimage.label(f <= 0)
porosity = 1 - abs(mesh.volume) / (LX * LY * LZ)
print(f"Box:          {LX:g} x {LY:g} x {LZ:g} mm")
print(f"Triangles:    {mesh.n_cells}, open edges {mesh.n_open_edges}")
print(f"Solid bodies: {n_bodies}")
print(f"Porosity:     {porosity:.6f}")
print(f"Wrote {STL_FILE}")

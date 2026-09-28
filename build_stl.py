#!/usr/bin/env python3
"""
Read data.json and write geometry.stl (in mm) for 3D printing: the packing
clipped to the bore plus a solid wall, not mirrored. Bore, height and sphere
file are taken from build_geometry.py. Marching cubes on the exact signed
distance field of spheres + wall.
"""

import math
import numpy as np
import pyvista as pv
from scipy import ndimage

from build_geometry import HERE, SPHERES_FILE, SAMPLE_HEIGHT, BORE_RADIUS, load_spheres

MM = 1e3                        # meters -> mm
CX, CY = 0.0, 0.0                # bore axis
R_IN = BORE_RADIUS * MM         # bore radius
R_OUT = 14.5                    # mm, outer radius of the wall (29 mm diameter)
LZ = SAMPLE_HEIGHT * MM         # height
H = 0.2                         # mm, marching-cubes grid spacing
PAD = 2                         # grid points outside the part on each side
STL_FILE = HERE / "geometry.stl"


def signed_distance(centers, radii):
    """Grid axes and signed distance field (negative inside the solid)."""
    n_xy = math.ceil(2 * R_OUT / H) + 2 * PAD + 1
    x = CX - R_OUT + (np.arange(n_xy) - PAD) * H
    y = CY - R_OUT + (np.arange(n_xy) - PAD) * H
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

    # add the wall, then clip to the outer cylinder and to [0, LZ]
    r_xy = np.hypot(x[:, None] - CX, y[None, :] - CY)[:, :, None]
    f = np.minimum(f, np.maximum(R_IN - r_xy, r_xy - R_OUT))
    f = np.maximum(f, r_xy - R_OUT)
    f = np.maximum(f, np.maximum(-z, z - LZ)[None, None, :]).astype(np.float32)
    # grid points exactly on the surface give non-manifold marching-cubes edges
    f[f == 0] = 1e-6
    return x, y, z, f


def main():
    spheres = np.array(load_spheres(SPHERES_FILE)) * MM
    print(f"Loaded {len(spheres)} spheres from {SPHERES_FILE.name}")

    x, y, z, f = signed_distance(spheres[:, :3], spheres[:, 3])
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

    _, n_bodies = ndimage.label(f <= 0)
    v_bore = math.pi * R_IN**2 * LZ
    v_wall = math.pi * (R_OUT**2 - R_IN**2) * LZ
    porosity = 1 - (abs(mesh.volume) - v_wall) / v_bore
    print(f"Cylinder:     d_in {2 * R_IN:g} mm, d_out {2 * R_OUT:g} mm, height {LZ:g} mm")
    print(f"Triangles:    {mesh.n_cells}, open edges {mesh.n_open_edges}")
    print(f"Solid bodies: {n_bodies}")
    print(f"Porosity:     {porosity:.6f}  (bore only)")
    print(f"Wrote {STL_FILE.name}")


if __name__ == "__main__":
    main()

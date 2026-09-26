#!/usr/bin/env python3
"""
Read spheres.txt (cx cy cz r per line, meters) and write geometry.raw
in the format POREMAPS expects: uint8 voxels (0 = fluid, 1 = solid),
Fortran memory order (x varies fastest).

Grid size and voxel size MUST match the .inp file.
"""

import numpy as np
from pathlib import Path

# --- must match input.inp ---
NX, NY, NZ = 30, 30, 50
VOXEL_SIZE = 1e-5           # meters
MIRROR_Z   = True           # double along z by reflecting, so face z=0 == face z=Lz
# ----------------------------

HERE = Path(__file__).resolve().parent
SPHERES_FILE = HERE / "spheres.txt"
OUTPUT_FILE  = HERE / "geometry.raw"


def load_spheres(path):
    spheres = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        cx, cy, cz, r = (float(v) for v in line.split())
        spheres.append((cx, cy, cz, r))
    return spheres


def voxelize(spheres, nx, ny, nz, vs):
    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    x = (i + 0.5) * vs
    y = (j + 0.5) * vs
    z = (k + 0.5) * vs

    geom = np.zeros((nx, ny, nz), dtype=np.uint8)
    for cx, cy, cz, r in spheres:
        inside = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r ** 2
        geom[inside] = 1
    return geom


def main():
    spheres = load_spheres(SPHERES_FILE)
    print(f"Loaded {len(spheres)} spheres from {SPHERES_FILE.name}")

    geom = voxelize(spheres, NX, NY, NZ, VOXEL_SIZE)

    if MIRROR_Z:
        geom = np.concatenate([geom, geom[:, :, ::-1]], axis=2)
        nz_out = 2 * NZ
        print(f"Mirror-doubled along z: geometry is now {NX} x {NY} x {nz_out}")
    else:
        nz_out = NZ

    # numpy's tofile() always writes C-order; use tobytes(order='F') for true Fortran order
    OUTPUT_FILE.write_bytes(geom.tobytes(order='F'))

    total = NX * NY * nz_out
    solid = int((geom == 1).sum())
    print(f"Grid: {NX} x {NY} x {nz_out} = {total} voxels, vs = {VOXEL_SIZE} m")
    print(f"Solid voxels: {solid}  ({solid/total:.1%})")
    print(f"Porosity:     {(geom == 0).mean():.4f}")
    print(f"Wrote {OUTPUT_FILE.name} ({OUTPUT_FILE.stat().st_size} bytes)")
    if MIRROR_Z:
        print(f"\nMake sure input.inp has:  size_x_y_z  {NX} {NY} {nz_out}")


if __name__ == "__main__":
    main()

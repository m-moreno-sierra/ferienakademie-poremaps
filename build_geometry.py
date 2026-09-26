#!/usr/bin/env python3
"""
Read data.json, clip the packing to a cylindrical bore (everything outside is
solid), and write BOTH:
  - geometry.raw  (voxel grid for POREMAPS: uint8, Fortran memory order)
  - input.inp     (POREMAPS input file)

Both files are regenerated every run — do not edit them by hand; edit this
script and rerun.

JSON schema (see data.json for an example):
{
  "scale": s,        # optional, defaults to SPHERE_COORD_SCALE
  "spheres": [
    {"center": [cx, cy, cz], "radius": r},
    ...
  ]
}
Every center + radius is multiplied by scale to get meters. Any extra
top-level keys (e.g. "notes") are ignored.
"""

import json
import numpy as np
from pathlib import Path

# --- geometry grid ---
NX, NY, NZ = 70, 70, 135
VOXEL_SIZE = 4e-4              # meters
MIRROR_Z   = True              # double along z by reflecting, so face z=0 == face z=Lz

# --- cylindrical sample (axis along z, centered in the x/y domain) ---
# The packing fills the bore; the printed wall (r = 11 .. 14.5 mm) lies outside
# it, so every voxel outside the bore is simply solid.
CYL_CENTER  = (NX * VOXEL_SIZE / 2, NY * VOXEL_SIZE / 2)   # meters, (14 mm, 14 mm)
BORE_RADIUS = 11e-3            # meters, 22 mm inner diameter

# --- sphere coordinate scaling ---
# Multiplies every sphere center + radius from the JSON (result: meters).
# If the JSON has its own "scale" field, that overrides this constant.
# Default: NX * VOXEL_SIZE assumes JSON is in a unit cube [0, 1] fitted to the x-domain.
# Set to 1.0 if the JSON is already in meters.
SPHERE_COORD_SCALE = NX * VOXEL_SIZE

# --- POREMAPS solver settings (see README) ---
BOUNDARY_METHOD    = 0                    # 0 = periodic all around
MAX_ITER           = 100_000
IT_EVAL            = 100
IT_WRITE           = 100
SOLVING_ALGORITHM  = 2
EPS                = 1e-6
DOM_DECOMPOSITION  = (0, 0, 0)            # 0 0 0 = let MPI decide
DOM_INTEREST       = (0, 0, 0, 0, 0, 0)   # all zeros = whole domain
WRITE_OUTPUT       = (1, 1, 0, 0)         # velocity, pressure, neighborhood, decomp
LOG_FILE_NAME      = "permeability_spheres.log"

HERE          = Path(__file__).resolve().parent
SPHERES_FILE  = HERE / "data.json"
GEOMETRY_FILE = HERE / "geometry.raw"
INPUT_FILE    = HERE / "input.inp"


def load_spheres(path):
    data = json.loads(path.read_text())
    scale = data.get("scale", SPHERE_COORD_SCALE)   # JSON overrides script default
    return [(*(c * scale for c in s["center"]), s["radius"] * scale) for s in data["spheres"]]


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


def bore_mask(nx, ny, vs):
    """2D mask (nx, ny): True where the voxel center lies inside the bore."""
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    cx, cy = CYL_CENTER
    return ((i + 0.5) * vs - cx) ** 2 + ((j + 0.5) * vs - cy) ** 2 <= BORE_RADIUS ** 2


def write_input_file(nx, ny, nz, vs, porosity):
    dd = " ".join(str(x) for x in DOM_DECOMPOSITION)
    di = " ".join(str(x) for x in DOM_INTEREST)
    wo = " ".join(str(x) for x in WRITE_OUTPUT)
    INPUT_FILE.write_text(
        f"dom_decomposition {dd}\n"
        f"boundary_method {BOUNDARY_METHOD}\n"
        f"geometry_file_name {GEOMETRY_FILE.name}\n"
        f"size_x_y_z  {nx} {ny} {nz}\n"
        f"voxel_size  {vs}\n"
        f"max_iter    {MAX_ITER}\n"
        f"it_eval    {IT_EVAL}\n"
        f"it_write   {IT_WRITE}\n"
        f"log_file_name {LOG_FILE_NAME}\n"
        f"solving_algorithm {SOLVING_ALGORITHM}\n"
        f"eps {EPS}\n"
        f"porosity {porosity:.6f}\n"
        f"dom_interest {di}\n"
        f"write_output {wo}\n"
    )


def main():
    spheres = load_spheres(SPHERES_FILE)
    print(f"Loaded {len(spheres)} spheres from {SPHERES_FILE.name}")

    geom = voxelize(spheres, NX, NY, NZ, VOXEL_SIZE)
    mask = bore_mask(NX, NY, VOXEL_SIZE)
    geom[~mask] = 1                 # everything outside the bore is solid

    if MIRROR_Z:
        geom = np.concatenate([geom, geom[:, :, ::-1]], axis=2)
        nz_out = 2 * NZ
        print(f"Mirror-doubled along z: geometry is now {NX} x {NY} x {nz_out}")
    else:
        nz_out = NZ

    GEOMETRY_FILE.write_bytes(geom.tobytes(order="F"))

    total    = NX * NY * nz_out
    bore     = int(mask.sum()) * nz_out
    fluid    = int((geom == 0).sum())   # fluid only exists inside the bore
    porosity = fluid / bore             # fluid_volume / bore_volume

    print(f"Grid: {NX} x {NY} x {nz_out} = {total} voxels, vs = {VOXEL_SIZE} m")
    print(f"Solid voxels: {total - fluid}  ({(total - fluid)/total:.1%} of grid)")
    print(f"Bore voxels:  {bore}  (d = {2 * BORE_RADIUS * 1e3:g} mm)")
    print(f"Fluid voxels: {fluid}")
    print(f"Porosity:     {porosity:.6f}  (fluid_volume / bore_volume)")
    print(f"Wrote {GEOMETRY_FILE.name} ({GEOMETRY_FILE.stat().st_size} bytes)")

    write_input_file(NX, NY, nz_out, VOXEL_SIZE, porosity)
    print(f"Wrote {INPUT_FILE.name}")


if __name__ == "__main__":
    main()

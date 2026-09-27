# Sphere packing → POREMAPS permeability

Voxelize a sphere packing in a cylindrical sample, compute its permeability with [POREMAPS](https://git.rwth-aachen.de/david.krach/poremaps), and export the same sample as an STL for 3D printing.

```
data.json ─► build_geometry.py ─► geometry.raw + input.inp ─► POREMAPS ─► permeability_spheres.log + field .raw files ─► fields2vtu.py ─► .vtu (ParaView)
          └► build_stl.py ─► geometry.stl (3D printing)
```

## Files

| File | Purpose |
|---|---|
| `data.json` | Sphere packing used by default (see [Sample and data](#sample-and-data)). |
| `data_old.json` | Earlier packing, kept for reference. |
| `build_geometry.py` | Voxelizes the packing into `geometry.raw` (uint8, Fortran order) and writes the matching POREMAPS `input.inp`. Holds all grid, sample and solver settings. |
| `build_stl.py` | Writes `geometry.stl` in mm: the packing clipped to the bore plus a solid wall. Takes the sample dimensions from `build_geometry.py`. |
| `fields2vtu.py` | Combines the POREMAPS output fields into one `.vtu`. Verbatim copy from upstream POREMAPS (MIT, David Krach & Matthias Ruf). |

All generated files (`geometry.raw`, `input.inp`, `geometry.stl`, `*.log`, field `.raw` files, `*.vtu`) are git-ignored. Do not edit them by hand; rerun the scripts.

## Requirements

- POREMAPS binary (build from source) and an MPI runtime (`mpiexec`/`mpirun`)
- Python 3 with `numpy`; `scipy` + `pyvista` for the STL; `pyevtk` for the VTU; [ParaView](https://www.paraview.org/) to view it

## Run

```bash
python build_geometry.py                                    # geometry.raw + input.inp
python build_stl.py                                         # optional: geometry.stl
mpiexec -n 4 <path>/POREMAPS input.inp                      # mpirun -np 4 on Linux/macOS
python fields2vtu.py geometry.raw 70 70 270 4e-4            # optional: fields_geometry.vtu
```

The `fields2vtu.py` arguments are the grid size and voxel size that `build_geometry.py` prints. The grid is `NX NY 2*NZ` when `MIRROR_Z` is on.

## Sample and data

The sample is a cylinder with its axis along z, centered at (14, 14) mm, 54 mm tall:

- **Bore**, 22 mm diameter: holds the packing, and spheres are cut at its edge.
- **Wall**, 22 to 29 mm diameter: solid. It only exists in the STL; in the voxel grid everything outside the bore is simply solid.

Sphere file schema: `{"scale": s, "spheres": [{"center": [cx, cy, cz], "radius": r}, ...]}`. Centers and radii times `scale` give meters. Without `scale`, `NX * VOXEL_SIZE` (28 mm) is used, i.e. x/y are unit-normalized. Other keys are ignored.

- **`data.json`**: 177 beads of about 5 mm diameter, settled under gravity into the bore (soft-sphere energy minimization). Afterwards the radii were grown by 6 % so touching beads fuse through necks into one printable body. The bed was filled higher than 54 mm and cut at 54 mm, so the spheres at the top are cut on purpose. The `params` field records the generator settings.
- **`data_old.json`**: earlier 75-sphere packing without `scale`. To use it, set `SPHERES_FILE` in `build_geometry.py`; `build_stl.py` follows automatically.

Coordinates use a corner origin: the domain spans `[0, NX*vs] x [0, NY*vs] x [0, NZ*vs]`, and voxel `(i, j, k)` has its center at `((i+0.5)*vs, (j+0.5)*vs, (k+0.5)*vs)`. The `.vtu` from `fields2vtu.py` is centered on the domain instead (shifted by half the domain size).

## Configuration

All settings live at the top of `build_geometry.py`: grid (`NX`, `NY`, `NZ`, `VOXEL_SIZE`, `MIRROR_Z`), sample (`CYL_CENTER`, `BORE_RADIUS`), sphere file and solver (`BOUNDARY_METHOD`, `MAX_ITER`, `EPS`, `WRITE_OUTPUT`, ...). `build_stl.py` imports these and only defines its own wall radius and mesh resolution.

## Simulation notes

**Mirroring in z.** `MIRROR_Z = True` appends a mirrored copy along z, so the last slice equals the first and periodic flow in z has no jump at the wrap. The sample becomes twice as long. Because of the mirror symmetry, the off-diagonal components `wk13` and `wk23` are zero by construction; only `wk33` is meaningful.

**Boundary conditions** (`BOUNDARY_METHOD`):

| Value | z (flow) | x, y (sides) |
|:---:|---|---|
| 0 | periodic | periodic |
| 1 | periodic | slip |
| 2 | periodic | no-slip |
| 3 | non-periodic | slip |
| 4 | non-periodic | no-slip |

The default is 0. The sides are already closed by the solid voxels outside the bore, so periodic sides do not connect anything. For methods 3 and 4, turn `MIRROR_Z` off; the mirror only helps with periodic z.

**Output.** `permeability_spheres.log` has one row every `IT_WRITE` iterations. The `wk33` column of the last row is the permeability along z in m².

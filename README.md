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
| `data_old.json` | Earlier packing, kept for reference in its original format (unit-normalized to 28 mm, no `scale`). It does not fit the current scripts; its `notes` say how to convert it. |
| `build_geometry.py` | Voxelizes the packing into `geometry.raw` (uint8, Fortran order) and writes the matching POREMAPS `input.inp`. Holds all grid, sample and solver settings. |
| `build_stl.py` | Writes `geometry.stl` in mm: the packing clipped to the bore plus a solid wall. Takes the sample dimensions from `build_geometry.py`. |
| `fields2vtu.py` | Combines the POREMAPS output fields into one `.vtu`. Verbatim copy from upstream POREMAPS (MIT, David Krach & Matthias Ruf). |

All generated files (`geometry.raw`, `input.inp`, `geometry.stl`, `*.log`, field `.raw` files, `*.vtu`) are git-ignored in the project root; finished runs are copied to `results/<run>/` and committed there (see its README). Do not edit them by hand; rerun the scripts.

## Requirements

- POREMAPS binary (build from source) and an MPI runtime (`mpiexec`/`mpirun`)
- Python 3 with `numpy`; `scipy` + `pyvista` for the STL and `view_results.py`; `pyevtk` for the VTU; [ParaView](https://www.paraview.org/) to view it

## Run

```bash
python build_geometry.py                                    # geometry.raw + input.inp
python build_stl.py                                         # optional: geometry.stl
mpiexec -n 4 <path>/POREMAPS input.inp                      # mpirun -np 4 on Linux/macOS
python fields2vtu.py geometry.raw 76 76 360 3e-4            # optional: fields_geometry.vtu
python view_results.py [results/<run>]                      # optional: quick PyVista view of a saved run
```

The `fields2vtu.py` arguments are the grid size and voxel size that `build_geometry.py` prints. The grid is `NX NY 2*NZ` when `MIRROR_Z` is on.

## Sample and data

The sample is a cylinder with its axis along z, its axis at x = y = 0, floor at z = 0, 54 mm tall:

- **Bore**, 22 mm diameter: holds the packing, and spheres are cut at its edge.
- **Wall**, 22 to 29 mm diameter: solid. It only exists in the STL; in the voxel grid everything outside the bore is simply solid.

Sphere file schema: `{"scale": s, "spheres": [{"center": [cx, cy, cz], "radius": r}, ...]}`. Centers and radii times `scale` give meters. Both files are in mm (`"scale": 0.001`) with the bore axis at x = y = 0; without `scale`, `SPHERE_COORD_SCALE` (mm) is used. Other keys are ignored.

- **`data.json`**: 177 beads of about 5 mm diameter, settled under gravity into the bore (soft-sphere energy minimization). Afterwards the radii were grown by 6 % so touching beads fuse through necks into one printable body. The bed was filled higher than 54 mm and cut at 54 mm, so the spheres at the top are cut on purpose. The `params` field records the generator settings.

The voxel grid is cropped in x/y to the bore plus a solid rim of at least one voxel, centered on the bore axis: `NX = NY = ceil(22 mm / vs) + 2`, `NZ = 54 mm / vs`. Voxel `(i, j, k)` has its center at `GRID_ORIGIN + ((i+0.5)*vs, (j+0.5)*vs, (k+0.5)*vs)`, with `GRID_ORIGIN = (-NX*vs/2, -NY*vs/2, 0)`. The STL uses the same frame (in mm). The `.vtu` from `fields2vtu.py` is centered on the domain in x/y, so it matches too; in z it is centered as well, i.e. shifted down by half the (mirrored) length.

## Configuration

All settings live at the top of `build_geometry.py`: grid (`VOXEL_SIZE`, `MIRROR_Z`; `NX`, `NY`, `NZ` follow from the sample size), sample (`BORE_RADIUS`, `SAMPLE_HEIGHT`), sphere file and solver (`BOUNDARY_METHOD`, `MAX_ITER`, `EPS`, `WRITE_OUTPUT`, ...). `build_stl.py` imports these and only defines its own wall radius and mesh resolution.

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

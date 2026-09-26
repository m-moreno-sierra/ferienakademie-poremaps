# Sphere-pack example for POREMAPS

A minimal end-to-end example: define a few spheres in a text file, voxelize them, and compute the permeability with POREMAPS.

## Files

| File | Purpose |
|---|---|
| `spheres.txt` | Sphere data. One line per sphere: `cx cy cz r` (all in meters). Comments start with `#`. |
| `build_geometry.py` | Reads `spheres.txt`, voxelizes it, writes `geometry.raw` in the format POREMAPS expects (uint8, Fortran byte order). Optionally mirror-doubles along z. |
| `input.inp` | POREMAPS input file. Grid size and voxel size must match the constants in `build_geometry.py`. |

Files generated at runtime (not tracked in git):

| File | Produced by |
|---|---|
| `geometry.raw` | `build_geometry.py` |
| `permeability_spheres.log` | POREMAPS — convergence history + permeability values |
| `velx_geometry.raw`, `vely_geometry.raw`, `velz_geometry.raw` | POREMAPS — velocity components (float64), when `write_output` flag 1 is set |
| `press_geometry.raw` | POREMAPS — pressure (float64), when flag 2 is set |
| `voxel_neighborhood_geometry.raw` | POREMAPS — voxel neighbor classification (uint32), when flag 3 is set |
| `domain_decomp_geometry.raw` | POREMAPS — MPI rank ownership (uint32), when flag 4 is set |
| `fields_geometry.vtu` | `fields2vtu.py` — combined VTU for ParaView |

## Requirements

- **POREMAPS binary** — build from source. Upstream: [git.rwth-aachen.de/david.krach/poremaps](https://git.rwth-aachen.de/david.krach/poremaps).
- **MPI runtime** with `mpirun` on `PATH` (OpenMPI, MPICH, or your cluster's vendor MPI).
- **Python 3** with `numpy`.
- Optional (for visualization): `pyevtk` (`pip install pyevtk`) and [ParaView](https://www.paraview.org/).

## Setup

Set an environment variable pointing to your POREMAPS install directory — the folder containing `bin/POREMAPS` and `fields2vtu.py`:

```bash
export POREMAPS_DIR=/absolute/path/to/poremaps
```

Add that line to your `~/.bashrc` / `~/.zshrc` to make it persistent.

## Run

```bash
# 1. build the geometry
python build_geometry.py

# 2. run POREMAPS (adjust -np to your physical core count)
mpirun -np 4 "$POREMAPS_DIR/bin/POREMAPS" input.inp

# 3. convert output fields to VTU for ParaView (optional)
python "$POREMAPS_DIR/fields2vtu.py" geometry.raw 30 30 100 1e-5
```

Then open `fields_geometry.vtu` in ParaView.

## Configuration

Grid resolution and physical scale live in two places that MUST match:

**`build_geometry.py` (top of file):**
```python
NX, NY, NZ = 30, 30, 50
VOXEL_SIZE = 1e-5           # meters
MIRROR_Z   = True           # doubles nz along z via mirror reflection
```

**`input.inp`:**
```
size_x_y_z  30 30 100        # NX, NY, NZ (or 2*NZ if MIRROR_Z is True)
voxel_size  1e-05
```

If `MIRROR_Z = True`, use `2*NZ` in `size_x_y_z`. The `fields2vtu.py` call in Step 3 also needs the *effective* nz (`100` here, not `50`).

## Coordinate convention

Corner origin. The domain spans `[0, NX*vs] x [0, NY*vs] x [0, NZ*vs]`. Voxel `(i, j, k)` has its center at `((i+0.5)*vs, (j+0.5)*vs, (k+0.5)*vs)`. Sphere coordinates in `spheres.txt` use the same convention.

## The mirror-z trick

`MIRROR_Z = True` duplicates the geometry along z with a mirror reflection: the result is `NX x NY x (2*NZ)`, and the last slice equals the first slice by construction. This makes periodic boundary conditions in z (`boundary_method 0`) natural for geometries that are not intrinsically periodic — no discontinuity at the wrap. The physical sample is now twice as long in the flow direction; the permeability tensor is unchanged by the symmetry.

Turn it off (`MIRROR_Z = False`) if your geometry is already periodic, or if you want to use non-periodic boundary conditions (`boundary_method 3` or `4`).

## Boundary conditions

Set in `input.inp`:

| `boundary_method` | z (flow) | x, y (sides) |
|:---:|---|---|
| 0 | periodic | periodic |
| 1 | periodic | slip |
| 2 | periodic | no-slip |
| 3 | non-periodic | slip |
| 4 | non-periodic | no-slip |

Default is `0`. Use `4` for a realistic finite sample with rigid walls.

## Reading the output

`permeability_spheres.log` — one row per `it_write` iterations. The last row's `wk33` column is the permeability in m² along the pressure gradient direction (z). `wk13`, `wk23` are the off-diagonal components; rotate the geometry and re-run to get the other rows of the full permeability tensor.

## Modifying the example

- **Different sphere pack**: edit `spheres.txt`. No other changes needed.
- **Higher/lower resolution**: change `NX, NY, NZ` and `VOXEL_SIZE` in `build_geometry.py`, and update `size_x_y_z` + `voxel_size` in `input.inp` to match.
- **Different geometry (non-spheres)**: modify `voxelize()` in `build_geometry.py` to mark solid voxels however you want — the write path and mirror logic stay the same.

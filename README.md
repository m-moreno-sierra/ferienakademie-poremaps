# Sphere-pack example for POREMAPS

A minimal end-to-end example: define a few spheres in a text file, voxelize them, and compute the permeability with POREMAPS.

## Files

| File | Purpose |
|---|---|
| `data.json` | Sphere data as JSON. Schema: `{"scale": s, "spheres": [{"center": [cx, cy, cz], "radius": r}, ...]}`. Every center + radius is multiplied by `scale` to get meters — set `scale = 1.0` if your data is already in meters. Missing `scale` defaults to `SPHERE_COORD_SCALE = NX * VOXEL_SIZE` (unit-normalized in x/y). The shipped `data.json` has no `scale`: x/y span `[0, 1]`, z spans `[0, 54/28]`, i.e. a 28 x 28 x 54 mm column. Extra top-level keys (like `"notes"`) are ignored. |
| `build_geometry.py` | Reads `data.json`, voxelizes it, writes `geometry.raw` **and** `input.inp` in one shot. All grid/solver constants live at the top of this script — single source of truth. |
| `build_stl.py` | Optional. Reads `data.json` and writes `geometry.stl` for 3D printing: marching cubes on the exact signed distance field, in mm, clipped to the 28 x 28 x 54 mm box, not mirrored. Box size and resolution (0.2 mm) are hard-coded at the top. Needs `pyvista` + `scipy`. |
| `fields2vtu.py` | Reads POREMAPS's output `.raw` fields (velocity, pressure, geometry, etc.) and writes a single `.vtu` for ParaView. Copied verbatim from upstream POREMAPS (MIT license, David Krach & Matthias Ruf). |

Files generated at runtime (not tracked in git):

| File | Produced by |
|---|---|
| `geometry.raw` | `build_geometry.py` |
| `input.inp` | `build_geometry.py` — POREMAPS input file, regenerated each run from constants in the script (also contains the computed porosity = fluid_voxels / total_voxels) |
| `geometry.stl` | `build_stl.py` — printable part in mm (unmirrored) |
| `permeability_spheres.log` | POREMAPS — convergence history + permeability values |
| `velx_geometry.raw`, `vely_geometry.raw`, `velz_geometry.raw` | POREMAPS — velocity components (float64), when `write_output` flag 1 is set |
| `press_geometry.raw` | POREMAPS — pressure (float64), when flag 2 is set |
| `voxel_neighborhood_geometry.raw` | POREMAPS — voxel neighbor classification (uint32), when flag 3 is set |
| `domain_decomp_geometry.raw` | POREMAPS — MPI rank ownership (uint32), when flag 4 is set |
| `fields_geometry.vtu` | `fields2vtu.py` — combined VTU for ParaView |

## How it works

```
data.json
     │
     ▼
build_geometry.py  ──►  geometry.raw  +  input.inp
build_stl.py       ──►  geometry.stl  (optional, for 3D printing)
                                                │
                                                ▼
                                          POREMAPS  ──►  permeability_spheres.log
                                                    ──►  velx / vely / velz / press .raw files
                                                                                │
                                                                                ▼
                                                                      fields2vtu.py  ──►  fields_geometry.vtu  (ParaView)
```

1. **Describe the geometry** in `data.json` — a list of sphere centers and radii (unit-normalized, scaled to meters by `scale`).
2. **Run `python build_geometry.py`.** This does two things in one shot:
   - voxelizes the spheres into `geometry.raw` (uint8 voxels, Fortran memory order — POREMAPS's format).
   - writes `input.inp` with the matching grid size, voxel size, and the actual porosity computed as `fluid_voxels / total_voxels`.

   All grid + solver knobs (`NX`, `NY`, `NZ`, `VOXEL_SIZE`, `MIRROR_Z`, `BOUNDARY_METHOD`, `MAX_ITER`, `EPS`, `WRITE_OUTPUT`, …) live at the top of `build_geometry.py`. **This is the single source of truth.** `input.inp` is regenerated every run and should not be edited by hand.

3. **Run POREMAPS**: `mpirun -np <cores> "$POREMAPS_DIR/bin/POREMAPS" input.inp`. It solves the Stokes flow through your voxel geometry and writes:
   - `permeability_spheres.log` — convergence history and the permeability values in each iteration.
   - Optional `.raw` output fields (velocity, pressure, etc.) depending on the `WRITE_OUTPUT` flags.

4. **Optionally visualize**: `python fields2vtu.py geometry.raw <NX> <NY> <NZ> <VOXEL_SIZE>` bundles POREMAPS's output `.raw` files into one `.vtu`. Open it in ParaView.

## Requirements

- **POREMAPS binary** — build from source. Upstream: [git.rwth-aachen.de/david.krach/poremaps](https://git.rwth-aachen.de/david.krach/poremaps).
- **MPI runtime** with `mpirun` on `PATH` (OpenMPI, MPICH, or your cluster's vendor MPI).
- **Python 3** with `numpy`.
- Optional (for the STL): `pyvista`, `scipy`.
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
python build_stl.py        # optional: geometry.stl for 3D printing

# 2. run POREMAPS (adjust -np to your physical core count)
mpirun -np 4 "$POREMAPS_DIR/bin/POREMAPS" input.inp

# 3. convert output fields to VTU for ParaView (optional)
python fields2vtu.py geometry.raw 30 30 100 1e-5
```

Then open `fields_geometry.vtu` in ParaView.

## Configuration

All grid and POREMAPS solver settings live at the top of `build_geometry.py` — this is the **single source of truth**. Running the script regenerates both `geometry.raw` and `input.inp` consistently.

**Geometry:**
```python
NX, NY, NZ = 30, 30, 50
VOXEL_SIZE = 1e-5           # meters
MIRROR_Z   = True           # doubles nz along z via mirror reflection
```

**POREMAPS solver:**
```python
BOUNDARY_METHOD    = 0                    # 0 = periodic all around (see table below)
MAX_ITER           = 100_000
EPS                = 1e-6
WRITE_OUTPUT       = (1, 1, 0, 0)         # velocity, pressure, neighborhood, decomp
# ... etc — see the file
```

Do not edit `input.inp` by hand: `build_geometry.py` will overwrite it on the next run. The `fields2vtu.py` call in Step 3 needs the **effective** nz (`2*NZ` if `MIRROR_Z = True`).

## Coordinate convention

Corner origin. The domain spans `[0, NX*vs] x [0, NY*vs] x [0, NZ*vs]`. Voxel `(i, j, k)` has its center at `((i+0.5)*vs, (j+0.5)*vs, (k+0.5)*vs)`. Sphere coordinates in `data.json` (after scaling) use the same convention.

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

## Visualizing the results (`fields2vtu.py`)

Combines POREMAPS's `.raw` output files into a single `.vtu` you can open in ParaView.

**Usage** — 5 positional arguments:
```
python fields2vtu.py <geometry_file> <NX> <NY> <NZ> <VOXEL_SIZE>
```

For this example (mirror-doubled to 30×30×100 at 10 µm voxels):
```bash
python fields2vtu.py geometry.raw 30 30 100 1e-5
```

**When and where to run it:** after POREMAPS finishes, from the folder that contains `geometry.raw` and POREMAPS's output files. The script looks for them in the current working directory.

**What it reads:**
- `geometry.raw` (required — errors out if missing)
- Any of these that are present (missing ones are skipped with a printed note):
  - `press_geometry.raw`
  - `velx_geometry.raw`, `vely_geometry.raw`, `velz_geometry.raw`
  - `voxel_neighborhood_geometry.raw`
  - `domain_decomp_geometry.raw`

**What it writes:** `fields_geometry.vtu` in the same folder.

**In ParaView:** File → Open → pick the `.vtu` → click **Apply** → change **Representation** to *Point Gaussian* → in **Coloring** pick e.g. `z_velocity [m/s]` → click the rainbow **Rescale to Data Range** icon.

**Gotcha:** the `NZ` argument is the *effective* grid size on disk. With `MIRROR_Z = True` in `build_geometry.py`, that's `2*NZ` (e.g. `100`), not the original (`50`).

**Dependency:** `pyevtk` (`pip install pyevtk`). Without it you get `ModuleNotFoundError: No module named 'pyevtk'`.

## Modifying the example

- **Different sphere pack**: edit `data.json`. If the domain size changes, also update `NX, NY, NZ` in `build_geometry.py` and the box size in `build_stl.py`.
- **Higher/lower resolution**: change `NX, NY, NZ` and `VOXEL_SIZE` in `build_geometry.py`, and update `size_x_y_z` + `voxel_size` in `input.inp` to match.
- **Different geometry (non-spheres)**: modify `voxelize()` in `build_geometry.py` to mark solid voxels however you want — the write path and mirror logic stay the same.

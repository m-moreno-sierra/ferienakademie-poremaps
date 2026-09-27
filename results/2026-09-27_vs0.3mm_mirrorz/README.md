# POREMAPS run: 0.3 mm voxels, mirrored z

Permeability of the `data.json` bead packing (177 beads, ~5 mm) in the 22 mm bore, 54 mm tall.

## Result

**k_zz = 1.245 × 10⁻⁸ m² (≈ 12.6 kD)**, from the `wk33` column of the last row of `permeability_spheres.log`.

- Referred to the bore cross-section: POREMAPS multiplies the mean fluid velocity by the `porosity` from `input.inp`, which is the bore porosity (0.3647), not the whole-grid porosity.
- Converged: stopped at iteration 29,600 when `conv` fell below `eps = 1e-6`. k changed by 0.009 % over the last 5,000 iterations.
- `wk13` and `wk23` are 0.04 % and 0.02 % of `wk33`. They are zero by construction (z mirror), so this is numerical noise.
- The `k13`, `k23`, `k33` columns are 0 because no domain of interest was set.

For water (20 °C) this is a hydraulic conductivity of K = kρg/μ ≈ 0.12 m/s. Darcy's law only holds at bead Reynolds numbers ≲ 1, i.e. a head difference of about 0.1 mm of water over the 54 mm sample; use a viscous liquid for larger pressure differences.

Kozeny–Carman at φ = 0.365 gives 1.68e-8 m² (d = 5 mm) or 1.88e-8 m² (d = 5.3 mm, radii after the 6 % neck growth). The simulation is 25–35 % lower, plausible given the necks and the coarse resolution (narrowest gaps ≈ 3 voxels). A resolution study has not been done yet.

## Setup

| | |
|---|---|
| Code | commit `dac71c0` (`build_geometry.py`, `data.json`) |
| Grid | 76 × 76 × 360 (54 mm mirrored to 108 mm), voxel size 0.3 mm |
| Porosity | 0.364714 (bore only) |
| Boundary | method 0 (periodic everywhere; the bore is enclosed by solid voxels) |
| Solver | algorithm 2, `eps 1e-6`, `max_iter 100000` |
| Run | 8 MPI ranks, ~11.5 iterations/s, ~43 min |

## Files

| File | Content |
|---|---|
| `input.inp` | POREMAPS input |
| `geometry.raw` | voxel geometry (uint8, Fortran order, 1 = solid) |
| `permeability_spheres.log` | convergence log, one row per 100 iterations |
| `press_geometry.raw`, `vel{x,y,z}_geometry.raw` | pressure and velocity fields (float64, Fortran order, 16.6 MB each) |

To view the fields in ParaView, run from this folder:

```
python ../../fields2vtu.py geometry.raw 76 76 360 3e-4
```

To reproduce: check out the commit above, run `python build_geometry.py`, then `mpiexec -n 8 <path>/POREMAPS.exe input.inp`.

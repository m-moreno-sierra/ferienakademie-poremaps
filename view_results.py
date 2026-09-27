#!/usr/bin/env python3
"""
Quick PyVista view of a POREMAPS run in results/<run>/: geometry, v_z and pressure,
each cut open along x. Sizes and file names come from the run's input.inp.
Coordinates in mm, bore axis at x = y = 0 (same frame as build_geometry.py / build_stl.py).

Usage:
    python view_results.py [results/<run>] [--full] [--screenshot view.png]

Without a folder the newest run in results/ is shown. For a z-mirrored geometry only
the lower half (the real sample) is shown; --full shows the whole domain.
"""

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent

# (field name, file prefix) as written by POREMAPS, float64
FIELDS = [
    ("v_x [m/s]", "velx_"),
    ("v_y [m/s]", "vely_"),
    ("v_z [m/s]", "velz_"),
    ("pressure [-]", "press_"),
]


def read_inp(path):
    return {p[0]: p[1:] for p in (line.split() for line in path.read_text().splitlines()) if p}


def final_k(folder, params):
    log = folder / params["log_file_name"][0]
    if not log.exists():
        return None
    rows = [line for line in log.read_text().splitlines() if line and not line.startswith("#")]
    return float(rows[-1].split(",")[10]) if rows else None     # wk33


def build_grid(folder, full=False):
    params = read_inp(folder / "input.inp")
    shape = tuple(int(v) for v in params["size_x_y_z"])
    vs = float(params["voxel_size"][0]) * 1e3                   # mm
    geom_name = params["geometry_file_name"][0]

    def load(name, dtype):
        return np.fromfile(folder / name, dtype=dtype).reshape(shape, order="F")

    geom = load(geom_name, np.uint8)
    nz = shape[2]
    mirrored = nz % 2 == 0 and np.array_equal(geom[:, :, : nz // 2], geom[:, :, nz // 2 :][:, :, ::-1])
    kz = nz // 2 if mirrored and not full else nz
    print(f"grid {shape}, mirrored in z: {mirrored}, showing {kz} of {nz} layers")

    grid = pv.ImageData(dimensions=(shape[0] + 1, shape[1] + 1, kz + 1), spacing=(vs, vs, vs),
                        origin=(-shape[0] * vs / 2, -shape[1] * vs / 2, 0.0))
    grid.cell_data["solid"] = geom[:, :, :kz].ravel(order="F")
    for name, prefix in FIELDS:
        if (folder / f"{prefix}{geom_name}").exists():
            grid.cell_data[name] = load(f"{prefix}{geom_name}", np.float64)[:, :, :kz].ravel(order="F")
        else:
            print(f"missing {prefix}{geom_name}")
    return grid, final_k(folder, params)


def show(grid, title, screenshot=None):
    fluid = grid.threshold((0, 0), scalars="solid")
    solid = grid.threshold((1, 1), scalars="solid")

    p = pv.Plotter(shape=(1, 3), window_size=(1600, 700), off_screen=screenshot is not None)

    p.subplot(0, 0)
    p.add_text("geometry (solid, clipped)", font_size=10)
    p.add_mesh(solid.clip(normal="x"), color="lightgray")

    p.subplot(0, 1)
    if "v_z [m/s]" in grid.cell_data:
        p.add_text("v_z in fluid (clipped)", font_size=10)
        p.add_mesh(fluid.clip(normal="x"), scalars="v_z [m/s]", cmap="viridis",
                   scalar_bar_args=dict(fmt="%.1e", n_labels=3))
        p.add_mesh(fluid.slice(normal="z"), scalars="v_z [m/s]", cmap="viridis", show_scalar_bar=False)

    p.subplot(0, 2)
    if "pressure [-]" in grid.cell_data:
        p.add_text("pressure in fluid (clipped)", font_size=10)
        p.add_mesh(fluid.clip(normal="x"), scalars="pressure [-]", cmap="turbo",
                   scalar_bar_args=dict(n_labels=3))

    p.link_views()
    p.subplot(0, 0)
    p.add_text(title, position="lower_left", font_size=9)
    if screenshot:
        p.screenshot(screenshot)
        print(f"saved {screenshot}")
    else:
        p.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", type=Path, nargs="?", help="run folder (default: newest in results/)")
    ap.add_argument("--full", action="store_true", help="show the whole z-mirrored domain")
    ap.add_argument("--screenshot", help="save a PNG instead of opening a window")
    args = ap.parse_args()

    folder = args.folder or max((d for d in (HERE / "results").iterdir() if d.is_dir()),
                                key=lambda d: d.stat().st_mtime)
    grid, k = build_grid(folder, args.full)
    title = folder.name + (f"   k_zz = {k:.3e} m^2" if k is not None else "")
    show(grid, title, args.screenshot)

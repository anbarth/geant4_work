#!/usr/bin/env python3
"""Heatmap of the dose in the active layer of the radiochromic film, voxel by voxel.

Reads voxel_dose.csv written by the Sr90Film simulation (columns ix,iy,x_mm,y_mm,area_mm2,
edep_per_decay_MeV,dose_per_decay_Gy,rel_error,dose_rate_Gy_per_h_at_1MBq). Seen from above:
x runs left to right, y bottom to top, (0, 0) is the film centre.

Examples
    python plot_film_heatmap.py                              # build/voxel_dose.csv -> film_heatmap.png
    python plot_film_heatmap.py voxel_dose.csv -a 5 -t min   # dose rate in Gy/min for a 5 MBq source
    python plot_film_heatmap.py voxel_dose.csv -q dose_per_decay
    python plot_film_heatmap.py voxel_dose.csv --log         # log colour scale
    python plot_film_heatmap.py voxel_dose.csv --with-error  # add a map of the relative statistical error
    python plot_film_heatmap.py voxel_dose.csv --source-radius 8.2   # dashed outline of the source disc
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle

QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", "Gy / decay", False),
    "edep_per_decay": ("edep_per_decay_MeV", "MeV / decay", False),
}
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}


def si_prefix(vmax, unit):
    if not unit.startswith("Gy") or vmax <= 0:
        return 1.0, unit
    if vmax < 1e-3:
        return 1e6, "µ" + unit
    if vmax < 1:
        return 1e3, "m" + unit
    return 1.0, unit


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/voxel_dose.csv")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_rate")
    ap.add_argument("-a", "--activity", type=float, default=500.0, help="Sr-90 activity in MBq (dose_rate only)")
    ap.add_argument("-t", "--time-unit", choices=TIME_FACTORS, default="h")
    ap.add_argument("--log", action="store_true", help="logarithmic colour scale")
    ap.add_argument("--with-error", action="store_true", help="also draw the relative statistical error")
    ap.add_argument("--source-radius", type=float, help="draw the source disc outline (mm)")
    ap.add_argument("--source-x", type=float, default=0.0)
    ap.add_argument("--source-y", type=float, default=0.0)
    ap.add_argument("-o", "--output", default="film_heatmap.png")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    column, unit, is_rate = QUANTITIES[args.quantity]
    with open(args.csv, newline="") as f:
        rd = csv.DictReader(f)
        for n in ("ix", "iy", "x_mm", "y_mm", "rel_error", column):
            if n not in (rd.fieldnames or []):
                sys.exit(f"{args.csv}: column '{n}' not found. Found: {rd.fieldnames}")
        recs = [(int(r["ix"]), int(r["iy"]), float(r["x_mm"]), float(r["y_mm"]),
                 float(r[column]), float(r["rel_error"])) for r in rd]
    if not recs:
        sys.exit("No data rows.")
    a = np.array(recs)
    ix, iy = a[:, 0].astype(int), a[:, 1].astype(int)
    nx, ny = ix.max() + 1, iy.max() + 1
    vs = float(np.min(np.diff(np.unique(np.round(a[:, 2], 6)))))
    val = np.full((ny, nx), np.nan); err = np.full((ny, nx), np.nan)
    val[iy, ix] = a[:, 4]; err[iy, ix] = a[:, 5]
    x0, y0 = a[:, 2].min() - vs / 2, a[:, 3].min() - vs / 2
    extent = (x0, x0 + nx * vs, y0, y0 + ny * vs)

    if is_rate:
        val = val * args.activity * TIME_FACTORS[args.time_unit]
        unit = f"Gy/{args.time_unit}"
    scale, unit = si_prefix(np.nanmax(val), unit)
    val = val * scale
    label = {"dose_rate": "Dose rate", "dose_per_decay": "Dose per decay", "edep_per_decay": "Energy per decay"}[args.quantity]
    if is_rate:
        label += f" ({args.activity:g} MBq)"

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ncol = 2 if args.with_error else 1
    fig, axes = plt.subplots(1, ncol, figsize=(6.2 * ncol, 5.4), squeeze=False)
    positive = val[np.isfinite(val) & (val > 0)]
    if args.log:
        if positive.size == 0:
            sys.exit("No positive values for a log scale.")
        norm = LogNorm(vmin=positive.min(), vmax=positive.max())
    else:
        norm = Normalize(vmin=0, vmax=np.nanmax(val))
    ax = axes[0, 0]
    im = ax.imshow(val, origin="lower", extent=extent, cmap="viridis", norm=norm, interpolation="nearest")
    fig.colorbar(im, ax=ax, label=f"{label} [{unit}]", fraction=0.046, pad=0.04)
    ax.set_title(f"Active layer: {label}")
    if args.with_error:
        ax2 = axes[0, 1]
        im2 = ax2.imshow(np.where(val > 0, err, np.nan), origin="lower", extent=extent, cmap="magma",
                         interpolation="nearest")
        fig.colorbar(im2, ax=ax2, label="relative statistical error", fraction=0.046, pad=0.04)
        ax2.set_title("Statistical error")
    for a_ in axes[0]:
        a_.set_xlabel("x (mm)"); a_.set_ylabel("y (mm)")
        if args.source_radius:
            a_.add_patch(Circle((args.source_x, args.source_y), args.source_radius, fill=False,
                                ec="white", ls="--", lw=1))
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Wrote {args.output}  (grid {nx} x {ny}, voxel {vs:g} mm)")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Heatmap of the dose in the cell layer of one well, voxel by voxel.

Reads voxel_dose.csv written by the Sr90Plate simulation when fVoxelScoring is on
(columns: ix,iy,x_mm,y_mm,area_mm2,edep_per_decay_MeV,dose_per_decay_Gy,rel_error,
dose_rate_Gy_per_h_at_1MBq) and draws the well as seen from above:
x runs left to right (towards column 8), y bottom to top (towards row A).
Voxels outside the well are left blank and the well wall is drawn as a circle.

Examples
    python plot_voxel_heatmap.py                                # build/voxel_dose.csv -> voxel_dose_heatmap.png
    python plot_voxel_heatmap.py voxel_dose.csv -a 5 -t min     # dose rate in Gy/min for a 5 MBq source
    python plot_voxel_heatmap.py voxel_dose.csv --log           # log colour scale
    python plot_voxel_heatmap.py voxel_dose.csv --with-error    # add a map of the statistical error
    python plot_voxel_heatmap.py voxel_dose.csv --min-area 0.5  # hide edge voxels less than half inside the well
    python plot_voxel_heatmap.py voxel_dose.csv -q dose_per_decay --show
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle

# quantity name -> (csv column, unit label, scales with activity?)
QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", "Gy / decay", False),
    "edep_per_decay": ("edep_per_decay_MeV", "MeV / decay", False),
}
# The CSV dose rate is per hour; value in Gy/<unit> = value in Gy/h * factor
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}

NEEDED = ("ix", "iy", "x_mm", "y_mm", "area_mm2", "rel_error", "dose_per_decay_Gy")


def load(path, column):
    """Return dict of 2-D arrays indexed [iy, ix] plus the voxel size (mm)."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in NEEDED + (column,) if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"{path} is missing columns {missing}. Found: {reader.fieldnames}")
        recs = list(reader)
    if not recs:
        sys.exit(f"{path} has no rows.")

    ix = np.array([int(r["ix"]) for r in recs])
    iy = np.array([int(r["iy"]) for r in recs])
    n = int(max(ix.max(), iy.max())) + 1
    grids = {name: np.full((n, n), np.nan) for name in ("x", "y", "area", "value", "rel", "dose")}
    for r, i, j in zip(recs, ix, iy):
        grids["x"][j, i] = float(r["x_mm"])
        grids["y"][j, i] = float(r["y_mm"])
        grids["area"][j, i] = float(r["area_mm2"])
        grids["value"][j, i] = float(r[column])
        grids["rel"][j, i] = float(r["rel_error"])
        grids["dose"][j, i] = float(r["dose_per_decay_Gy"])

    # voxel size from the spacing of the centres (the CSV only lists voxels overlapping the well)
    xs = np.unique(np.round(grids["x"][~np.isnan(grids["x"])], 6))
    size = float(np.min(np.diff(xs))) if len(xs) > 1 else 1.0
    return grids, n, size


def pick_unit(values, unit):
    """For Gy/<time>, switch to mGy/<time> or µGy/<time> so the numbers are readable."""
    if not unit.startswith("Gy/"):
        return values, unit, 1.0
    per = unit[len("Gy"):]  # e.g. "/min"
    vmax = np.nanmax(values)
    if not np.isfinite(vmax) or vmax == 0:
        return values, unit, 1.0
    if vmax < 1e-3:
        return values * 1e6, "µGy" + per, 1e6
    if vmax < 1:
        return values * 1e3, "mGy" + per, 1e3
    return values, unit, 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/voxel_dose.csv", help="input CSV (default: build/voxel_dose.csv)")
    ap.add_argument("-o", "--output", default="voxel_dose_heatmap.png", help="output image (png/pdf/svg)")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_rate",
                    help="what to plot (default: dose_rate)")
    ap.add_argument("-a", "--activity", type=float, default=1.0,
                    help="Sr-90 activity in MBq; dose_rate only (CSV is for 1 MBq). Default 1")
    ap.add_argument("-t", "--time-unit", choices=TIME_FACTORS, default="h",
                    help="time unit for dose rate: h (hours), min (minutes) or s (seconds). Default h")
    ap.add_argument("--well-diameter", type=float, default=11.0,
                    help="well diameter in mm, for the outline (default 11.0)")
    ap.add_argument("--min-area", type=float, default=0.0, metavar="FRAC",
                    help="hide voxels with less than this fraction of a full voxel inside the well (0-1). Default 0")
    ap.add_argument("--max-rel-error", type=float, default=None, metavar="FRAC",
                    help="hide voxels whose statistical error is larger than this fraction (e.g. 0.2 = 20%%)")
    ap.add_argument("--with-error", action="store_true", help="add a second map with the statistical error")
    ap.add_argument("--log", action="store_true", help="logarithmic colour scale (empty voxels shown as grey)")
    ap.add_argument("--interp", default="nearest", help="matplotlib interpolation (default nearest = sharp voxels)")
    ap.add_argument("--cmap", default="inferno", help="matplotlib colormap (default: inferno)")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    args = ap.parse_args()

    if not args.show:
        matplotlib.use("Agg")

    column, unit, scales_with_activity = QUANTITIES[args.quantity]
    g, n, size = load(args.csv, column)
    values = g["value"].copy()
    if scales_with_activity:  # i.e. the dose rate
        values = values * args.activity * TIME_FACTORS[args.time_unit]
        unit = f"Gy/{args.time_unit}"

    # Area-weighted mean over the whole layer (equal density, so also mass-weighted). With layer-only
    # scoring this should agree with the dose of this well in well_dose.csv.
    valid = ~np.isnan(values) & (g["area"] > 0)
    mean_all = np.nansum(values[valid] * g["area"][valid]) / np.nansum(g["area"][valid])

    # voxels to hide
    hide = ~valid
    if args.min_area > 0:
        hide |= g["area"] < args.min_area * size * size
    if args.max_rel_error is not None:
        hide |= (g["rel"] > args.max_rel_error) | (g["rel"] == 0)
    shown = np.where(hide, np.nan, values)

    shown, unit_label, factor = pick_unit(shown, unit)
    mean_all *= factor

    vis = shown[~np.isnan(shown)]
    if vis.size == 0 or not np.any(vis > 0):
        print("Warning: no voxel has a positive dose, so the heatmap will be empty. "
              "Run more events, or check the generator and the voxel well.", file=sys.stderr)

    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("0.88")  # hidden / outside the well
    if args.log and np.any(vis > 0):
        norm = LogNorm(vmin=vis[vis > 0].min(), vmax=vis.max())
        shown = np.ma.masked_where(~(shown > 0), shown)
    else:
        vmax = float(np.nanmax(vis)) if vis.size else 1.0
        norm = Normalize(vmin=0, vmax=vmax if vmax > 0 else 1.0)
        shown = np.ma.masked_invalid(shown)

    xmin, xmax = np.nanmin(g["x"]) - size / 2, np.nanmax(g["x"]) + size / 2
    ymin, ymax = np.nanmin(g["y"]) - size / 2, np.nanmax(g["y"]) + size / 2
    extent = (xmin, xmax, ymin, ymax)

    title = {"dose_rate": "Dose rate", "dose_per_decay": "Dose per decay",
             "edep_per_decay": "Energy deposited per decay"}[args.quantity]
    if scales_with_activity:
        title += f" at {args.activity:g} MBq Sr-90"

    ncols = 2 if args.with_error else 1
    fig, axes = plt.subplots(1, ncols, figsize=(7.2 * ncols, 6.4), squeeze=False)
    ax = axes[0, 0]

    def finish(ax, im, label, title_text):
        ax.add_patch(Circle((0, 0), args.well_diameter / 2, fill=False, ec="white", lw=1.2))
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm, towards column 8)")
        ax.set_ylabel("y (mm, towards row A)")
        ax.set_title(title_text, fontsize=12)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(label)

    # origin="lower": row iy=0 is the most negative y, so +y (towards row A) points up
    im = ax.imshow(shown, cmap=cmap, norm=norm, origin="lower", extent=extent, interpolation=args.interp)
    finish(ax, im, unit_label, title)

    if args.with_error:
        err = np.where(hide | (g["rel"] == 0), np.nan, g["rel"] * 100.0)
        errshown = np.ma.masked_invalid(err)
        cmap_e = plt.get_cmap("viridis_r").copy()
        cmap_e.set_bad("0.88")
        emax = float(np.nanmax(err)) if np.any(~np.isnan(err)) else 100.0
        im2 = axes[0, 1].imshow(errshown, cmap=cmap_e, vmin=0, vmax=min(emax, 100.0),
                                origin="lower", extent=extent, interpolation=args.interp)
        finish(axes[0, 1], im2, "statistical error (%)", "Statistical error (1 sigma)")

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Wrote {args.output}")
    n_shown = int(np.sum(~np.isnan(np.asarray(shown.filled(np.nan) if np.ma.isMaskedArray(shown) else shown))))
    print(f"  voxels shown: {n_shown} of {int(np.sum(g['area'] > 0))}  (size {size:g} mm)")
    print(f"  layer mean (area-weighted): {mean_all:.4g} {unit_label}")
    if vis.size:
        print(f"  voxel range: {vis.min():.4g} .. {vis.max():.4g} {unit_label}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

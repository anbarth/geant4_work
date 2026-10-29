#!/usr/bin/env python3
"""Heatmaps of the dose in the cell layer of one or more wells, voxel by voxel.

Reads voxel_dose.csv written by the Sr90Plate simulation when fVoxelScoring is on
(columns: well,ix,iy,x_mm,y_mm,area_mm2,edep_per_decay_MeV,dose_per_decay_Gy,rel_error,
dose_rate_Gy_per_h_at_1MBq; a file without the 'well' column is treated as a single well).
Each well is drawn as seen from above: x runs left to right (towards column 8) and y bottom
to top (towards row A). Voxels outside the well are left blank, the well wall is a circle.
All wells share one colour scale so they can be compared directly.

Layouts
    panels : one small map per well
    plate  : every well drawn at its place on the 6 x 8 plate (best for many wells)
    auto   : one well -> single map, up to 6 wells -> panels, more -> plate  (default)

Examples
    python plot_voxel_heatmap.py                                  # build/voxel_dose.csv -> voxel_dose_heatmap.png
    python plot_voxel_heatmap.py voxel_dose.csv --wells C5,C4     # only these wells
    python plot_voxel_heatmap.py voxel_dose.csv --layout plate    # whole plate mosaic
    python plot_voxel_heatmap.py voxel_dose.csv -a 5 -t min       # dose rate in Gy/min for a 5 MBq source
    python plot_voxel_heatmap.py voxel_dose.csv --log             # log colour scale
    python plot_voxel_heatmap.py voxel_dose.csv --with-error      # also draw the statistical error
    python plot_voxel_heatmap.py voxel_dose.csv --min-area 0.5 --max-rel-error 0.3
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Circle

ROWS, COLS = 6, 8
ROW_NAMES = "ABCDEF"

# quantity name -> (csv column, unit label, scales with activity?)
QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", "Gy / decay", False),
    "edep_per_decay": ("edep_per_decay_MeV", "MeV / decay", False),
}
# The CSV dose rate is per hour; value in Gy/<unit> = value in Gy/h * factor
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}

NEEDED = ("ix", "iy", "x_mm", "y_mm", "area_mm2", "rel_error")


def well_position(name, pitch):
    """Centre (x, y) in mm of a well such as 'C5' on the plate (+x towards column 8, +y towards row A)."""
    row = ROW_NAMES.index(name[0])
    col = int(name[1:]) - 1
    return (col - (COLS - 1) / 2.0) * pitch, ((ROWS - 1) / 2.0 - row) * pitch


def load(path, column):
    """Return ({well: dict of 2-D arrays indexed [iy, ix]}, voxels per side, voxel size in mm)."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in NEEDED + (column,) if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"{path} is missing columns {missing}. Found: {reader.fieldnames}")
        has_well = "well" in reader.fieldnames
        recs = list(reader)
    if not recs:
        sys.exit(f"{path} has no rows.")

    n = int(max(max(int(r["ix"]), int(r["iy"])) for r in recs)) + 1
    wells = {}
    for r in recs:
        name = r["well"].strip().upper() if has_well else "?"
        g = wells.setdefault(name, {k: np.full((n, n), np.nan) for k in ("x", "y", "area", "value", "rel")})
        j, i = int(r["iy"]), int(r["ix"])
        g["x"][j, i] = float(r["x_mm"])
        g["y"][j, i] = float(r["y_mm"])
        g["area"][j, i] = float(r["area_mm2"])
        g["value"][j, i] = float(r[column])
        g["rel"][j, i] = float(r["rel_error"])

    first = next(iter(wells.values()))
    xs = np.unique(np.round(first["x"][~np.isnan(first["x"])], 6))
    size = float(np.min(np.diff(xs))) if len(xs) > 1 else 1.0
    return wells, n, size


def pick_unit(values, unit):
    """For Gy/<time>, switch to mGy/<time> or µGy/<time> so the numbers are readable."""
    if not unit.startswith("Gy/"):
        return unit, 1.0
    per = unit[len("Gy"):]  # e.g. "/min"
    vmax = np.nanmax(values) if values.size else 0
    if not np.isfinite(vmax) or vmax == 0:
        return unit, 1.0
    if vmax < 1e-3:
        return "µGy" + per, 1e6
    if vmax < 1:
        return "mGy" + per, 1e3
    return unit, 1.0


def draw_well(ax, arr, spec, extent, offset=(0.0, 0.0), well_diameter=11.0):
    """Draw one well's voxel map on ax. extent = (xmin, xmax, ymin, ymax) relative to the well axis."""
    ox, oy = offset
    im = ax.imshow(arr, cmap=spec["cmap"], norm=spec["norm"], origin="lower",   # +y (row A) points up
                   extent=(extent[0] + ox, extent[1] + ox, extent[2] + oy, extent[3] + oy),
                   interpolation=spec["interp"])
    ax.add_patch(Circle((ox, oy), well_diameter / 2, fill=False, ec="white", lw=1.0))
    return im


def render(specs, wells, names, extent, args, layout):
    """Return a figure showing one or two specs (dose and/or error)."""
    if len(names) == 1 and layout == "single":
        fig, axes = plt.subplots(1, len(specs), figsize=(7.2 * len(specs), 6.4), squeeze=False)
        for ax, spec in zip(axes[0], specs):
            im = draw_well(ax, spec["arrs"][names[0]], spec, extent, well_diameter=args.well_diameter)
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[2], extent[3])
            ax.set_aspect("equal")
            ax.set_xlabel("x (mm, towards column 8)")
            ax.set_ylabel("y (mm, towards row A)")
            ax.set_title(spec["title"] + (f"  -  well {names[0]}" if names[0] != "?" else ""), fontsize=12)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(spec["label"])
        fig.tight_layout()
        return fig

    spec = specs[0]
    if layout == "plate":
        pos = {w: well_position(w, args.pitch) for w in names}
        xs = [p[0] for p in pos.values()]
        ys = [p[1] for p in pos.values()]
        pad = args.pitch / 2
        fig, ax = plt.subplots(figsize=(11, 8))
        im = None
        for w in names:
            im = draw_well(ax, spec["arrs"][w], spec, extent, pos[w], args.well_diameter)
            ax.text(pos[w][0], pos[w][1] + args.well_diameter / 2 + 0.3, w, ha="center", va="bottom", fontsize=7)
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad + 1.0)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm from plate centre, towards column 8)")
        ax.set_ylabel("y (mm from plate centre, towards row A)")
        ax.set_title(spec["title"], fontsize=13)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(spec["label"])
        fig.tight_layout()
        return fig

    # panels
    ncols = min(len(names), 3)
    nrows = math.ceil(len(names) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols + 1.2, 4.3 * nrows), squeeze=False)
    im = None
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, w in zip(axes.ravel(), names):
        ax.set_visible(True)
        im = draw_well(ax, spec["arrs"][w], spec, extent, well_diameter=args.well_diameter)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("equal")
        ax.set_title(f"well {w}", fontsize=11)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
    fig.suptitle(spec["title"], fontsize=13)
    fig.colorbar(im, ax=[a for a in axes.ravel() if a.get_visible()], fraction=0.04, pad=0.02).set_label(spec["label"])
    return fig


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
    ap.add_argument("-w", "--wells", help="comma-separated wells to draw, e.g. C5,C4 (default: all in the file)")
    ap.add_argument("--layout", choices=("auto", "panels", "plate"), default="auto", help="see above (default auto)")
    ap.add_argument("--pitch", type=float, default=13.0, help="well pitch in mm for --layout plate (default 13.0)")
    ap.add_argument("--well-diameter", type=float, default=11.0,
                    help="well diameter in mm, for the outline (default 11.0)")
    ap.add_argument("--min-area", type=float, default=0.0, metavar="FRAC",
                    help="hide voxels with less than this fraction of a full voxel inside the well (0-1). Default 0")
    ap.add_argument("--max-rel-error", type=float, default=None, metavar="FRAC",
                    help="hide voxels whose statistical error is larger than this fraction (e.g. 0.2 = 20%%)")
    ap.add_argument("--with-error", action="store_true",
                    help="also draw the statistical error (side by side for one well, else a second image *_error.*)")
    ap.add_argument("--log", action="store_true", help="logarithmic colour scale (empty voxels shown as grey)")
    ap.add_argument("--interp", default="nearest", help="matplotlib interpolation (default nearest = sharp voxels)")
    ap.add_argument("--cmap", default="inferno", help="matplotlib colormap (default: inferno)")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    args = ap.parse_args()

    if not args.show:
        matplotlib.use("Agg")

    column, unit, scales_with_activity = QUANTITIES[args.quantity]
    wells, n, size = load(args.csv, column)

    names = list(wells)
    if args.wells:
        wanted = [w.strip().upper() for w in args.wells.split(",") if w.strip()]
        unknown = [w for w in wanted if w not in wells]
        if unknown:
            sys.exit(f"Wells not in {args.csv}: {', '.join(unknown)}. Available: {', '.join(names)}")
        names = wanted
    if "?" in names and len(names) > 1:
        sys.exit("The CSV has no 'well' column, so only one well can be drawn.")

    layout = args.layout
    if layout == "auto":
        layout = "single" if len(names) == 1 else ("panels" if len(names) <= 6 else "plate")
    elif len(names) == 1:
        layout = "single"
    if layout == "plate" and "?" in names:
        sys.exit("--layout plate needs the 'well' column in the CSV.")

    # ---- per-well values, hidden voxels, one common unit for everything ----
    raw, hide, mean_by_well = {}, {}, {}
    for w in names:
        g = wells[w]
        v = g["value"].copy()
        if scales_with_activity:  # i.e. the dose rate
            v = v * args.activity * TIME_FACTORS[args.time_unit]
        valid = ~np.isnan(v) & (g["area"] > 0)
        # Area-weighted mean over the layer (equal density, so also mass-weighted). With layer-only
        # scoring it should agree with this well's dose in well_dose.csv.
        mean_by_well[w] = float(np.nansum(v[valid] * g["area"][valid]) / np.nansum(g["area"][valid]))
        h = ~valid
        if args.min_area > 0:
            h |= g["area"] < args.min_area * size * size
        if args.max_rel_error is not None:
            h |= (g["rel"] > args.max_rel_error) | (g["rel"] == 0)
        raw[w], hide[w] = v, h

    if scales_with_activity:
        unit = f"Gy/{args.time_unit}"
    allvals = np.concatenate([np.where(hide[w], np.nan, raw[w]).ravel() for w in names])
    allvals = allvals[~np.isnan(allvals)]
    unit_label, factor = pick_unit(allvals, unit)
    allvals = allvals * factor

    if allvals.size == 0 or not np.any(allvals > 0):
        print("Warning: no voxel has a positive dose, so the heatmap will be empty. "
              "Run more events, or check the generator and the voxel wells.", file=sys.stderr)

    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("0.88")  # hidden / outside the well

    arrs = {}
    for w in names:
        a = np.where(hide[w], np.nan, raw[w] * factor)
        arrs[w] = np.ma.masked_where(~(a > 0), a) if args.log else np.ma.masked_invalid(a)
    if args.log and np.any(allvals > 0):
        norm = LogNorm(vmin=allvals[allvals > 0].min(), vmax=allvals.max())
    else:
        vmax = float(allvals.max()) if allvals.size else 1.0
        norm = Normalize(vmin=0, vmax=vmax if vmax > 0 else 1.0)

    title = {"dose_rate": "Dose rate", "dose_per_decay": "Dose per decay",
             "edep_per_decay": "Energy deposited per decay"}[args.quantity]
    if scales_with_activity:
        title += f" at {args.activity:g} MBq Sr-90"
    if layout in ("plate", "panels"):
        title += "  -  cell layer, per voxel"

    half = size * n / 2
    extent = (-half, half, -half, half)  # relative to the well axis (the grid is centred on it)
    dose_spec = dict(arrs=arrs, norm=norm, cmap=cmap, interp=args.interp, label=unit_label, title=title)
    specs = [dose_spec]

    err_spec = None
    if args.with_error:
        earrs, emax = {}, 0.0
        for w in names:
            e = np.where(hide[w] | (wells[w]["rel"] == 0), np.nan, wells[w]["rel"] * 100.0)
            earrs[w] = np.ma.masked_invalid(e)
            if np.any(~np.isnan(e)):
                emax = max(emax, float(np.nanmax(e)))
        ecmap = plt.get_cmap("viridis_r").copy()
        ecmap.set_bad("0.88")
        err_spec = dict(arrs=earrs, norm=Normalize(vmin=0, vmax=min(emax, 100.0) or 100.0), cmap=ecmap,
                        interp=args.interp, label="statistical error (%)", title="Statistical error (1 sigma)")

    out_files = []
    if err_spec is not None and layout == "single":
        fig = render([dose_spec, err_spec], wells, names, extent, args, layout)
        fig.savefig(args.output, dpi=200)
        out_files.append(args.output)
    else:
        fig = render([dose_spec], wells, names, extent, args, layout)
        fig.savefig(args.output, dpi=200)
        out_files.append(args.output)
        if err_spec is not None:
            stem, ext = os.path.splitext(args.output)
            err_out = f"{stem}_error{ext}"
            fig2 = render([err_spec], wells, names, extent, args, layout)
            fig2.savefig(err_out, dpi=200)
            out_files.append(err_out)

    for f in out_files:
        print(f"Wrote {f}")
    shown = sum(int(np.sum(~np.ma.getmaskarray(arrs[w]))) for w in names)
    print(f"  wells: {len(names)}   voxels shown: {shown}   voxel size {size:g} mm   layout: {layout}")
    if allvals.size:
        print(f"  voxel range: {allvals.min():.4g} .. {allvals.max():.4g} {unit_label}")
    if len(names) <= 12:
        for w in names:
            print(f"  well {w}: layer mean (area-weighted) {mean_by_well[w] * factor:.4g} {unit_label}")
    else:
        means = np.array([mean_by_well[w] * factor for w in names])
        print(f"  layer means over wells: {means.min():.4g} .. {means.max():.4g} {unit_label}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

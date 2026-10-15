#!/usr/bin/env python3
"""Heatmap of the per-well dose written by the Sr90Plate simulation.

Reads well_dose.csv (columns: well,row,col,edep_per_decay_MeV,dose_per_decay_Gy,
rel_error,dose_rate_Gy_per_h_at_1MBq) and draws the 6 x 8 plate layout
(rows A-F top to bottom, columns 1-8 left to right).

Examples
    python plot_heatmap.py                              # build/well_dose.csv -> well_dose_heatmap.png
    python plot_heatmap.py well_dose.csv --activity 5   # dose rate for a 5 MBq source
    python plot_heatmap.py well_dose.csv -t min         # dose rate per minute (also: s, h)
    python plot_heatmap.py well_dose.csv --log          # log colour scale (good for wide dynamic range)
    python plot_heatmap.py well_dose.csv --quantity dose_per_decay --show
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

ROWS, COLS = 6, 8
ROW_NAMES = "ABCDEF"

# quantity name -> (csv column, scale factor from csv units, unit label, scales with activity?)
QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", 1.0, "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", 1.0, "Gy / decay", False),
    "edep_per_decay": ("edep_per_decay_MeV", 1.0, "MeV / decay", False),
}


def load(path, column):
    values = np.full((ROWS, COLS), np.nan)
    errors = np.full((ROWS, COLS), np.nan)
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if column not in (reader.fieldnames or []):
            sys.exit(f"Column '{column}' not found in {path}. Found: {reader.fieldnames}")
        for rec in reader:
            r, c = int(rec["row"]), int(rec["col"])
            values[r, c] = float(rec[column])
            errors[r, c] = float(rec["rel_error"])
    if np.isnan(values).any():
        sys.exit(f"{path} does not contain all {ROWS * COLS} wells.")
    return values, errors


# The CSV dose rate is per hour; value in Gy/<unit> = value in Gy/h * factor
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}


def pick_unit(values, base_unit):
    """For Gy/<time>, switch to mGy/<time> or µGy/<time> so the numbers are readable."""
    if not base_unit.startswith("Gy/"):
        return values, base_unit
    per = base_unit[len("Gy"):]  # e.g. "/min"
    vmax = np.nanmax(values)
    if vmax == 0:
        return values, base_unit
    if vmax < 1e-3:
        return values * 1e6, "µGy" + per
    if vmax < 1:
        return values * 1e3, "mGy" + per
    return values, base_unit


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/well_dose.csv", help="input CSV (default: build/well_dose.csv)")
    ap.add_argument("-o", "--output", default="well_dose_heatmap.png", help="output image (png/pdf/svg)")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_per_decay",
                    help="what to plot (default: dose_per_decay)")
    ap.add_argument("-a", "--activity", type=float, default=1.0,
                    help="Sr-90 activity in MBq; dose_rate only (CSV is for 1 MBq). Default 1")
    ap.add_argument("-t", "--time-unit", choices=TIME_FACTORS, default="h",
                    help="time unit for dose rate: h (hours), min (minutes) or s (seconds). Default h")
    ap.add_argument("--log", action="store_true", help="logarithmic colour scale (zero wells shown as grey)")
    ap.add_argument("--no-values", action="store_true", help="don't print the numbers in each well")
    ap.add_argument("--cmap", default="inferno", help="matplotlib colormap (default: inferno)")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    args = ap.parse_args()

    if not args.show:
        matplotlib.use("Agg")

    column, scale, unit, scales_with_activity = QUANTITIES[args.quantity]
    values, rel_err = load(args.csv, column)
    values = values * scale
    if scales_with_activity:  # i.e. the dose rate
        values = values * args.activity * TIME_FACTORS[args.time_unit]
        unit = f"Gy/{args.time_unit}"
    values, unit = pick_unit(values, unit)

    title = {"dose_rate": "Dose rate", "dose_per_decay": "Dose per decay",
             "edep_per_decay": "Energy deposited per decay"}[args.quantity]
    if scales_with_activity:
        title += f" at {args.activity:g} MBq Sr-90"

    if not np.any(values > 0):
        print("Warning: every well has zero dose, so the heatmap will be empty. "
              "Check that the primary generator is really emitting Sr-90.", file=sys.stderr)

    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad("0.85")  # wells with zero dose when using a log scale

    if args.log and np.any(values > 0):
        positive = values[values > 0]
        norm = LogNorm(vmin=positive.min(), vmax=positive.max())
        shown = np.ma.masked_less_equal(values, 0)
    else:
        vmax = np.nanmax(values)
        norm = Normalize(vmin=0, vmax=vmax if vmax > 0 else 1)
        shown = values

    fig, ax = plt.subplots(figsize=(10, 7.5))
    im = ax.imshow(shown, cmap=cmap, norm=norm, aspect="equal")  # row 0 (A) at the top

    ax.set_xticks(range(COLS), labels=[str(c + 1) for c in range(COLS)])
    ax.set_yticks(range(ROWS), labels=list(ROW_NAMES))
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(-0.5, COLS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, ROWS, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0)

    if not args.no_values:
        for r in range(ROWS):
            for c in range(COLS):
                v = values[r, c]
                frac = norm(v) if (v > 0 or not args.log) else 0
                frac = 0 if np.ma.is_masked(frac) else float(frac)
                colour = "black" if frac > 0.6 else "white"
                ax.text(c, r, f"{v:.2e}" if v != 0 else "0", ha="center", va="center",
                        fontsize=8, color=colour)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(unit)
    ax.set_title(title, pad=30, fontsize=13)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Wrote {args.output}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
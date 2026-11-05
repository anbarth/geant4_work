#!/usr/bin/env python3
"""Dose versus distance from the centre of the source, from the film voxel results.

Reads voxel_dose.csv (columns ix,iy,x_mm,y_mm,area_mm2,edep_per_decay_MeV,dose_per_decay_Gy,
rel_error,dose_rate_Gy_per_h_at_1MBq). For every voxel the distance r from the chosen centre
(default: the film centre, which is the source axis when the source is not offset) is computed
from the voxel centre. Voxels at the same distance are then averaged together:

    default       voxels whose distances agree to within --tol (default 1e-6 mm) form one point.
                  On the square grid these are the symmetric copies of a voxel (up to 8 of them),
                  so each point is exactly one distinct radius.
    --bin-width W instead, voxels are grouped in radial bins of width W mm; the point sits at
                  the mean radius of the voxels in the bin.

The average is weighted by voxel area. The uncertainty of a point is the statistical errors of its
voxels added in quadrature (the voxels are treated as independent), divided by the number of voxels.
Optionally --spread uses the standard error from the scatter between the voxels instead.

Outputs a PNG and a CSV with columns  r_mm, dose, uncertainty, n_voxels  (the unit is in the header).

Examples
    python radial_profile.py voxel_dose.csv                       # exact-distance grouping
    python radial_profile.py voxel_dose.csv --bin-width 0.25      # 0.25 mm radial bins
    python radial_profile.py voxel_dose.csv -a 5 -t min --log     # Gy/min for 5 MBq, log axis
    python radial_profile.py voxel_dose.csv --center 1.0 -0.5     # centre not at the film centre
    python radial_profile.py voxel_dose.csv --source-radius 8.2   # mark the edge of the source disc
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib

QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", "Gy/decay", False),
}
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}


def si_prefix(vmax, unit):
    if not unit.startswith("Gy/") or vmax <= 0:
        return 1.0, unit
    if vmax < 1e-3:
        return 1e6, "µ" + unit
    if vmax < 1:
        return 1e3, "m" + unit
    return 1.0, unit


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/voxel_dose.csv")
    ap.add_argument("--center", nargs=2, type=float, default=[0.0, 0.0], metavar=("X", "Y"),
                    help="centre in mm, in the film coordinates of the CSV (default 0 0)")
    ap.add_argument("--bin-width", type=float, help="group voxels into radial bins of this width (mm)")
    ap.add_argument("--tol", type=float, default=1e-6, help="distance tolerance for exact grouping (mm)")
    ap.add_argument("--spread", action="store_true",
                    help="uncertainty from the scatter of the voxels (std. error of the mean) instead of "
                         "their statistical errors; needs 2 or more voxels per point")
    ap.add_argument("--rmax", type=float, help="ignore voxels further than this (mm)")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_rate")
    ap.add_argument("-a", "--activity", type=float, default=1.0, help="Sr-90 activity in MBq (dose_rate only)")
    ap.add_argument("-t", "--time-unit", choices=TIME_FACTORS, default="h")
    ap.add_argument("--log", action="store_true", help="logarithmic dose axis")
    ap.add_argument("--source-radius", type=float, help="draw a vertical line at the source radius (mm)")
    ap.add_argument("-o", "--output", default="radial_profile.png")
    ap.add_argument("--csv-out", help="output CSV (default: same name as the image, .csv)")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    column, unit, is_rate = QUANTITIES[args.quantity]
    with open(args.csv, newline="") as f:
        rd = csv.DictReader(f)
        for n in ("x_mm", "y_mm", "area_mm2", "rel_error", column):
            if n not in (rd.fieldnames or []):
                sys.exit(f"{args.csv}: column '{n}' not found. Found: {rd.fieldnames}")
        a = np.array([[float(r["x_mm"]), float(r["y_mm"]), float(r["area_mm2"]),
                       float(r[column]), float(r["rel_error"])] for r in rd])
    if len(a) == 0:
        sys.exit("No data rows.")
    x, y, area, val, rel = a.T
    r = np.hypot(x - args.center[0], y - args.center[1])
    if args.rmax is not None:
        keep = r <= args.rmax
        x, y, area, val, rel, r = (v[keep] for v in (x, y, area, val, rel, r))
    if r.size == 0:
        sys.exit("No voxels left after --rmax.")

    if is_rate:
        val = val * args.activity * TIME_FACTORS[args.time_unit]
        unit = f"Gy/{args.time_unit}"
    scale, unit = si_prefix(np.nanmax(val), unit)
    val = val * scale
    err = val * rel

    # group voxels
    order = np.argsort(r)
    r, area, val, err = r[order], area[order], val[order], err[order]
    if args.bin_width:
        group = np.floor(r / args.bin_width).astype(int)
    else:
        gaps = np.diff(r) > args.tol
        group = np.concatenate([[0], np.cumsum(gaps)])
    R, D, E, N = [], [], [], []
    for g in np.unique(group):
        m = group == g
        A = area[m].sum()
        n = int(m.sum())
        mean = float(np.sum(area[m] * val[m]) / A)
        if args.spread:
            if n < 2:
                continue
            w = area[m] / A
            unc = float(np.sqrt(np.sum(w * (val[m] - mean) ** 2) / (n - 1)))   # std. error of the mean
        else:
            unc = float(np.sqrt(np.sum((area[m] * err[m]) ** 2)) / A)
        R.append(float(np.average(r[m], weights=area[m]))); D.append(mean); E.append(unc); N.append(n)
    R, D, E, N = map(np.array, (R, D, E, N))

    out_csv = args.csv_out or (args.output.rsplit(".", 1)[0] + ".csv")
    ulab = unit.replace("/", "_per_").replace("µ", "u")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["r_mm", f"dose_{ulab}", f"uncertainty_{ulab}", "n_voxels"])
        for ri, di, ei, ni in zip(R, D, E, N):
            w.writerow([f"{ri:.6g}", f"{di:.6g}", f"{ei:.6g}", int(ni)])
    print(f"Wrote {out_csv} ({len(R)} points from {int(N.sum())} voxels)")

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    qname = "Dose rate" if is_rate else "Dose per decay"
    if is_rate:
        qname += f" ({args.activity:g} MBq)"
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.errorbar(R, D, yerr=E, fmt="o", ms=3, color="black", capsize=2, lw=0.8)
    if args.source_radius:
        ax.axvline(args.source_radius, color="tab:red", ls="--", lw=1, label="source radius")
        ax.legend()
    if args.log:
        ax.set_yscale("log")
    ax.set_xlabel(f"distance from ({args.center[0]:g}, {args.center[1]:g}) mm (mm)")
    ax.set_ylabel(f"{qname} [{unit}]")
    ax.set_title(f"Active layer: {qname} vs radius")
    ax.grid(alpha=0.3)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Wrote {args.output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

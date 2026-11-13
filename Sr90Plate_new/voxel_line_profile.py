#!/usr/bin/env python3
"""Dose along a straight line through the centre of a well, at any angle, across the whole plate.

Reads voxel_dose.csv (cell-layer voxel scoring, columns well,ix,iy,x_mm,y_mm,area_mm2,
edep_per_decay_MeV,dose_per_decay_Gy,rel_error,dose_rate_Gy_per_h_at_1MBq).

The line passes through the centre of the chosen well and runs from one edge of the plate to
the other. Angles are measured in degrees seen from above, anticlockwise from +x:
    0   = along the row (towards column 8)        90  = along the column (towards row A)
    180 = towards column 1                        270 = towards row F
Because the line is straight and has no width, a point's dose is that of the voxel it falls in.
Instead of a well you can give a point with --xy X Y (mm from the plate centre; +x towards column 8,
+y towards row A); the line then passes through that point and t is measured from it.
Positions are the signed distance t (mm) from the centre of the chosen well (or the point) along the line
(t increases in the direction of the angle). Where the line is not inside a well that has voxel
data (the wall, the gaps between wells, wells not in the CSV) there is no dose: the plot shows a
gap and the CSV leaves the dose columns empty.

Outputs: a PNG (dose vs position, every well the line crosses shaded and labelled, plus a small
map of the line on the plate) and a CSV with columns
    t_mm, x_mm, y_mm, well, dose, uncertainty        (x, y in plate coordinates, centred on the plate)

Examples
    python voxel_line_profile.py voxel_dose.csv --well C5 --angle 0
    python voxel_line_profile.py voxel_dose.csv --xy 6.5 6.5 --angle 30     # through a point, not a well centre
    python voxel_line_profile.py voxel_dose.csv --well C5 --angle 45 -a 5 -t min
    python voxel_line_profile.py voxel_dose.csv --well B3 --angle 90 --log --step 0.05
    python voxel_line_profile.py voxel_dose.csv --well C5 --angle 30 -o c5_30.png --csv-out c5_30.csv
"""
import argparse
import csv
import math
import sys

import numpy as np
import matplotlib

ROWS, COLS = 6, 8
ROW_NAMES = "ABCDEF"
QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", "Gy/decay", False),
}
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}


def well_centre(name, pitch):
    return ((int(name[1:]) - 1 - (COLS - 1) / 2.0) * pitch, ((ROWS - 1) / 2.0 - ROW_NAMES.index(name[0])) * pitch)


def parse_well(s):
    s = s.strip().upper()
    if len(s) < 2 or s[0] not in ROW_NAMES or not s[1:].isdigit() or not 1 <= int(s[1:]) <= COLS:
        sys.exit(f"Bad well name '{s}' - use e.g. C5 (rows A-F, columns 1-8).")
    return s


def si_prefix(values, unit):
    if not unit.startswith("Gy/"):
        return 1.0, unit
    vmax = np.nanmax(np.abs(values)) if np.any(np.isfinite(values)) else 0
    if vmax == 0:
        return 1.0, unit
    if vmax < 1e-3:
        return 1e6, "µ" + unit
    if vmax < 1:
        return 1e3, "m" + unit
    return 1.0, unit


def load_grids(path, column):
    """Return {well: (dose[iy, ix], rel_err[iy, ix])} and the voxel size / grid origin."""
    need = ("well", "ix", "iy", "x_mm", "y_mm", "rel_error", column)
    recs = []
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        miss = [n for n in need if n not in (rd.fieldnames or [])]
        if miss:
            sys.exit(f"{path}: missing column(s) {miss}. Found: {rd.fieldnames}")
        for r in rd:
            recs.append((r["well"], int(r["ix"]), int(r["iy"]), float(r["x_mm"]), float(r["y_mm"]),
                         float(r[column]), float(r["rel_error"])))
    if not recs:
        sys.exit(f"{path} has no data rows.")
    ix = np.array([r[1] for r in recs]); iy = np.array([r[2] for r in recs])
    x = np.array([r[3] for r in recs]); y = np.array([r[4] for r in recs])
    # voxel size and the position of the lower-left corner of voxel (0, 0), from the voxel centres
    ux = np.unique(np.round(x, 6)); d = np.diff(ux)
    vs = float(d[d > 1e-6].min())
    x0 = float(np.median(x - ix * vs)) - vs / 2
    y0 = float(np.median(y - iy * vs)) - vs / 2
    nx, ny = ix.max() + 1, iy.max() + 1
    grids = {}
    for w, i, j, _, _, v, e in recs:
        if w not in grids:
            grids[w] = (np.full((ny, nx), np.nan), np.full((ny, nx), np.nan))
        grids[w][0][j, i] = v
        grids[w][1][j, i] = e
    return grids, vs, x0, y0


def clip_line(c, d, hx, hy):
    """t range of the line c + t d inside the box |x|<=hx, |y|<=hy."""
    lo, hi = -math.inf, math.inf
    for ci, di, h in ((c[0], d[0], hx), (c[1], d[1], hy)):
        if abs(di) < 1e-12:
            if abs(ci) > h:
                sys.exit("The line misses the plate.")
            continue
        t1, t2 = (-h - ci) / di, (h - ci) / di
        lo, hi = max(lo, min(t1, t2)), min(hi, max(t1, t2))
    return lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/voxel_dose.csv", help="voxel CSV (default: build/voxel_dose.csv)")
    pt = ap.add_mutually_exclusive_group(required=True)
    pt.add_argument("-w", "--well", help="pass the line through the centre of this well, e.g. C5")
    pt.add_argument("--xy", nargs=2, type=float, metavar=("X", "Y"),
                    help="pass the line through this point, in mm from the plate centre "
                         "(+x towards column 8, +y towards row A)")
    ap.add_argument("--angle", type=float, required=True,
                    help="line angle in degrees, anticlockwise from +x (0 = along the row, 90 = along the column)")
    ap.add_argument("--step", type=float, default=0.5, help="sampling step along the line in mm (default 0.5)")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_rate")
    ap.add_argument("-a", "--activity", type=float, default=1000.0, help="Sr-90 activity in MBq (dose_rate only)")
    ap.add_argument("-t", "--time-unit", choices=TIME_FACTORS, default="h")
    ap.add_argument("--pitch", type=float, default=13.0, help="well pitch in mm (default 13.0)")
    ap.add_argument("--well-diameter", type=float, default=11.0, help="inner well diameter in mm (default 11.0)")
    ap.add_argument("--plate-x", type=float, default=127.76, help="plate length in x, mm")
    ap.add_argument("--plate-y", type=float, default=85.48, help="plate width in y, mm")
    ap.add_argument("--log", action="store_true", help="logarithmic dose axis")
    ap.add_argument("-o", "--output", help="output image (default line_<well>_<angle>deg.png)")
    ap.add_argument("--csv-out", help="output CSV (default: same name as the image, .csv)")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    args = ap.parse_args()

    well = parse_well(args.well) if args.well else None
    column, unit, is_rate = QUANTITIES[args.quantity]
    grids, vs, gx0, gy0 = load_grids(args.csv, column)
    if well is not None and well not in grids:
        sys.exit(f"Well {well} is not in {args.csv}. Wells present: {', '.join(sorted(grids))}")

    R = args.well_diameter / 2.0
    th = math.radians(args.angle)
    d = np.array([math.cos(th), math.sin(th)])
    c0 = np.array(well_centre(well, args.pitch)) if well else np.array(args.xy, dtype=float)
    if well:
        origin, tag = f"the centre of {well}", well
    else:
        origin, tag = f"({c0[0]:g}, {c0[1]:g}) mm", f"x{c0[0]:g}_y{c0[1]:g}"
    tmin, tmax = clip_line(c0, d, args.plate_x / 2, args.plate_y / 2)
    n = int(round((tmax - tmin) / args.step)) + 1
    t = np.linspace(tmin, tmax, n)
    px, py = c0[0] + t * d[0], c0[1] + t * d[1]

    # which well does each sample lie in, and which voxel
    dose = np.full(n, np.nan)
    err = np.full(n, np.nan)
    wname = np.full(n, "", dtype=object)
    colw = np.rint(px / args.pitch + (COLS - 1) / 2.0).astype(int)
    roww = np.rint((ROWS - 1) / 2.0 - py / args.pitch).astype(int)
    for k in range(n):
        if not (0 <= colw[k] < COLS and 0 <= roww[k] < ROWS):
            continue
        name = f"{ROW_NAMES[roww[k]]}{colw[k] + 1}"
        cx, cy = well_centre(name, args.pitch)
        lx, ly = px[k] - cx, py[k] - cy
        if math.hypot(lx, ly) > R:
            continue
        wname[k] = name
        if name not in grids:
            continue
        vals, rels = grids[name]
        i, j = int(math.floor((lx - gx0) / vs)), int(math.floor((ly - gy0) / vs))
        if 0 <= j < vals.shape[0] and 0 <= i < vals.shape[1] and np.isfinite(vals[j, i]):
            dose[k] = vals[j, i]
            err[k] = vals[j, i] * rels[j, i]

    if is_rate:
        f = args.activity * TIME_FACTORS[args.time_unit]
        dose, err = dose * f, err * f
        unit = f"Gy/{args.time_unit}"
    scale, unit = si_prefix(dose, unit)
    dose, err = dose * scale, err * scale
    if not np.any(np.isfinite(dose)):
        sys.exit("The line does not cross any voxel data - check the well and angle.")

    ang = f"{args.angle:g}"
    out_png = args.output or f"line_{tag}_{ang}deg.png"
    out_csv = args.csv_out or (out_png.rsplit(".", 1)[0] + ".csv")
    ulab = unit.replace("/", "_per_").replace("µ", "u")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_mm", "x_mm", "y_mm", "well", f"dose_{ulab}", f"uncertainty_{ulab}"])
        for k in range(n):
            ok = np.isfinite(dose[k])
            w.writerow([f"{t[k]:.4f}", f"{px[k]:.4f}", f"{py[k]:.4f}", wname[k],
                        f"{dose[k]:.6g}" if ok else "", f"{err[k]:.6g}" if ok else ""])
    print(f"Wrote {out_csv} ({n} points, t = {tmin:.1f} to {tmax:.1f} mm)")

    # wells the line crosses (chords through the circles), whether or not they are in the CSV
    chords = []
    for r in range(ROWS):
        for cidx in range(COLS):
            name = f"{ROW_NAMES[r]}{cidx + 1}"
            c = np.array(well_centre(name, args.pitch))
            t0 = float(np.dot(c - c0, d))
            perp = abs(float((c[0] - c0[0]) * d[1] - (c[1] - c0[1]) * d[0]))
            if perp < R:
                h = math.sqrt(R * R - perp * perp)
                chords.append((name, t0 - h, t0 + h))

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle

    fig = plt.figure(figsize=(12, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.2, 1], wspace=0.18)
    ax = fig.add_subplot(gs[0])
    for name, a, b in chords:
        target = well is not None and name == well
        ax.axvspan(a, b, color="0.8" if target else "0.92", zorder=0)
        ax.axvline(a, color="0.45", lw=0.8, ls="--", zorder=1)
        ax.axvline(b, color="0.45", lw=0.8, ls="--", zorder=1)
        ax.text(0.5 * (a + b), 1.01, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=8, fontweight="bold" if target else "normal")
    ax.plot(t, dose, color="tab:blue", lw=1.2, label="dose")
    ax.fill_between(t, dose - err, dose + err, color="tab:blue", alpha=0.25, lw=0, label="± uncertainty")
    ax.set_xlim(tmin, tmax)
    if args.log:
        ax.set_yscale("log")
    qname = "Cell-layer dose rate" if is_rate else "Cell-layer dose per decay"
    if is_rate:
        qname += f" ({args.activity:g} MBq)"
    ax.set_ylabel(f"{qname} [{unit}]")
    ax.set_xlabel(f"position along the line (mm from {origin})")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"Line through {well if well else origin} at {ang}°", pad=22)

    # map of the line on the plate
    am = fig.add_subplot(gs[1])
    am.add_patch(Rectangle((-args.plate_x / 2, -args.plate_y / 2), args.plate_x, args.plate_y,
                           fill=False, ec="0.3"))
    for r in range(ROWS):
        for cidx in range(COLS):
            name = f"{ROW_NAMES[r]}{cidx + 1}"
            cx, cy = well_centre(name, args.pitch)
            am.add_patch(Circle((cx, cy), R, fc="0.8" if (well and name == well) else "none", ec="0.5", lw=0.6))
    am.plot(px[[0, -1]], py[[0, -1]], color="tab:red", lw=1.5)
    am.plot(*c0, "o", color="tab:red", ms=3)
    am.annotate("", xy=(c0[0] + 12 * d[0], c0[1] + 12 * d[1]), xytext=(c0[0], c0[1]),
                arrowprops=dict(arrowstyle="->", color="tab:red"))
    am.set_xlim(-args.plate_x / 2 - 2, args.plate_x / 2 + 2)
    am.set_ylim(-args.plate_y / 2 - 2, args.plate_y / 2 + 2)
    am.set_aspect("equal")
    am.set_xticks([]); am.set_yticks([])
    am.set_title("top view (A at top)", fontsize=8)

    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"Wrote {out_png}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fit a dose profile through the cell layer using the voxel results (voxel_dose.csv).

fit_profile.py has one point per well. This script uses the voxels, so a profile has many
points per well. Pick a line, and every voxel in a strip around it contributes:

    --well C5            profile through ONE well along x (or --axis y); x is mm from the well centre
    --row C              profile along x through the centres of all wells in row C (plate coordinates)
    --col 5              profile along y through the centres of all wells in column 5

Voxels with the same position along the line are combined (area-weighted mean dose, errors added
in quadrature), so a strip several voxels wide gives one point per voxel column. Voxels partly
outside the well wall are included but weighted by their area; --min-area drops the small ones.

The models are the ones in fit_profile.py (same directory): `--list-models`, `-m <name>`, and
anything you register there is available here too.

Examples
    python fit_voxel_profile.py voxel_dose.csv --well C5                 # Gaussian across well C5
    python fit_voxel_profile.py voxel_dose.csv --well C5 --axis y -m lorentzian
    python fit_voxel_profile.py voxel_dose.csv --row C                   # whole row, voxel resolution
    python fit_voxel_profile.py voxel_dose.csv --row C --strip 2.0       # strip 2 mm wide (default: 2 voxels)
    python fit_voxel_profile.py voxel_dose.csv --col 5 --offset 1.5      # line shifted by 1.5 mm
    python fit_voxel_profile.py voxel_dose.csv --well C5 -a 5 -t min     # Gy/min for 5 MBq
    python fit_voxel_profile.py voxel_dose.csv --row C --csv-out         # also write x, y, uncertainty to a CSV
    python fit_voxel_profile.py voxel_dose.csv --row C --csv-out pts.csv -m none
    python fit_voxel_profile.py voxel_dose.csv --well C5 --min-area 0.1 --max-rel-error 0.3
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib
from scipy.optimize import curve_fit

import fit_profile as fp   # models, units helpers (must sit next to this script)

ROWS, COLS = fp.ROWS, fp.COLS
ROW_NAMES = fp.ROW_NAMES
QUANTITIES = {k: v for k, v in fp.QUANTITIES.items() if k != "edep_per_decay"}  # per-voxel energy is not intensive


def well_position(name, pitch):
    """Centre (x, y) in mm of a well such as 'C5' (+x towards column 8, +y towards row A)."""
    return (int(name[1:]) - 1 - (COLS - 1) / 2.0) * pitch, ((ROWS - 1) / 2.0 - ROW_NAMES.index(name[0])) * pitch


def parse_well(s):
    s = s.strip().upper()
    if len(s) < 2 or s[0] not in ROW_NAMES or not s[1:].isdigit() or not 1 <= int(s[1:]) <= COLS:
        sys.exit(f"Bad well name '{s}' - use e.g. C5 (rows A-F, columns 1-8).")
    return s


def load(path, column):
    need = ("x_mm", "y_mm", "area_mm2", "rel_error", column)
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        miss = [n for n in need if n not in (rd.fieldnames or [])]
        if miss:
            sys.exit(f"{path}: missing column(s) {miss}. Found: {rd.fieldnames}")
        has_well = "well" in rd.fieldnames
        rows = [(r["well"] if has_well else "", float(r["x_mm"]), float(r["y_mm"]), float(r["area_mm2"]),
                 float(r[column]), float(r["rel_error"])) for r in rd]
    if not rows:
        sys.exit(f"{path} has no data rows.")
    a = np.array([r[1:] for r in rows])
    return np.array([r[0] for r in rows]), a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]


def voxel_size(v):
    u = np.unique(np.round(v, 6))
    d = np.diff(u)
    return float(d[d > 1e-6].min()) if len(d) else 1.0


def combine(s, area, val, err):
    """Group voxels by position s along the line: area-weighted mean and its error."""
    key = np.round(s, 4)
    xs, ys, es = [], [], []
    for k in np.unique(key):
        m = key == k
        A = area[m].sum()
        xs.append(float(s[m].mean()))
        ys.append(float(np.sum(area[m] * val[m]) / A))
        es.append(float(np.sqrt(np.sum((area[m] * err[m]) ** 2)) / A))
    return np.array(xs), np.array(ys), np.array(es)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/voxel_dose.csv", help="voxel CSV (default: build/voxel_dose.csv)")
    line = ap.add_mutually_exclusive_group()
    line.add_argument("-w", "--well", help="one well, e.g. C5 (x = mm from the well centre)")
    line.add_argument("-r", "--row", help="row A-F: profile along x through all its wells")
    line.add_argument("-c", "--col", type=int, help="column 1-8: profile along y through all its wells")
    ap.add_argument("--axis", choices=["x", "y"], default="x", help="direction of the profile for --well (default x)")
    ap.add_argument("--strip", type=float, help="full width of the strip around the line in mm (default 2 voxels)")
    ap.add_argument("--offset", type=float, default=0.0, help="shift the line sideways by this many mm")
    ap.add_argument("-m", "--model", default="gaussian", help="model name (see --list-models), default gaussian; 'none' = plot only, no fit")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_rate")
    ap.add_argument("-a", "--activity", type=float, default=1000.0, help="Sr-90 activity in MBq (dose_rate only)")
    ap.add_argument("-t", "--time-unit", choices=fp.TIME_FACTORS, default="h")
    ap.add_argument("--pitch", type=float, default=13.0, help="well pitch in mm (default 13.0)")
    ap.add_argument("--well-diameter", type=float, default=11.0, help="inner well diameter in mm (for the shading)")
    ap.add_argument("--min-area", type=float, default=0.0, help="ignore voxels with less than this area in mm2")
    ap.add_argument("--max-rel-error", type=float, default=None, help="ignore points with a larger relative error")
    ap.add_argument("--no-weights", action="store_true", help="plain least squares, ignore the errors")
    ap.add_argument("-o", "--output", help="output image")
    ap.add_argument("--csv-out", nargs="?", const="", metavar="PATH",
                    help="also write the profile points (x, y, uncertainty) to a CSV; "
                         "PATH optional (default voxprofile_<line>_<quantity>.csv)")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    ap.add_argument("--list-models", action="store_true")
    args = ap.parse_args()

    if args.list_models:
        for m in fp.MODELS.values():
            print(f"{m.name:16s} {m.doc}\n{'':16s} parameters: {', '.join(m.params)}")
        return
    if args.well is None and args.row is None and args.col is None:
        ap.error("give --well C5, --row A-F or --col 1-8")
    if args.model not in fp.MODELS:
        sys.exit(f"Unknown model '{args.model}'. Available: {', '.join(fp.MODELS)}")
    model = fp.MODELS[args.model]

    column, unit, is_rate = QUANTITIES[args.quantity]
    well, vx, vy, area, val, rel = load(args.csv, column)
    vs = voxel_size(vx)
    strip = args.strip if args.strip else 2 * vs
    half = strip / 2.0
    R = args.well_diameter / 2.0

    # absolute voxel positions and the wells that are present
    centres, sel = {}, np.zeros(len(vx), bool)
    names = sorted(set(well)) if well[0] != "" else []
    if args.well:
        w = parse_well(args.well)
        if w not in names and names:
            sys.exit(f"Well {w} is not in {args.csv}. Wells present: {', '.join(names)}")
        along, across = (vx, vy) if args.axis == "x" else (vy, vx)
        sel = (np.abs(across - args.offset) <= half + 1e-9) & ((well == w) | (well == ""))
        s = along
        label = f"well {w}, along {args.axis} (strip {strip:g} mm)"
        axis_label = f"{args.axis} from the well centre (mm)"
        shade = [0.0]
    else:
        if not names:
            sys.exit("--row/--col need a voxel CSV with a 'well' column.")
        ax_i = 0 if args.row is not None else 1
        if args.row is not None:
            r = args.row.upper()
            if len(r) != 1 or r not in ROW_NAMES:
                sys.exit("--row must be one of A-F")
        elif not 1 <= args.col <= COLS:
            sys.exit("--col must be 1-8")
        sabs = np.full(len(vx), np.nan)
        across_abs = np.full(len(vx), np.nan)
        shade = []
        for w in names:
            cx, cy = well_position(w, args.pitch)
            inline = (w[0] == args.row.upper()) if args.row is not None else (int(w[1:]) == args.col)
            if not inline:
                continue
            m = well == w
            sabs[m] = (vx[m] + cx) if args.row is not None else (vy[m] + cy)
            across_abs[m] = (vy[m] - 0.0) if args.row is not None else (vx[m] - 0.0)  # local = relative to centre
            shade.append(cx if args.row is not None else cy)
        sel = np.isfinite(sabs) & (np.abs(np.nan_to_num(across_abs) - args.offset) <= half + 1e-9)
        s = sabs
        if args.row is not None:
            label = f"row {args.row.upper()} (strip {strip:g} mm)"
            axis_label = "x position (mm from plate centre)"
        else:
            label = f"column {args.col} (strip {strip:g} mm)"
            axis_label = "y position (mm from plate centre, +y toward row A)"

    sel &= area >= max(args.min_area, 1e-12)
    if not np.any(sel):
        sys.exit("No voxels selected - check the well/row/column, --strip, --offset and --min-area.")

    v = val[sel].copy()
    if is_rate:
        v *= args.activity * fp.TIME_FACTORS[args.time_unit]
        unit = f"Gy/{args.time_unit}"
    scale, unit = fp.si_prefix(v, unit)
    v *= scale
    e = v * rel[sel]
    x, y, sigma = combine(s[sel], area[sel], v, e)
    order = np.argsort(x)
    x, y, sigma = x[order], y[order], sigma[order]

    if args.csv_out is not None:
        tagc = (args.well or args.row or f"col{args.col}").upper() if (args.well or args.row) else f"col{args.col}"
        path = args.csv_out or f"voxprofile_{tagc}_{args.quantity}.csv"
        ulab = unit.replace("/", "_per_").replace("µ", "u")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["x_mm", f"y_{ulab}", f"uncertainty_{ulab}"])
            for xi, yi, si in zip(x, y, sigma):
                w.writerow([f"{xi:.6g}", f"{yi:.6g}", f"{si:.6g}"])
        print(f"Wrote {path} ({len(x)} points)")

    qname = {"dose_rate": "Cell-layer dose rate", "dose_per_decay": "Cell-layer dose per decay"}[args.quantity]
    if is_rate:
        qname += f" ({args.activity:g} MBq)"

    if model.name == "none":   # no fit: just show the dose distribution
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4.5))
        for c in shade:
            ax.axvspan(c - R, c + R, color="0.9", zorder=0)
        ax.errorbar(x, y, yerr=sigma, fmt="o", ms=3.5, color="black", capsize=2)
        ax.set_ylabel(f"{qname} [{unit}]")
        ax.set_xlabel(axis_label + "   (grey = inside the well)")
        ax.set_title(f"{qname}: {label}", fontsize=12)
        ax.grid(alpha=0.3)

        ax.set_ylim(0,350)

        tag = (args.well or args.row or f"col{args.col}").upper() if (args.well or args.row) else f"col{args.col}"
        out = args.output or f"voxprofile_{tag}.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"Wrote {out} (no fit, {len(x)} points)")

        if args.show:
            plt.show()
        return

    # drop empty / noisy points
    good = y > 0
    if args.max_rel_error is not None:
        good &= (sigma <= args.max_rel_error * np.where(y > 0, y, np.inf))
    weighted = not args.no_weights
    if weighted:
        good &= sigma > 0
    dropped = int((~good).sum())
    if dropped:
        print(f"Note: {dropped} of {len(x)} points ignored (zero dose/error or above --max-rel-error).")
    xf_pts, yf_pts, sf_pts = x[good], y[good], sigma[good]
    n, k = len(xf_pts), len(model.params)
    if n == 0:
        sys.exit("No usable points left to fit.")
    if n <= k:
        print(f"Warning: {n} points for {k} parameters - no degrees of freedom left.", file=sys.stderr)


    kwargs = {"p0": model.guess(xf_pts, yf_pts), "maxfev": 20000}
    if model.bounds is not None:
        kwargs["bounds"] = model.bounds(xf_pts, yf_pts)
    if weighted:
        kwargs.update(sigma=sf_pts, absolute_sigma=True)
    try:
        popt, pcov = curve_fit(model.func, xf_pts, yf_pts, **kwargs)
    except (RuntimeError, ValueError) as exc:
        sys.exit(f"Fit failed: {exc}\nTry a different model or check the data.")
    perr = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.full(k, np.nan)

    # replace fit parameters with initial guess
    popt = model.guess(xf_pts, yf_pts)

    resid = yf_pts - model.func(xf_pts, *popt)
    ss_tot = float(np.sum((yf_pts - yf_pts.mean()) ** 2))
    r2 = 1 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    ndf = n - k

    print(f"\n{qname}, {label}, model '{model.name}': {model.doc}")
    print(f"Voxel size {vs:g} mm, points fitted: {n}, weights: {'1/sigma (absolute)' if weighted else 'none'}")
    for nm, pv, pe in zip(model.params, popt, perr):
        print(f"  {nm:>8s} = {pv:.5g} +/- {pe:.2g}")
    if model.derived is not None:
        for nm, (pv, pe) in model.derived(popt, perr).items():
            print(f"  {nm:>8s} = {pv:.5g} +/- {pe:.2g}   (derived)")
    print(f"  R^2 = {r2:.5f}")
    if weighted and ndf > 0:
        chi2 = float(np.sum((resid / sf_pts) ** 2))
        print(f"  chi2/ndf = {chi2:.3g}/{ndf} = {chi2 / ndf:.3g}")

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pad = 0.5 * args.pitch
    xf = np.linspace(x.min() - 0.1 * pad, x.max() + 0.1 * pad, 600)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    for a in (a1, a2):
        for c in shade:
            a.axvspan(c - R, c + R, color="0.9", zorder=0)
    a1.errorbar(x[good], y[good], yerr=sigma[good], fmt="o", ms=3.5, color="black", capsize=2, label="voxels")
    if dropped:
        a1.plot(x[~good], y[~good], "x", color="0.6", label="ignored")

    a1.plot(xf, model.func(xf, *popt), color="tab:red", label=f"{model.name} fit")
    a1.set_ylabel(f"{qname} [{unit}]")
    a1.set_title(f"{qname}: {label}", fontsize=12)
    txt = "\n".join(f"{nm} = {pv:.4g} ± {pe:.2g}" for nm, pv, pe in zip(model.params, popt, perr)) + f"\nR² = {r2:.4f}"
    a1.text(0.02, 0.97, txt, transform=a1.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    a1.legend(loc="upper right")
    a1.grid(alpha=0.3)
    a2.axhline(0, color="tab:red", lw=1)
    a2.errorbar(xf_pts, resid, yerr=sf_pts, fmt="o", ms=3.5, color="black", capsize=2)
    a2.set_ylabel("data − fit")
    a2.set_xlabel(axis_label + "   (grey = inside the well)")
    a2.grid(alpha=0.3)

    tag = (args.well or (args.row or f"col{args.col}")).upper() if args.well or args.row else f"col{args.col}"
    out = args.output or f"voxfit_{tag}_{model.name}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Wrote {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

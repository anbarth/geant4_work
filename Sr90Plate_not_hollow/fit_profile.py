#!/usr/bin/env python3
"""Fit the dose profile along one row or one column of the well plate.

Reads well_dose.csv (written by the Sr90Plate simulation), takes either a row
(A-F, profile along x) or a column (1-8, profile along y), and fits it with a
model function. Positions are in mm from the plate centre, which is also the
source axis when the applicator is not offset (fSourceX/fSourceY = 0).

Examples
    python fit_profile.py --row C                          # Gaussian fit along row C
    python fit_profile.py --col 4 -m gaussian_offset       # column 4, Gaussian + constant
    python fit_profile.py --row D -q dose_per_decay        # fit dose per decay instead
    python fit_profile.py --row C -t min -a 5              # dose rate in Gy/min for 5 MBq
    python fit_profile.py --list-models

Adding a model: see the "MODELS" section below. Write the function, then call
register(...) with parameter names and a function for the initial guess. It then
shows up in --list-models and works with -m <name> straight away.
"""

"""
Human notes on this code:
- Since we need real X and Y positions, make sure you know the pitch! (can be supplied as an extra argument)
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Sequence

import numpy as np
import matplotlib
from scipy.optimize import curve_fit

ROWS, COLS = 6, 8
ROW_NAMES = "ABCDEF"

# quantity -> (csv column, unit label, is it a dose rate?)
QUANTITIES = {
    "dose_rate": ("dose_rate_Gy_per_h_at_1MBq", "Gy/h", True),
    "dose_per_decay": ("dose_per_decay_Gy", "Gy/decay", False),
    "edep_per_decay": ("edep_per_decay_MeV", "MeV/decay", False),
}
TIME_FACTORS = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}  # per-hour -> per-unit


# =============================================================================
# MODELS  -  edit / extend this section to change what gets fitted
# =============================================================================
@dataclass
class Model:
    name: str
    func: Callable                      # func(x, *params) -> y
    params: Sequence[str]               # parameter names, in the order func takes them
    guess: Callable                     # guess(x, y) -> list of initial values
    bounds: Optional[Callable] = None   # bounds(x, y) -> (lower_list, upper_list); None = unbounded
    derived: Optional[Callable] = None  # derived(popt, perr) -> {label: (value, error)}
    doc: str = ""


MODELS: Dict[str, Model] = {}


def register(model: Model):
    MODELS[model.name] = model


def _moments(x, y):
    """Centre and RMS width of the profile, treating y as a weight."""
    w = np.clip(y, 0, None)
    if w.sum() <= 0:
        return float(np.mean(x)), float(np.ptp(x)) / 4 or 1.0
    x0 = float(np.sum(w * x) / w.sum())
    s = float(np.sqrt(np.sum(w * (x - x0) ** 2) / w.sum()))
    return x0, max(s, float(np.ptp(x)) / 20)


def _span(x):
    return float(np.ptp(x)) or 1.0


# ---- Gaussian: A * exp(-(x-x0)^2 / (2 sigma^2)) ------------------------------
def gaussian(x, A, x0, sigma):
    return A * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def _gauss_guess(x, y):
    x0, s = _moments(x, y)
    return [float(np.max(y)), x0, s]


FWHM_GAUSS = 2.0 * np.sqrt(2.0 * np.log(2.0))  # 2.3548...

register(Model(
    "gaussian", gaussian, ["A", "x0", "sigma"], _gauss_guess,
    bounds=lambda x, y: ([0, x.min() - _span(x), 1e-6], [np.inf, x.max() + _span(x), 10 * _span(x)]),
    derived=lambda p, e: {"FWHM": (FWHM_GAUSS * p[2], FWHM_GAUSS * e[2])},
    doc="A * exp(-(x-x0)^2 / (2 sigma^2))",
))


# ---- Gaussian on a constant background ---------------------------------------
def gaussian_offset(x, A, x0, sigma, c):
    return c + A * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


register(Model(
    "gaussian_offset", gaussian_offset, ["A", "x0", "sigma", "c"],
    lambda x, y: _gauss_guess(x, y)[:3] + [float(np.min(y))],
    bounds=lambda x, y: ([0, x.min() - _span(x), 1e-6, -np.inf],
                         [np.inf, x.max() + _span(x), 10 * _span(x), np.inf]),
    derived=lambda p, e: {"FWHM": (FWHM_GAUSS * p[2], FWHM_GAUSS * e[2])},
    doc="c + A * exp(-(x-x0)^2 / (2 sigma^2))",
))


# ---- Lorentzian (heavier tails than a Gaussian) -------------------------------
def lorentzian(x, A, x0, gamma):
    return A / (1.0 + ((x - x0) / gamma) ** 2)


register(Model(
    "lorentzian", lorentzian, ["A", "x0", "gamma"],
    lambda x, y: [float(np.max(y)), _moments(x, y)[0], _moments(x, y)[1]],
    bounds=lambda x, y: ([0, x.min() - _span(x), 1e-6], [np.inf, x.max() + _span(x), 10 * _span(x)]),
    derived=lambda p, e: {"FWHM": (2 * p[2], 2 * e[2])},
    doc="A / (1 + ((x-x0)/gamma)^2)      (gamma = half width at half maximum)",
))


# ---- Exponential fall-off from the centre ------------------------------------
def exp_tail(x, A, x0, L):
    return A * np.exp(-np.abs(x - x0) / L)


register(Model(
    "exponential", exp_tail, ["A", "x0", "L"],
    lambda x, y: [float(np.max(y)), _moments(x, y)[0], _moments(x, y)[1]],
    bounds=lambda x, y: ([0, x.min() - _span(x), 1e-6], [np.inf, x.max() + _span(x), 10 * _span(x)]),
    doc="A * exp(-|x-x0| / L)",
))


# ---- Geometric: point source at height h above the plane ---------------------
# Flux on a flat plane from a point source, no absorption in between:
# cos(theta)/r^2 -> (1 + ((x-x0)/h)^2)^(-3/2). Ignores beta attenuation in air.
def point_source(x, A, x0, h):
    return A / (1.0 + ((x - x0) / h) ** 2) ** 1.5


register(Model(
    "point_source", point_source, ["A", "x0", "h"],
    lambda x, y: [float(np.max(y)), _moments(x, y)[0], _moments(x, y)[1] * 2],
    bounds=lambda x, y: ([0, x.min() - _span(x), 1e-3], [np.inf, x.max() + _span(x), 100 * _span(x)]),
    doc="A / (1 + ((x-x0)/h)^2)^(3/2)    (geometry only, h = effective source height)",
))
# To add your own: define func(x, *params), then
#   register(Model("name", func, ["p1", "p2", ...], guess=lambda x, y: [...]))
# =============================================================================


def read_profile(path, column, row=None, col=None, pitch=13.0):
    """Return x (mm from plate centre), values, relative errors for the chosen line."""
    vals, errs, pos = [], [], []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for name in (column, "rel_error", "row", "col"):
            if name not in (reader.fieldnames or []):
                sys.exit(f"Column '{name}' not found in {path}. Found: {reader.fieldnames}")
        for rec in reader:
            r, c = int(rec["row"]), int(rec["col"])
            if row is not None and r == row:
                pos.append((c - (COLS - 1) / 2.0) * pitch)          # x, column 1 at -x
            elif col is not None and c == col:
                pos.append(((ROWS - 1) / 2.0 - r) * pitch)          # y, row A at +y
            else:
                continue
            vals.append(float(rec[column]))
            errs.append(float(rec["rel_error"]))
    if not vals:
        sys.exit("No wells selected - check --row / --col.")
    order = np.argsort(pos)
    return (np.array(pos)[order], np.array(vals)[order], np.array(errs)[order])


def si_prefix(values, unit):
    """Rescale Gy/<time> values to mGy or µGy when that is easier to read."""
    if not unit.startswith("Gy/"):
        return 1.0, unit
    vmax = np.nanmax(np.abs(values))
    if vmax == 0:
        return 1.0, unit
    if vmax < 1e-3:
        return 1e6, "µ" + unit
    if vmax < 1:
        return 1e3, "m" + unit
    return 1.0, unit


def main():
    # parse arguments
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/well_dose.csv", help="input CSV (default: build/well_dose.csv)")
    line = ap.add_mutually_exclusive_group()
    line.add_argument("-r", "--row", help="row to fit, A-F (profile along x)")
    line.add_argument("-c", "--col", type=int, help="column to fit, 1-8 (profile along y)")
    ap.add_argument("-m", "--model", default="gaussian", help="model name (see --list-models), default gaussian")
    ap.add_argument("-q", "--quantity", choices=QUANTITIES, default="dose_per_decay",
                    help="what to fit (default: dose_per_decay)")
    ap.add_argument("-a", "--activity", type=float, default=1.0, help="Sr-90 activity in MBq (dose_rate only)")
    ap.add_argument("-t", "--time-unit", choices=TIME_FACTORS, default="h", help="dose-rate time unit (h, min, s)")
    ap.add_argument("--pitch", type=float, default=13.0, help="well pitch in mm (default 13.0, as in the geometry)")
    ap.add_argument("--no-weights", action="store_true", help="ignore rel_error; plain least squares")
    ap.add_argument("-o", "--output", help="output image (default fit_<row|col>_<model>.png)")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    ap.add_argument("--list-models", action="store_true", help="list available models and exit")
    args = ap.parse_args()

    # handle model argument
    if args.list_models:
        for m in MODELS.values():
            print(f"{m.name:16s} {m.doc}\n{'':16s} parameters: {', '.join(m.params)}")
        return
    if args.row is None and args.col is None:
        ap.error("give --row A-F or --col 1-8")
    if args.model not in MODELS:
        sys.exit(f"Unknown model '{args.model}'. Available: {', '.join(MODELS)}")
    model = MODELS[args.model]

    # handle row/column argument
    row_idx = col_idx = None
    if args.row is not None:
        if args.row.upper() not in ROW_NAMES or len(args.row) != 1:
            sys.exit("--row must be one of A-F")
        row_idx = ROW_NAMES.index(args.row.upper())
        label, axis_label = f"row {args.row.upper()}", "x position along row (mm from plate centre)"
    else:
        if not 1 <= args.col <= COLS:
            sys.exit("--col must be 1-8")
        col_idx = args.col - 1
        label = f"column {args.col}"
        axis_label = "y position along column (mm from plate centre, +y toward row A)"

    # read in the desired quantity (eg dose/event, dose rate) from the csv
    column, unit, is_rate = QUANTITIES[args.quantity]
    x, y, rel = read_profile(args.csv, column, row_idx, col_idx, args.pitch)
    if is_rate:
        y = y * args.activity * TIME_FACTORS[args.time_unit]
        unit = f"Gy/{args.time_unit}"
    scale, unit = si_prefix(y, unit)
    y = y * scale
    sigma = y * rel  # statistical error of each point, in the same units

    if not np.any(y > 0):
        sys.exit("Every well on this line has zero dose - nothing to fit. "
                 "Check that the primary generator is emitting Sr-90.")

    weighted = not args.no_weights
    if weighted and (np.any(y <= 0) or np.any(sigma <= 0)):
        print("Note: some wells have zero dose or zero error, so the fit is unweighted.")
        weighted = False

    n, k = len(x), len(model.params)
    if n <= k:
        print(f"Warning: {n} points for {k} parameters - no degrees of freedom left.", file=sys.stderr)

    p0 = model.guess(x, y)
    kwargs = {"p0": p0, "maxfev": 20000}
    if model.bounds is not None:
        kwargs["bounds"] = model.bounds(x, y)
    if weighted:
        kwargs.update(sigma=sigma, absolute_sigma=True)
    try:
        popt, pcov = curve_fit(model.func, x, y, **kwargs)
    except (RuntimeError, ValueError) as exc:
        sys.exit(f"Fit failed: {exc}\nTry a different model or check the data.")
    perr = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.full(k, np.nan)

    resid = y - model.func(x, *popt)
    ss_res, ss_tot = float(np.sum(resid ** 2)), float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    ndf = n - k

    qname = {"dose_rate": "Dose rate", "dose_per_decay": "Dose per decay",
             "edep_per_decay": "Energy per decay"}[args.quantity]
    if is_rate:
        qname += f" ({args.activity:g} MBq)"
    print(f"\n{qname}, {label}, model '{model.name}': {model.doc}")
    print(f"Weights: {'1/rel_error (absolute sigma)' if weighted else 'none'}   points: {n}")
    for name, v, e in zip(model.params, popt, perr):
        print(f"  {name:>8s} = {v:.5g} +/- {e:.2g}")
    if model.derived is not None:
        for name, (v, e) in model.derived(popt, perr).items():
            print(f"  {name:>8s} = {v:.5g} +/- {e:.2g}   (derived)")
    print(f"  R^2 = {r2:.5f}")
    if weighted and ndf > 0:
        chi2 = float(np.sum((resid / sigma) ** 2))
        print(f"  chi2/ndf = {chi2:.3g}/{ndf} = {chi2 / ndf:.3g}")

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xf = np.linspace(x.min() - 0.5 * args.pitch, x.max() + 0.5 * args.pitch, 400)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    ax.errorbar(x, y, yerr=sigma, fmt="o", color="black", capsize=3, label="simulation")
    ax.plot(xf, model.func(xf, *popt), color="tab:red", label=f"{model.name} fit")
    ax.set_ylabel(f"{qname} [{unit}]")
    ax.set_title(f"{qname}: {label}", fontsize=12)
    txt = "\n".join(f"{nm} = {v:.4g} ± {e:.2g}" for nm, v, e in zip(model.params, popt, perr))
    txt += f"\nR² = {r2:.4f}"
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    axr.axhline(0, color="tab:red", lw=1)
    axr.errorbar(x, resid, yerr=sigma, fmt="o", color="black", capsize=3)
    axr.set_ylabel("data − fit")
    axr.set_xlabel(axis_label)
    axr.grid(alpha=0.3)

    out = args.output or f"fit_{label.replace(' ', '')}_{model.name}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Wrote {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
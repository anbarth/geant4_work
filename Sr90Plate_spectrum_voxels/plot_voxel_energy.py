#!/usr/bin/env python3
"""Map the mean electron energy per cell-layer voxel (voxel_energy.csv from the simulation).

Examples
  python3 plot_voxel_energy.py                          # whole plate, mean energy
  python3 plot_voxel_energy.py --well C5                # one well
  python3 plot_voxel_energy.py -q error                 # statistical error map
  python3 plot_voxel_energy.py -q fluence --log         # electron fluence per decay
  python3 plot_voxel_energy.py --min-fluence 1e-6       # hide voxels with too few electrons
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

PITCH = 13.0  # mm, well pitch (must match DetectorConstruction)
ROWS, COLS = "ABCDEF", 8


def well_centre(name):
    r, c = ROWS.index(name[0].upper()), int(name[1:]) - 1
    return (c - (COLS - 1) / 2) * PITCH, ((len(ROWS) - 1) / 2 - r) * PITCH   # A is at +y


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="voxel_energy.csv")
    ap.add_argument("-o", "--output", default=None, help="output image (default: voxel_energy_<quantity>.png)")
    ap.add_argument("-q", "--quantity", choices=["mean_energy", "error", "rel_error", "fluence"],
                    default="mean_energy", help="what to plot (default: mean_energy)")
    ap.add_argument("--well", default=None, help="plot a single well, e.g. C5 (default: all wells in the file)")
    ap.add_argument("--vmin", type=float, default=None)
    ap.add_argument("--vmax", type=float, default=None)
    ap.add_argument("--log", action="store_true", help="logarithmic colour scale")
    ap.add_argument("--min-fluence", type=float, default=0.0,
                    help="blank voxels whose fluence per decay [1/mm2] is below this (their mean energy is noisy)")
    ap.add_argument("--cmap", default="viridis")
    args = ap.parse_args()

    d = np.genfromtxt(args.csv, delimiter=",", names=True, dtype=None, encoding=None)
    wells = np.array(d["well"])
    if args.well:
        keep = wells == args.well.upper()
        if not keep.any():
            raise SystemExit(f"Well {args.well} is not in {args.csv}")
        d, wells = d[keep], wells[keep]

    cx = np.array([well_centre(w)[0] for w in wells])
    cy = np.array([well_centre(w)[1] for w in wells])
    x, y = cx + d["x_mm"], cy + d["y_mm"]
    # voxel size from the grid spacing inside one well
    xs = np.unique(np.round(d["x_mm"][wells == wells[0]], 6))
    size = np.min(np.diff(xs)) if len(xs) > 1 else 0.5

    q = args.quantity
    if q == "mean_energy":
        val, label = d["mean_energy_keV"].astype(float), "Mean electron energy [keV]"
    elif q == "error":
        val, label = d["error_keV"].astype(float), "Statistical error of mean energy [keV]"
    elif q == "rel_error":
        with np.errstate(divide="ignore", invalid="ignore"):
            val = d["error_keV"] / d["mean_energy_keV"]
        label = "Relative statistical error"
    else:
        val, label = d["fluence_per_decay_per_mm2"].astype(float), "Electron fluence [per decay per mm$^2$]"

    bad = (d["track_length_per_decay_mm"] <= 0) | (d["fluence_per_decay_per_mm2"] < args.min_fluence)
    val = np.where(bad, np.nan, val)
    if args.log:
        val = np.where(val > 0, val, np.nan)

    h = size / 2
    quads = np.stack([np.stack([x - h, y - h], 1), np.stack([x + h, y - h], 1),
                      np.stack([x + h, y + h], 1), np.stack([x - h, y + h], 1)], 1)
    ok = np.isfinite(val)

    fig, ax = plt.subplots(figsize=(11, 7.5) if not args.well else (7, 6))
    norm = matplotlib.colors.LogNorm(args.vmin, args.vmax) if args.log else \
        matplotlib.colors.Normalize(args.vmin, args.vmax)
    pc = PolyCollection(quads[ok], array=val[ok], cmap=args.cmap, norm=norm, edgecolors="none")
    ax.add_collection(pc)
    if args.vmin is None and args.vmax is None and ok.any():
        pc.set_clim(np.nanmin(val[ok]), np.nanmax(val[ok]))
    ax.set_xlim(x.min() - size, x.max() + size)
    ax.set_ylim(y.min() - size, y.max() + size)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    fig.colorbar(pc, ax=ax, label=label, shrink=0.85)
    ax.set_title(label + (f" - well {args.well.upper()}" if args.well else ""))

    out = args.output or f"voxel_energy_{q}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Wrote {out}  ({ok.sum()} voxels plotted, {(~ok).sum()} blank)")


if __name__ == "__main__":
    main()

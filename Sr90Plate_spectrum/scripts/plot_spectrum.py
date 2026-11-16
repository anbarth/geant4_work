#!/usr/bin/env python3
"""Plot the electron (beta) energy spectrum in the spectrum plane next to the steel filter, from spectrum.csv.

spectrum.csv is written by the simulation: the kinetic energy of every electron that crosses the
spectrum plane moving away from the source, per simulated Sr-90 decay (each decay includes its Y-90 daughter).
Columns: e_low_MeV, e_high_MeV, e_mid_MeV, counts_per_decay, rel_error, counts_per_decay_per_MeV.

Examples
    python plot_spectrum.py                          # build/spectrum.csv -> spectrum.png (counts / decay / MeV)
    python plot_spectrum.py spectrum.csv --log       # logarithmic y axis
    python plot_spectrum.py spectrum.csv --normalise # area under the curve = 1
    python plot_spectrum.py spectrum.csv --per-bin   # counts per decay per bin instead of per MeV
    python plot_spectrum.py spectrum.csv --emax 2.5 -a 1000 # electrons per second for 1000 MBq (1 GBq)
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?", default="build/spectrum.csv")
    ap.add_argument("--log", action="store_true", help="logarithmic y axis")
    ap.add_argument("--normalise", action="store_true", help="normalise to unit area")
    ap.add_argument("--per-bin", action="store_true", help="plot counts per decay per bin (default: per MeV)")
    ap.add_argument("-a", "--activity", type=float, help="Sr-90 activity in MBq: plot electrons per second per MeV")
    ap.add_argument("--emax", type=float, help="upper limit of the energy axis in MeV")
    ap.add_argument("-o", "--output", default="spectrum.png")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    with open(args.csv, newline="") as f:
        rd = csv.DictReader(f)
        need = ("e_low_MeV", "e_high_MeV", "e_mid_MeV", "counts_per_decay", "rel_error")
        miss = [n for n in need if n not in (rd.fieldnames or [])]
        if miss:
            sys.exit(f"{args.csv}: missing column(s) {miss}. Found: {rd.fieldnames}")
        a = np.array([[float(r[n]) for n in need] for r in rd])
    if len(a) == 0:
        sys.exit("No data rows.")
    lo, hi, mid, cnt, rel = a.T
    w = hi - lo
    err = cnt * rel

    total = cnt.sum()
    mean_e = float(np.sum(cnt * mid) / total) if total > 0 else float("nan")
    cum = np.cumsum(cnt) / total if total > 0 else cnt
    median = float(np.interp(0.5, cum, hi)) if total > 0 else float("nan")
    print(f"Electrons per decay (in range): {total:.4g}")
    print(f"Mean energy {mean_e * 1000:.1f} keV, median {median * 1000:.1f} keV")

    y, e = cnt.copy(), err.copy()
    ylabel = "electrons per decay"
    if not args.per_bin:
        y, e = y / w, e / w
        ylabel += " per MeV"
    else:
        ylabel += f" per {w[0] * 1000:g} keV bin"
    if args.activity:
        f = args.activity * 1e6
        y, e = y * f, e * f
        ylabel = ylabel.replace("per decay", "per second")
        ylabel += f" ({args.activity:g} MBq)"
    if args.normalise:
        norm = np.sum(y * (w if not args.per_bin else 1.0))
        y, e = y / norm, e / norm
        ylabel = "probability density" if not args.per_bin else "fraction per bin"
        ylabel += " (area = 1)" if not args.per_bin else ""

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.step(np.append(lo, hi[-1]), np.append(y, y[-1]), where="post", color="tab:blue", lw=1.2)
    ax.fill_between(mid, y - e, y + e, color="tab:blue", alpha=0.25, lw=0, step=None)
    if args.log:
        ax.set_yscale("log")
        pos = y[y > 0]
        if pos.size:
            ax.set_ylim(pos.min() / 2, y.max() * 2)
    ax.set_xlim(0, args.emax if args.emax else hi[-1])
    ax.set_xlabel("electron kinetic energy (MeV)")
    ax.set_ylabel(ylabel)
    ax.set_title("Electron spectrum at the spectrum plane")
    ax.text(0.97, 0.95, f"{total:.3g} e⁻/decay\nmean {mean_e * 1000:.0f} keV\nmedian {median * 1000:.0f} keV",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.grid(alpha=0.3)
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    print(f"Wrote {args.output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

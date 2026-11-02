#!/usr/bin/env python3
"""
Step 8: computed Raman-active frequencies (black sticks) over the measured YB4 spectra.

The DFT harmonic calculation gives mode frequencies only. Raman intensities need the
Raman tensor (derivative of the polarizability along each mode), which is not computed
here, so all sticks have equal height and only their positions mean anything. Degenerate
partners (Eg) are drawn as one stick. Measured spectra are normalised to their maximum.

    conda activate ELFnet
    python ../scripts/plot_raman_sticks.py --csv gamma_modes.csv --out raman_sticks
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wdf_min import read_wdf  # noqa: E402



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="gamma_modes.csv")
    ap.add_argument("--spot1", default=None, help="optional measured spectrum (.wdf), drawn in blue")
    ap.add_argument("--spot2", default=None, help="optional second measured spectrum (.wdf), drawn in orange")
    ap.add_argument("--out", default="raman_sticks")
    ap.add_argument("--xmin", type=float, default=100.0)
    ap.add_argument("--xmax", type=float, default=1400.0)
    args = ap.parse_args()

    peaks = set()
    for r in csv.DictReader(open(args.csv)):
        if r["activity"] == "Raman":
            peaks.add(float(r["freq_cm1"]))
    freqs = sorted(peaks)
    print(f"{len(freqs)} distinct Raman-active peaks")

    fig, ax = plt.subplots(figsize=(9, 4.4), constrained_layout=True)
    for path, color, label in [(args.spot1, "#0072B2", "Experiment, spot 1"),
                               (args.spot2, "#E69F00", "Experiment, spot 2")]:
        if path is None:
            continue
        x, y, _ = read_wdf(path)
        y = y[0]
        m = (x >= args.xmin) & (x <= args.xmax)
        y = y / y[m].max()
        ax.fill_between(x, 0, y, color=color, alpha=0.35, lw=0, zorder=1)
        ax.plot(x, y, color=color, alpha=0.6, lw=0.8, zorder=2)
        ax.fill_between([], [], color=color, alpha=0.35, label=label)
    for i, f in enumerate(freqs):
        ax.vlines(f, 0, 1.0, color="black", lw=1.4, zorder=4)
        ax.text(f, 1.02 + 0.07 * (i % 3), f"{f:.0f}", ha="center", va="bottom", fontsize=7.5, color="black")
    ax.vlines([], 0, 0, color="black", lw=1.4, label="DFT Raman-active modes (positions only)")

    ax.set_xlim(args.xmin, args.xmax)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([])
    ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
    has_exp = args.spot1 is not None or args.spot2 is not None
    ax.set_ylabel("Normalised intensity (experiment)" if has_exp else "Equal heights (positions only)")
    ax.set_title(("YB$_4$: measured Raman spectra and DFT Raman-active modes at $\\Gamma$" if has_exp
                  else "YB$_4$: DFT Raman-active modes at $\\Gamma$ (harmonic, PBE)"), fontsize=10)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17), fontsize=8)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(f"{args.out}.{ext}", dpi=300 if ext == "png" else None)
        print(f"wrote {args.out}.{ext}")


if __name__ == "__main__":
    main()

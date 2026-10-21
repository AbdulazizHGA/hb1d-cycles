"""Plot results/phase_diagram_k100.json as figures/phase_diagram_1d_k100.png.

    python scripts/plot_phase_diagram.py
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from hb1d.geometry import tuning_to_unrolled  # noqa: E402
from hb1d.tunings import CONV_GAP, FAN7, FAN9, PROP13  # noqa: E402
from scripts.phase_diagram import gap_top  # noqa: E402

SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
CLASSES = [  # (label, color); order = index in the colormap
    (r"no cycle of length $\leq 30$ in any dimension (numerical)", "#e4e3df"),
    (r"in $\Omega_{\rm aug}$: 1D cycle found (exact)", "#2a78d6"),
    (r"cycle in higher dimension (mixed condition, $K \leq 10$); no 1D cycle (exact)", "#eb6834"),
    (r"cycle in higher dimension only at $10 < K \leq 30$; 1D not searched", "#1baf7a"),
    ("on a region boundary: exact search unresolved", "#8a8985"),
]
IDX = {"omega1": 0, "nocycle": 0, "aug": 1, "fan": 2, "fan_long": 3, "unresolved": 4}


def main():
    d = json.load(open(os.path.join(HERE, "..", "results", "phase_diagram_k100.json")))
    cells, H = d["cells"], d["meta"]["H_LEVELS"]
    rhos = sorted({c["rho"] for c in cells})
    grid = np.full((len(H), len(rhos)), np.nan)
    tops = {r: gap_top(r) for r in rhos}           # cells store delta; h = delta / delta_top(rho)
    for c in cells:
        i = rhos.index(c["rho"])
        h = c["delta"] / tops[c["rho"]]
        j = min(range(len(H)), key=lambda k: abs(H[k] - h))
        if c["cls"] == "CANDIDATE":
            raise SystemExit(f"candidate cell present, inspect it first: {c}")
        grid[j, i] = IDX.get(c["cls"], np.nan)
    r_edges = np.r_[rhos[0] - (rhos[1] - rhos[0]) / 2, (np.array(rhos[1:]) + rhos[:-1]) / 2,
                    rhos[-1] + (rhos[1] - rhos[0]) / 2]
    Hm = np.array(H)
    h_edges = np.r_[Hm[0] - (Hm[1] - Hm[0]) / 2, (Hm[1:] + Hm[:-1]) / 2, Hm[-1] + (Hm[-1] - Hm[-2]) / 2]

    fig, ax = plt.subplots(figsize=(10, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    cmap = ListedColormap([c for _, c in CLASSES])
    ax.pcolormesh(r_edges, h_edges, np.ma.masked_invalid(grid), cmap=cmap, vmin=-0.5, vmax=len(CLASSES) - 0.5,
                  edgecolors=SURFACE, linewidth=0.6)
    for K in (3, 4, 5, 6):
        ax.axvline(1 / K, color=INK2, lw=0.8, ls=(0, (2, 3)))
        ax.text(1 / K, h_edges[-1] + 0.01, f"1/{K}", ha="center", va="bottom", fontsize=9, color=INK2)
    for K, lo, hi in ((3, 1 / 4, 1 / 3), (4, 1 / 5, 1 / 4), (5, 1 / 6, 1 / 5)):
        ax.text((lo + hi) / 2, 0.02, f"gap $I_{K}$", ha="center", va="bottom", fontsize=10, color=INK)
    for t, lab in ((PROP13, "Prop. 13"), (CONV_GAP, "conv-gap"), (FAN7, "FAN7"), (FAN9, "FAN9")):
        r, dl = tuning_to_unrolled(t)
        h = dl / gap_top(r)
        ax.plot(r, h, marker="o", ms=8, mfc="white", mec=INK, mew=1.5, zorder=5)
        ax.annotate(lab, (r, h), xytext=(8, -4), textcoords="offset points", fontsize=9, color=INK, zorder=6,
                    va="center", bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.85))
    ax.set_xlim(r_edges[0], r_edges[-1])
    ax.set_ylim(h_edges[0], h_edges[-1])
    ax.set_xlabel(r"rotation number $\rho$ of the tangency on $\partial\Omega_1$", color=INK)
    ax.set_ylabel(r"height $\delta/\delta_{\rm top}(\rho)$ above $\partial\Omega_1$", color=INK)
    ax.set_title(r"Heavy-ball cycles on one-dimensional functions, first-level gaps, $\kappa = L/\mu = 100$",
                 color=INK, fontsize=11, loc="left", pad=22)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2)
    present = {int(v) for v in np.unique(grid[~np.isnan(grid)])}
    ax.legend(handles=[Patch(facecolor=c, edgecolor="none", label=l) for i, (l, c) in enumerate(CLASSES) if i in present],
              loc="upper center",
              bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=8.5)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "..", "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "..", "figures", "phase_diagram_1d_k100.png"), dpi=160, facecolor=SURFACE)


if __name__ == "__main__":
    main()

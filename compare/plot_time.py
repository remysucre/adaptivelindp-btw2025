#!/usr/bin/env python3
"""Fig-15(b)-style companion to plot.py: raw planning time vs. query size.

Reads compare/results/{sms,dp-new}_*.csv (algo,type,n,m,trial,time), averages
trials, and plots planning time (log y) vs. number of relations (log x) for
size-first MCS (SMS) and the adaptive DP kernel ("LinDP" in Fig 15), with a
1-second budget reference line. Each panel spans the n-range where the kernel
has data (it stops once a run exceeds ~1 s).
"""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

# palette (dataviz reference instance, light mode; same entity->hue mapping
# as plot.py: SMS = blue, adaptive DP kernel = yellow)
BLUE = "#2a78d6"
YELLOW = "#eda100"
RED = "#e34948"      # 1-second budget reference line
TEXT1, TEXT2 = "#0b0b0b", "#52514e"
SURFACE = "#fcfcfb"

def load(algo, typ):
    d = defaultdict(list)
    path = os.path.join(RES, f"{algo}_{typ}.csv")
    with open(path) as fh:
        for row in csv.DictReader(fh):
            d[int(row["n"])].append(float(row["time"]))
    return dict(sorted((n, sum(v) / len(v)) for n, v in d.items()))

def fmt_count(v, _pos=None):
    if v >= 1e6:
        s = f"{v/1e6:.1f}".rstrip("0").rstrip(".")
        return f"{s}M"
    if v >= 1e3:
        s = f"{v/1e3:.1f}".rstrip("0").rstrip(".")
        return f"{s}K"
    return f"{int(v)}"

def main():
    types = ["chain", "clique", "star", "tree"]
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.4), dpi=160)
    fig.patch.set_facecolor(SURFACE)

    for ax, typ in zip(axes.flat, types):
        ax.set_facecolor(SURFACE)
        sms = load("sms", typ)
        dp = load("dp-new", typ)
        xmax = max(dp)  # kernel stops once a run exceeds ~1 s
        for label, color, d in [("size-first MCS", BLUE, sms), ("LinDP", YELLOW, dp)]:
            xs = [n for n in d if n <= xmax]
            ys = [d[n] for n in xs]
            ax.plot(xs, ys, color=color, linewidth=1.6, solid_capstyle="round")
            ax.annotate(f" {label}", (xs[-1], ys[-1]),
                        textcoords="offset points", xytext=(2, 0),
                        fontsize=7.5, color=TEXT1, va="center", ha="left")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(1.0, color=RED, linewidth=0.9, linestyle=(0, (1.5, 2)))
        ax.annotate("1 second", (3, 1.0), textcoords="offset points",
                    xytext=(0, 4), fontsize=7.5, color=RED, va="bottom")
        ax.set_xlim(2.6, xmax * 3.2)  # room for direct labels
        ax.set_ylim(top=3.0)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_count))
        ax.set_title(typ, fontsize=10, color=TEXT1)
        ax.tick_params(colors=TEXT2, labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#d8d6cf")
        ax.grid(True, which="major", color="#e8e6df", linewidth=0.6)
        ax.set_axisbelow(True)

    fig.supxlabel("Number of relations (log scale)", fontsize=10, color=TEXT1)
    fig.supylabel("Planning time in seconds (log scale)", fontsize=10, color=TEXT1)
    handles = [
        matplotlib.lines.Line2D([], [], color=BLUE, lw=1.6, label="size-first MCS"),
        matplotlib.lines.Line2D([], [], color=YELLOW, lw=1.6,
                                label='LinDP (adaptive DP kernel, Fig 15 "new")'),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.95))
    out = os.path.join(HERE, "fig15b_planning_time.png")
    fig.savefig(out, facecolor=SURFACE)
    print("wrote", out)

if __name__ == "__main__":
    main()

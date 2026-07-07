#!/usr/bin/env python3
"""Fig-15-style comparison: max query size optimizable within a time budget.

Reads compare/results/*.csv (algo,type,n,m,trial,time), aggregates trials by
mean, takes the best transfer variant per n for adaptive LinDP, and plots
number of nodes (log y) vs. time budget (linear x, 0..1s) for each graph type.
"""
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

# palette (dataviz reference instance, light mode)
BLUE = "#2a78d6"     # series 1: SMS (ours)
AQUA = "#1baf7a"     # series 2: adaptive LinDP (full optimizer)
YELLOW = "#eda100"   # series 3: adaptive DP kernel (Fig 15 "new")
GRAY = "#9a9891"     # recessive reference: old DP kernel (Fig 15 "old")
TEXT1, TEXT2 = "#0b0b0b", "#52514e"
SURFACE = "#fcfcfb"

def load():
    # data[(algo, type)][n] = mean time
    raw = defaultdict(lambda: defaultdict(list))
    for f in glob.glob(os.path.join(RES, "*.csv")):
        with open(f) as fh:
            for row in csv.DictReader(fh):
                raw[(row["algo"], row["type"])][int(row["n"])].append(float(row["time"]))
    return {k: {n: sum(v) / len(v) for n, v in d.items()} for k, d in raw.items()}

def curve(series):
    """(n -> time) dict -> monotone step curve (budget, max n)."""
    pts = sorted(series.items())  # by n
    out = []
    best_t = None
    for n, t in pts:
        out.append((t, n))
    # enforce monotonicity in time (keep the cheapest time seen for growing n)
    out.sort()
    xs, ys, ymax = [], [], 0
    for t, n in out:
        ymax = max(ymax, n)
        xs.append(t)
        ys.append(ymax)
    return xs, ys

def best_lindp(data, typ):
    """min time across transfer variants per n"""
    merged = defaultdict(list)
    for algo in ("lindp-new-sorted", "lindp-new-none", "lindp-new-basic"):
        for n, t in data.get((algo, typ), {}).items():
            merged[n].append(t)
    return {n: min(ts) for n, ts in merged.items()}

def fmt_count(v, _pos=None):
    if v >= 1e6:
        s = f"{v/1e6:.1f}".rstrip("0").rstrip(".")
        return f"{s}M"
    if v >= 1e3:
        s = f"{v/1e3:.1f}".rstrip("0").rstrip(".")
        return f"{s}K"
    return f"{int(v)}"

def main():
    data = load()
    types = ["chain", "clique", "star", "tree"]
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.4), dpi=160)
    fig.patch.set_facecolor(SURFACE)

    for ax, typ in zip(axes.flat, types):
        ax.set_facecolor(SURFACE)
        series = [
            ("SMS (ours)", BLUE, data.get(("sms", typ), {}), "-"),
            ("adaptive LinDP", AQUA, best_lindp(data, typ), "-"),
            ("adaptive DP kernel", YELLOW, data.get(("dp-new", typ), {}), "-"),
            ("old DP kernel", GRAY, data.get(("dp-old", typ), {}), "--"),
        ]
        for label, color, d, ls in series:
            if not d:
                continue
            xs, ys = curve(d)
            # clip to the 0..1.05 s budget window
            xs2 = [x for x in xs if x <= 1.05]
            ys2 = ys[: len(xs2)]
            if not xs2:
                continue
            lw = 1.6 if color != GRAY else 1.1
            ax.plot(xs2, ys2, ls, color=color, linewidth=lw, solid_capstyle="round")
            ax._end_labels = getattr(ax, "_end_labels", [])
            ax._end_labels.append((xs2[-1], ys2[-1], label, color))
        ax.set_yscale("log")
        # direct end labels, dodged apart vertically in log space
        import math
        labels = sorted(getattr(ax, "_end_labels", []), key=lambda t: math.log10(t[1]))
        MIN_GAP = 0.22  # decades
        adj = [math.log10(y) for _, y, _, _ in labels]
        for i in range(1, len(adj)):
            if adj[i] - adj[i - 1] < MIN_GAP:
                adj[i] = adj[i - 1] + MIN_GAP
        if adj:  # keep labels inside the axis
            top = math.log10(ax.get_ylim()[1]) - 0.08
            shift = max(0.0, adj[-1] - top)
            adj = [a - shift for a in adj]
        for (x, y, label, color), ly in zip(labels, adj):
            ax.annotate(
                f" {label}: {fmt_count(y)}",
                (x, 10 ** ly),
                textcoords="offset points", xytext=(2, 0),
                fontsize=7.5, color=TEXT1 if color != GRAY else TEXT2,
                va="center", ha="left",
            )
        ax.set_xlim(-0.02, 1.45)  # room for direct labels
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0 s", "250 ms", "500 ms", "750 ms", "1 s"])
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_count))
        ax.set_title(typ, fontsize=10, color=TEXT1)
        ax.tick_params(colors=TEXT2, labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#d8d6cf")
        ax.grid(True, which="major", axis="y", color="#e8e6df", linewidth=0.6)
        ax.set_axisbelow(True)

    fig.supxlabel("Time budget for optimization", fontsize=10, color=TEXT1)
    fig.supylabel("Number of relations (log scale)", fontsize=10, color=TEXT1)
    handles = [
        matplotlib.lines.Line2D([], [], color=BLUE, lw=1.6, label="SMS (ours)"),
        matplotlib.lines.Line2D([], [], color=AQUA, lw=1.6, label="adaptive LinDP (full optimizer)"),
        matplotlib.lines.Line2D([], [], color=YELLOW, lw=1.6, label='adaptive DP kernel (Fig 15 "new")'),
        matplotlib.lines.Line2D([], [], color=GRAY, lw=1.1, ls="--", label='old DP kernel (Fig 15 "old")'),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.95))
    out = os.path.join(HERE, "fig15_comparison.png")
    fig.savefig(out, facecolor=SURFACE)
    print("wrote", out)

    # summary table: max n within 1s budget
    cols = [
        ("SMS (ours)", lambda t: data.get(("sms", t), {})),
        ("adaptLinDP", best_lindp),
        ("DP-new", lambda t: data.get(("dp-new", t), {})),
        ("DP-old", lambda t: data.get(("dp-old", t), {})),
    ]
    print(f"\nMax relations within a 1s budget:")
    print(f"{'type':<8}" + "".join(f"{name:>14}" for name, _ in cols)
          + f"{'vs LinDP':>12}{'vs DP-new':>12}")
    for typ in types:
        vals = []
        for _, get in cols:
            d = get(data, typ) if get is best_lindp else get(typ)
            xs, ys = curve(d) if d else ([], [])
            vals.append(max((y for x, y in zip(xs, ys) if x <= 1.0), default=0))
        r1 = vals[0] / vals[1] if vals[1] else float("nan")
        r2 = vals[0] / vals[2] if vals[2] else float("nan")
        print(f"{typ:<8}" + "".join(f"{fmt_count(v):>14}" for v in vals)
              + f"{r1:>11.0f}x{r2:>11.1f}x")

if __name__ == "__main__":
    main()

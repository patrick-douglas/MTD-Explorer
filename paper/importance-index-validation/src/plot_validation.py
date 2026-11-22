"""Export editable Supplementary Figure S2 and its complete panel tables."""
from pathlib import Path
import csv
import gzip
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)


def read(name):
    p = DATA / name
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write(name, rows):
    with (OUT / name).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)


scenarios = read("scenario_metrics.csv")
metrics = read("bootstrap_metrics.csv.gz")
taxon_summary = read("bootstrap_taxon_summary.csv")
relative_top = sorted([r for r in taxon_summary if r["representation"] == "relative_percent"
                       and r["model"] == "default" and int(r["reference_rank"]) <= 20],
                      key=lambda r: int(r["reference_rank"]))
by_key = {(r["representation"], r["model"], r["taxon"]): r for r in taxon_summary}
COLORS = {"relative_percent": "#0C6D7E", "estimated_reads": "#4D94C4"}
LABELS = {"relative_percent": "Relative abundance", "estimated_reads": "Estimated read counts"}
TEXT, MUTED, GRID = "#20303B", "#52616D", "#E5EAED"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "text.color": TEXT,
    "axes.labelcolor": TEXT, "axes.edgecolor": "#84929A", "axes.linewidth": .65,
    "xtick.color": MUTED, "ytick.color": MUTED, "svg.fonttype": "none",
    "svg.hashsalt": "MTD-Explorer-S2-importance-20261120", "pdf.fonttype": 42,
    "figure.facecolor": "white", "savefig.facecolor": "white"})
fig = plt.figure(figsize=(7.4, 8.3))
fig.text(.055, .978, "Sensitivity and internal stability of taxon prioritization", fontsize=11, weight="bold")
handles = [Line2D([], [], marker="o" if r == "relative_percent" else "D", color=COLORS[r], linestyle="none", markersize=4, label=LABELS[r]) for r in COLORS]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.055, .958), frameon=False,
           ncol=2, fontsize=8, handlelength=1, columnspacing=1.8)


def heading(letter, title, x, y):
    fig.text(x, y, letter, fontsize=11, weight="bold")
    fig.text(x + .032, y, title, fontsize=9, weight="bold")


def clean(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=3, width=.65, labelsize=7.5)
    ax.set_axisbelow(True)


heading("A", "Joint changes in component weights", .055, .905)
heading("B", "Stability after library resampling", .555, .905)
a = fig.add_axes([.105, .670, .365, .200])
b = fig.add_axes([.605, .670, .365, .200])
panel_a, panel_b = [], []
for j, rep in enumerate(COLORS):
    rows = [r for r in scenarios if r["representation"] == rep and r["family"] == "joint_weight_grid"]
    values = np.array([int(r["top20_retained"]) for r in rows])
    bins = np.arange(14, 21)
    fractions = np.array([(values == k).mean() for k in bins])
    a.bar(bins + (j - .5) * .32, fractions * 100, width=.30, color=COLORS[rep], zorder=3)
    panel_a += [{"representation": rep, "top20_retained": int(k), "n_weight_vectors": int((values == k).sum()),
                 "fraction_of_79_vectors": float(f)} for k, f in zip(bins, fractions)]
    for model, style in (("default", "-"), ("no_total", "--")):
        rows_b = [r for r in metrics if r["representation"] == rep and r["model"] == model]
        values_b = np.array([int(r["top20_retained"]) for r in rows_b])
        bins_b = np.arange(13, 21)
        frac_b = np.array([(values_b == k).mean() for k in bins_b])
        b.plot(bins_b, frac_b * 100, color=COLORS[rep], linestyle=style,
               marker=("o" if rep == "relative_percent" else "D") if model == "default" else None, markersize=3,
               linewidth=1.1 if model == "default" else .9)
        panel_b += [{"representation": rep, "model": model, "top20_retained": int(k),
                     "n_bootstrap_draws": int((values_b == k).sum()), "fraction_of_5000_draws": float(f)}
                    for k, f in zip(bins_b, frac_b)]
a.set(xlim=(13.5, 20.5), xticks=np.arange(14, 21), ylim=(0, 65), xlabel="Original top-20 features retained", ylabel="Weight vectors (%)")
b.set(xlim=(12.8, 20.2), xticks=np.arange(13, 21), ylim=(0, 65), xlabel="Model's original top-20 features retained", ylabel="Bootstrap draws (%)")
for ax in (a, b):
    ax.set_yticks([0, 20, 40, 60]); ax.yaxis.grid(color=GRID, linewidth=.5); clean(ax)
fig.text(.105, .885, "79 distinct normalized weight vectors", fontsize=7.6, color=MUTED)
fig.text(.605, .885, "5,000 stratified draws per model", fontsize=7.6, color=MUTED)
fig.text(.605, .615, "Solid: default weights   Dashed: total omitted", fontsize=7.1, color=MUTED)

fig.text(.055, .570, "Original relative-abundance top 20", fontsize=8.5, weight="bold")
heading("C", "Top-20 inclusion", .450, .570)
heading("D", "Bootstrap rank range", .735, .570)
c = fig.add_axes([.465, .085, .215, .450])
d = fig.add_axes([.765, .085, .205, .450])
y = np.arange(20)[::-1]
panel_c, panel_d = [], []
for i, (row, yi) in enumerate(zip(relative_top, y)):
    for ax in (c, d):
        ax.axhline(yi, color=GRID, linewidth=.4, zorder=0)
    fy = .085 + .450 * ((yi + .7) / 20.4)
    fig.text(.038, fy, row["reference_rank"], ha="right", va="center", fontsize=7.4, color=MUTED)
    fig.text(.055, fy, row["taxon"], ha="left", va="center", fontsize=7.4)
    for rep, shift in (("relative_percent", .16), ("estimated_reads", -.16)):
        r = by_key[(rep, "default", row["taxon"])]
        f = float(r["top20_frequency"])
        c.scatter(f * 100, yi + shift, s=15, marker="o" if rep == "relative_percent" else "D", color=COLORS[rep], edgecolor="white", linewidth=.3, zorder=3)
        panel_c.append({"taxon": row["taxon"], "original_relative_rank": i + 1, "representation": rep,
                        "top20_inclusion_frequency": f, "eligibility_frequency": float(r["eligibility_frequency"])})
    lo, med, hi = [float(row[k]) for k in ("rank_p025_given_eligible", "rank_median_given_eligible", "rank_p975_given_eligible")]
    d.hlines(yi, lo, hi, color=COLORS["relative_percent"], linewidth=1)
    d.scatter(med, yi, s=14, color=COLORS["relative_percent"], zorder=3)
    d.scatter(i + 1, yi, s=22, marker="|", color=TEXT, linewidth=.7, zorder=4)
    panel_d.append({"taxon": row["taxon"], "original_relative_rank": i + 1,
                    "rank_p025_given_eligible": lo, "rank_median_given_eligible": med,
                    "rank_p975_given_eligible": hi})
c.set(xlim=(-2, 103), xticks=[0, 50, 100], xlabel="Inclusion in 5,000 draws (%)", ylim=(-.7, 19.7))
d.set(xlim=(0, 40), xticks=[1, 10, 20, 30, 40], xlabel="Relative-abundance rank", ylim=(-.7, 19.7))
for ax in (c, d):
    ax.set_yticks([]); ax.xaxis.grid(color=GRID, linewidth=.5); clean(ax)
d.axvline(20, color="#9DA9AF", linewidth=.7, linestyle="--", zorder=1)
fig.text(.465, .028, "Default weights in both representations", fontsize=7.3, color=MUTED)
fig.text(.735, .025, "Range: 2.5th to 97.5th percentiles", fontsize=7.0, color=MUTED)
fig.text(.735, .010, "Dot: median   Tick: original rank", fontsize=7.0, color=MUTED)

# Preserve editable text and avoid SVG clipping paths for Illustrator imports.
fig.canvas.draw()
for artist in fig.findobj():
    if hasattr(artist, "set_clip_on"):
        artist.set_clip_on(False)
for extension in ("svg", "pdf", "png"):
    fig.savefig(OUT / f"Supplementary_Figure_S2_importance_robustness.{extension}", dpi=450)
write("panel_A_weight_grid.csv", panel_a)
write("panel_B_bootstrap_retention.csv", panel_b)
write("panel_C_top20_inclusion.csv", panel_c)
write("panel_D_rank_stability.csv", panel_d)
(OUT / "plot_software.json").write_text(json.dumps({"matplotlib": matplotlib.__version__, "numpy": np.__version__}, indent=2) + "\n")
print(f"Exported Supplementary Figure S2 and four panel tables to {OUT}")

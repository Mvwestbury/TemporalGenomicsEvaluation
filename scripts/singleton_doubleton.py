"""Singleton:doubleton ratios per substitution class from the HEARTY 5-read draws (Liu birds + black robin).

Input: hearty/raw/<species>/ds5_r1_s1/<sample>.directional.tsv (hearty/directional_singletons.py: each site downsampled to
5 reads with HEARTY's own RNG; k=1 = one non-reference read (singleton), k=2 = two identical non-reference reads (doubleton)).
Complementary directions are pooled (C>T = C>T + G>A, etc.), since read strand is not recorded.
A true heterozygote at 5 reads usually gives k>=2; a single-read error gives k=1, so S:D rises with per-read error/damage.
Output: results/figures/FigS_singleton_doubleton.{pdf,png}, results/figures/figS_singleton_doubleton_data.tsv
"""
import glob, os
import numpy as np, pandas as pd

B = "/path/to/TemporalGenomicsEvaluation"
CLASSES = [("C→T", ["C>T", "G>A"]), ("T→C", ["T>C", "A>G"]),                      # transitions
           ("G→T", ["G>T", "C>A"]), ("G→C", ["G>C", "C>G"]), ("T→A", ["T>A", "A>T"]), ("T→G", ["T>G", "A>C"])]
TS = {"C→T", "T→C"}
SPECIES = [("Fringilla_montifringilla", "Fringilla montifringilla"), ("Actinodura_cyanouroptera", "Actinodura cyanouroptera"),
           ("Trochalopteron_erythrocephalum", "Trochalopteron erythrocephalum"), ("Yuhina_occipitalis", "Yuhina occipitalis"),
           ("Petroica_traversi", "Petroica traversi")]

man = pd.concat([pd.read_csv(f"{B}/hearty/{m}", sep="\t", usecols=[0, 1, 2], names=["species", "sample", "period"], comment=None)
                 for m in ("manifest.tsv", "manifest_robin.tsv")])
man = man[man.species != "species"].drop_duplicates(["species", "sample"])
rows = []
for sp, _ in SPECIES:
    for f in sorted(glob.glob(f"{B}/hearty/raw/{sp}/ds5_r1_s1/*.directional.tsv")):
        s = os.path.basename(f).replace(".directional.tsv", "")
        period = man[(man.species == sp) & (man["sample"] == s)].period.iloc[0]
        t = pd.read_csv(f, sep="\t", comment="#")
        sites = int(open(f).readline().split()[1])
        for name, types in CLASSES:
            k1 = t[(t.k == 1) & t.type.isin(types)]["count"].sum()
            k2 = t[(t.k == 2) & t.type.isin(types)]["count"].sum()
            rows.append(dict(species=sp, sample=s, period=period, cls=name, sites=sites,
                             singletons=k1, doubletons=k2, ratio=k1 / k2))
D = pd.DataFrame(rows)
os.makedirs(f"{B}/results/figures", exist_ok=True)
D.to_csv(f"{B}/results/figures/figS_singleton_doubleton_data.tsv", sep="\t", index=False)
M = D.groupby(["species", "period", "cls"], sort=False).ratio.mean().unstack("cls")[[c for c, _ in CLASSES]]
FC = (M.xs("historical", level="period") / M.xs("modern", level="period"))
if __name__ == "__main__" and os.environ.get("SD_PRINT"):
    pd.set_option("display.width", 200)
    print(M.round(3).to_string()); print("\nfold change historical/modern"); print(FC.round(2).to_string())

# ---------------------------------------------------------------- figure
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
HIST, MOD, INK, MUTED, GRID = "#eb6834", "#2a78d6", "#222222", "#6b6b66", "#e6e6e3"
COL = {"historical": HIST, "modern": MOD}
rng = np.random.default_rng(1)
cls = [c for c, _ in CLASSES]


fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 5, height_ratios=[1.25, 1], hspace=0.62, wspace=0.12, left=0.075, right=0.99, top=0.9, bottom=0.11)
axs = [fig.add_subplot(gs[0, j]) for j in range(5)]
ymin, ymax = 0.5, 100
for ax, (sp, lab) in zip(axs, SPECIES):
    ax.axvspan(1.5, 5.6, color="#f4f4f2", zorder=0, lw=0)
    for i, c in enumerate(cls):
        for p, dx in (("historical", -0.18), ("modern", 0.18)):
            v = D[(D.species == sp) & (D.period == p) & (D.cls == c)].ratio.values
            x = i + dx + rng.uniform(-0.06, 0.06, len(v))
            ax.scatter(x, v, s=9, color=COL[p], edgecolor="white", linewidth=0.4, zorder=3)
            ax.plot([i + dx - 0.13, i + dx + 0.13], [v.mean()] * 2, color=INK, lw=1, zorder=4)
    ax.set_yscale("log"); ax.set_ylim(ymin, ymax); ax.set_xlim(-0.6, 5.6)
    ax.set_xticks(range(6)); ax.set_xticklabels(cls, rotation=90, fontsize=6.5)
    ax.axhline(1, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_title(lab.split()[0][0] + ". " + lab.split()[1], fontsize=7, style="italic", pad=4)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0); ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=0.6, colors=INK)
    if ax is not axs[0]: ax.tick_params(labelleft=False)
    ax.text(0.5, ymax * 0.8, "transitions", ha="center", va="top", fontsize=5.5, color=MUTED)
    ax.text(3.55, ymax * 0.8, "transversions", ha="center", va="top", fontsize=5.5, color=MUTED)
axs[0].set_ylabel("singletons : doubletons\n(5-read draw)")
axs[0].text(-0.32, 1.13, "a", transform=axs[0].transAxes, fontsize=9, fontweight="bold")
fig.legend(handles=[Line2D([], [], marker="o", ls="", color=HIST, mec="white", ms=5, label="historical"),
                    Line2D([], [], marker="o", ls="", color=MOD, mec="white", ms=5, label="modern"),
                    Line2D([], [], color=INK, lw=1, label="mean")],
           loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.53, 1.0), fontsize=7)


# c: historical / modern fold change of S:D, per class (species means)
bot = gs[1, :].subgridspec(1, 3, width_ratios=[1, 1.5, 1], wspace=0)
axc = fig.add_subplot(bot[1])
fc = FC.loc[[s for s, _ in SPECIES], cls]
im = axc.imshow(np.log2(fc.values), cmap="Oranges", vmin=0, vmax=3.5, aspect="auto")
for (r, c), v in np.ndenumerate(fc.values):
    axc.text(c, r, f"{v:.1f}", ha="center", va="center", fontsize=6, color="white" if np.log2(v) > 2.2 else INK)
axc.set_xticks(range(6)); axc.set_xticklabels(cls, rotation=90, fontsize=6.5)
axc.set_yticks(range(5)); axc.set_yticklabels([l.split()[0] for _, l in SPECIES], fontsize=6.5, style="italic")
axc.tick_params(length=0); [s.set_visible(False) for s in axc.spines.values()]
axc.set_title("historical / modern S:D (×)", fontsize=7, pad=4)
axc.text(-0.36, 1.08, "b", transform=axc.transAxes, fontsize=9, fontweight="bold")

for ext in ("pdf", "png"):
    fig.savefig(f"{B}/results/figures/FigS_singleton_doubleton.{ext}", dpi=300)

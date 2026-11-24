"""Main figures for the manuscript.

Fig_main_liu    Liu et al. birds: (a) per-individual heterozygosity relative to the modern mean, by method;
                (b) SCAR degradation of three modern F. montifringilla: apparent change if degraded reads were historical.
Fig_main_robin  Chatham Island black robin: (a) per-individual heterozygosity by method; (b) SCAR degradation of five modern robins.
Fig_main_empirical  alternative split: Liu (a) + robin (b) empirical; Fig_main_scar: SCAR for both species.
Fig_main_roh    copy of Fig3_ROH_modern (scripts/figures.py), as a standalone figure.

Simulation convention (both figures): apparent change = 100 * (H_original / H_degraded - 1), i.e. the change one would
report if the degraded reads were the historical sample and the unmodified reads the modern sample
(the same direction as the empirical historical -> modern comparison).
Writes results/figures/Fig_main_{liu,robin,roh}.{pdf,png} and fig_main_data.tsv.
Run: /path/to/software/miniconda3/bin/python scripts/main_figures.py
"""
import glob, gzip, os, re, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

B = "/path/to/TemporalGenomicsEvaluation"
OUT = f"{B}/results/figures"
SP = ["Fringilla_montifringilla", "Actinodura_cyanouroptera", "Trochalopteron_erythrocephalum", "Yuhina_occipitalis"]
SPL = {s: f"{s.split('_')[0][0]}. {s.split('_')[1]}" for s in SP}
HIST, MOD = "#eb6834", "#2a78d6"
INK, INK2, GRID = "#1f1f1e", "#5f5e58", "#e4e3df"
ARMS = [("frag", "fragmentation", "#b3b2ab"), ("frag_deam", "+ deamination", "#77766f"),
        ("frag_deam_tvex", "+ transversion excess", "#2f2f2c")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})
rng = np.random.default_rng(3)
rows_out = []


def style(ax, ygrid=True):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5)
    if ygrid: ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)


def pct(h, m):
    return 100 * (np.mean(m) / np.mean(h) - 1)


def label(ax, s, x=-0.16, y=1.08):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def dots(ax, x, v, col, w=0.09, s=11):
    ax.scatter(x + rng.uniform(-w / 2, w / 2, len(v)), v, s=s, color=col, edgecolor="white", linewidth=0.4, zorder=3)
    ax.plot([x - 0.14, x + 0.14], [np.mean(v)] * 2, color=INK, lw=1, zorder=4)


# ------------------------------------------------------------------ data
E = pd.read_csv(f"{OUT}/fig1_per_individual.tsv", sep="\t")          # Liu birds, methods a-g (scripts/figures.py)
LIU_METHODS = [("a", "Liu et al., as published"), ("b", "ANGSD, native depth"),
               ("c", "ANGSD, depth-matched\n+ mapDamage-rescaled"), ("d", "ANGSD, depth-matched,\ntransversions only"),
               ("e", "HEARTY"), ("g", "ATLAS (5×)")]

# simulations: per-bird estimates (HEARTY, ATLAS, ANGSD all) from the composite-figure table; ANGSD tv from the .het files
S = pd.read_csv(f"{OUT}/fig_liu_robin_data.tsv", sep="\t")
S = S[S.panel == "sim"].rename(columns={"period": "arm", "metric": "method"})[["species", "sample", "arm", "method", "value"]]
tv = []
for sp in ("Fringilla_montifringilla", "Petroica_traversi"):
    for f in glob.glob(f"{B}/scar_sim/{sp}/*/het/tv/*.het"):
        x = open(f).read().split("\t")
        tv.append(dict(species=sp, sample=x[1], arm=x[2], method="angsd_tv", value=float(x[-2])))
S = pd.concat([S, pd.DataFrame(tv)])
SIM_METHODS = [("angsd", "ANGSD"), ("angsd_tv", "ANGSD\ntransversions"), ("hearty", "HEARTY"), ("atlas", "ATLAS")]


def sim_change(sp):
    W = S[S.species == sp].pivot_table(index=["sample", "method"], columns="arm", values="value")
    R = {}
    for a, _, _ in ARMS:
        R[a] = 100 * (W["original"] / W[a] - 1)
    return pd.DataFrame(R).reset_index()


def sim_panel(ax, sp, ref=None):
    C = sim_change(sp)
    for i, (m, _) in enumerate(SIM_METHODS):
        for j, (a, _, col) in enumerate(ARMS):
            v = C[C.method == m][a].values
            dots(ax, i + (j - 1) * 0.26, v, col, w=0.08)
            rows_out.extend(dict(figure=sp, panel="sim", method=m, arm=a, value=x) for x in v)
        mfull = C[C.method == m]["frag_deam_tvex"].mean()
        ax.text(i + 0.26, mfull + (4 if mfull >= 0 else -4), f"{mfull:+.0f}%".replace("-", "−"), ha="center",
                va="bottom" if mfull >= 0 else "top", fontsize=6.5, fontweight="bold", color=INK)
    ax.axhline(0, color=INK2, lw=0.6, zorder=1)
    if ref is not None:
        ax.axhline(ref[0], color=HIST, lw=0.8, ls=(0, (3, 2)), zorder=1)
        ax.text(3.62, ref[0] + 1.5, ref[1], ha="right", va="bottom", fontsize=6, color=HIST)
    ax.set_xticks(range(len(SIM_METHODS))); ax.set_xticklabels([l for _, l in SIM_METHODS])
    ax.set_xlim(-0.55, len(SIM_METHODS) - 0.45)
    ax.set_ylabel("apparent change in heterozygosity (%)\ndegraded (as historical) → unmodified")
    style(ax)
    return C


def arm_legend(ax, loc="lower left", anchor=None):
    h = [Line2D([], [], marker="o", ls="", color=c, mec="white", ms=5, label=l) for _, l, c in ARMS]
    ax.legend(handles=h, loc=loc, bbox_to_anchor=anchor, frameon=True, framealpha=0.9, edgecolor="none", fontsize=6.5, handletextpad=0.2, borderaxespad=0)


# ------------------------------------------------------------------ Fig: Liu
fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.15], hspace=0.62, wspace=0.12,
                      left=0.1, right=0.985, top=0.93, bottom=0.07)
axs = []
for k, (m, title) in enumerate(LIU_METHODS):
    ax = fig.add_subplot(gs[k // 3, k % 3], sharey=axs[0] if axs else None); axs.append(ax)
    for i, sp in enumerate(SP):
        q = E[(E.method == m) & (E.species == sp)]
        ref = q[q.period == "modern"].h.mean()
        h, mo = q[q.period == "historical"].h.values / ref, q[q.period == "modern"].h.values / ref
        dots(ax, i - 0.18, h, HIST); dots(ax, i + 0.18, mo, MOD)
        c = pct(h, mo)
        ax.text(i, 3.05, f"{c:+.0f}%".replace("-", "−"), ha="center", va="top", fontsize=6.5, fontweight="bold", color=INK)
        rows_out.extend(dict(figure="liu", panel=m, species=sp, period=p, value=x)
                        for p, vv in (("historical", h), ("modern", mo)) for x in vv)
    ax.axhline(1, color=INK2, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_title(title, fontsize=7, pad=4, va="bottom")
    ax.set_xticks(range(4)); ax.set_xticklabels([SPL[s] for s in SP], rotation=35, ha="right", style="italic", fontsize=6.3)
    ax.set_xlim(-0.6, 3.6); ax.set_ylim(0, 3.15)
    style(ax)
    if k % 3: ax.tick_params(labelleft=False)
    else: ax.set_ylabel("heterozygosity\n(relative to modern mean)")
label(axs[0], "a", x=-0.3, y=1.12)
fig.legend(handles=[Line2D([], [], marker="o", ls="", color=HIST, mec="white", ms=5, label="historical"),
                    Line2D([], [], marker="o", ls="", color=MOD, mec="white", ms=5, label="modern"),
                    Line2D([], [], color=INK, lw=1, label="mean")],
           loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 1.0), fontsize=7)

sub = gs[2, :].subgridspec(1, 3, width_ratios=[0.22, 2.6, 0.18], wspace=0)
axb = fig.add_subplot(sub[1])
sim_panel(axb, "Fringilla_montifringilla", ref=(-30, "Liu et al. reported −30%"))
axb.set_ylim(-50, 20)
axb.set_title("SCAR: historical-like degradation of three modern F. montifringilla", fontsize=7, pad=4)
arm_legend(axb, loc="lower right")
label(axb, "b", x=-0.2, y=1.06)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/Fig_main_liu.{ext}", dpi=300)
plt.close(fig)

# ------------------------------------------------------------------ Fig: robin
RSP = "Petroica_traversi"
man = pd.read_csv(f"{B}/hearty/manifest_robin.tsv", sep="\t", usecols=[0, 1, 2], names=["species", "sample", "period"])
man = man[man.species == RSP].drop_duplicates("sample")
pat = re.compile(r": (\d+) sites, ([0-9.]+) sites with")


def het_file(f):
    x = open(f).read().split()
    return float(x[-1]) if float(x[-1]) < 1 else float(x[-2])


R = []
for s, per in zip(man["sample"], man.period):
    R.append(("angsd", s, per, het_file(f"{B}/het/all/{RSP}/{s}.het")))
    R.append(("angsd_tv", s, per, het_file(f"{B}/het/tv/{RSP}/{s}.het")))
    f = glob.glob(f"{B}/hearty/raw/{RSP}/ds5_exact/{s}_het0.25*windows.txt.gz")[0]
    n = h = 0.0
    for l in gzip.open(f, "rt"):
        mm = pat.search(l); n += int(mm.group(1)); h += float(mm.group(2))
    R.append(("hearty", s, per, h / n))
    for lev, key in (("native", "atlas_native"), ("5", "atlas")):
        R.append((key, s, per, float(pd.read_csv(f"{B}/atlas/results/{RSP}/{s}/{lev}/heterozygosity_stats.txt.gz",
                                                 sep="\t").HKY85_het.iloc[0])))
R = pd.DataFrame(R, columns=["method", "sample", "period", "h"])
ROBIN_METHODS = [("angsd", "ANGSD"), ("angsd_tv", "ANGSD\ntransversions"), ("hearty", "HEARTY"),
                 ("atlas_native", "ATLAS\n(native)"), ("atlas", "ATLAS\n(5×)")]

fig = plt.figure(figsize=(7.2, 3.3))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.42, left=0.085, right=0.985, top=0.86, bottom=0.17)
axa = fig.add_subplot(gs[0])
for i, (m, _) in enumerate(ROBIN_METHODS):
    q = R[R.method == m]
    h, mo = q[q.period == "historical"].h.values, q[q.period == "modern"].h.values
    dots(axa, i - 0.18, h, HIST); dots(axa, i + 0.18, mo, MOD)
    axa.text(i, 1.3e-2, f"{pct(h, mo):+.0f}%".replace("-", "−"), ha="center", va="top", fontsize=6.5, fontweight="bold")
    rows_out.extend(dict(figure="robin", panel=m, species=RSP, period=p, value=x)
                    for p, vv in (("historical", h), ("modern", mo)) for x in vv)
axa.set_yscale("log"); axa.set_ylim(3e-5, 1.4e-2)
axa.set_xticks(range(len(ROBIN_METHODS))); axa.set_xticklabels([l for _, l in ROBIN_METHODS], fontsize=6.3)
axa.set_xlim(-0.6, len(ROBIN_METHODS) - 0.4)
axa.set_ylabel("heterozygosity (per site)")
axa.set_title("Empirical: 5 historical vs 5 modern P. traversi", fontsize=7, pad=4)
style(axa); label(axa, "a", x=-0.2, y=1.04)
axa.legend(handles=[Line2D([], [], marker="o", ls="", color=HIST, mec="white", ms=5, label="historical"),
                    Line2D([], [], marker="o", ls="", color=MOD, mec="white", ms=5, label="modern")],
           loc="lower left", frameon=False, fontsize=6.5, handletextpad=0.2)

axb = fig.add_subplot(gs[1])
sim_panel(axb, RSP)
axb.set_ylim(-60, 30)
axb.set_title("SCAR: degradation of 5 modern P. traversi", fontsize=7, pad=4)
arm_legend(axb, loc="center right", anchor=(1.0, 0.42))
label(axb, "b", x=-0.24, y=1.04)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/Fig_main_robin.{ext}", dpi=300)
plt.close(fig)


# ------------------------------------------------------------------ alternative split: empirical (Liu + robin) and SCAR (both)
HM_LEG = [Line2D([], [], marker="o", ls="", color=HIST, mec="white", ms=5, label="historical"),
          Line2D([], [], marker="o", ls="", color=MOD, mec="white", ms=5, label="modern"),
          Line2D([], [], color=INK, lw=1, label="mean")]


def liu_emp(fig, cells):
    axs = []
    for k, (m, title) in enumerate(LIU_METHODS):
        ax = fig.add_subplot(cells[k // 3, k % 3], sharey=axs[0] if axs else None); axs.append(ax)
        for i, sp in enumerate(SP):
            q = E[(E.method == m) & (E.species == sp)]
            ref = q[q.period == "modern"].h.mean()
            h, mo = q[q.period == "historical"].h.values / ref, q[q.period == "modern"].h.values / ref
            dots(ax, i - 0.18, h, HIST); dots(ax, i + 0.18, mo, MOD)
            ax.text(i, 3.05, f"{pct(h, mo):+.0f}%".replace("-", "−"), ha="center", va="top", fontsize=6.5, fontweight="bold")
        ax.axhline(1, color=INK2, lw=0.5, ls=(0, (2, 2)), zorder=1)
        ax.set_title(title, fontsize=7, pad=4, va="bottom")
        ax.set_xticks(range(4)); ax.set_xticklabels([SPL[s] for s in SP], rotation=35, ha="right", style="italic", fontsize=6.3)
        ax.set_xlim(-0.6, 3.6); ax.set_ylim(0, 3.15); style(ax)
        if k % 3: ax.tick_params(labelleft=False)
        else: ax.set_ylabel("heterozygosity\n(relative to modern mean)")
    return axs


def robin_emp(ax):
    for i, (m, _) in enumerate(ROBIN_METHODS):
        q = R[R.method == m]
        h, mo = q[q.period == "historical"].h.values, q[q.period == "modern"].h.values
        dots(ax, i - 0.18, h, HIST); dots(ax, i + 0.18, mo, MOD)
        ax.text(i, 1.3e-2, f"{pct(h, mo):+.0f}%".replace("-", "−"), ha="center", va="top", fontsize=6.5, fontweight="bold")
    ax.set_yscale("log"); ax.set_ylim(3e-5, 1.4e-2)
    ax.set_xticks(range(len(ROBIN_METHODS))); ax.set_xticklabels([l for _, l in ROBIN_METHODS], fontsize=6.3)
    ax.set_xlim(-0.6, len(ROBIN_METHODS) - 0.4)
    ax.set_ylabel("heterozygosity (per site)"); style(ax)


# empirical figure
fig = plt.figure(figsize=(7.2, 8.0))
gs = fig.add_gridspec(2, 1, height_ratios=[2.15, 1], hspace=0.42, left=0.1, right=0.985, top=0.93, bottom=0.065)
top = gs[0].subgridspec(2, 3, hspace=0.62, wspace=0.12)
axs = liu_emp(fig, top)
label(axs[0], "a", x=-0.3, y=1.12)
bot = gs[1].subgridspec(1, 3, width_ratios=[0.04, 2.2, 0.5], wspace=0)
axr = fig.add_subplot(bot[1]); robin_emp(axr)
axr.set_title("P. traversi", fontsize=7.5, style="italic", pad=6)
label(axr, "b", x=-0.115, y=1.04)
fig.legend(handles=HM_LEG, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 1.0), fontsize=7)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/Fig_main_empirical.{ext}", dpi=300)
plt.close(fig)

# SCAR figure
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
fig.subplots_adjust(left=0.1, right=0.985, top=0.84, bottom=0.15, wspace=0.08)
sim_panel(a1, "Fringilla_montifringilla", ref=(-30, "Liu et al. reported −30%"))
sim_panel(a2, RSP)
a2.set_ylabel(""); a2.tick_params(labelleft=False)
a1.set_ylim(-60, 30)
a1.set_title("F. montifringilla (3 modern birds)", fontsize=7.5, style="italic", pad=4)
a2.set_title("P. traversi (5 modern birds)", fontsize=7.5, style="italic", pad=4)
label(a1, "a", x=-0.22, y=1.04); label(a2, "b", x=-0.06, y=1.04)
h = [Line2D([], [], marker="o", ls="", color=c, mec="white", ms=5, label=l) for _, l, c in ARMS] + [Line2D([], [], color=INK, lw=1, label="mean")]
fig.legend(handles=h, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.55, 1.0), fontsize=7)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/Fig_main_scar.{ext}", dpi=300)
plt.close(fig)

# ------------------------------------------------------------------ ROH (unchanged figure, standalone name)
for ext in ("pdf", "png"):
    shutil.copy(f"{OUT}/Fig3_ROH_modern.{ext}", f"{OUT}/Fig_main_roh.{ext}")

pd.DataFrame(rows_out).to_csv(f"{OUT}/fig_main_data.tsv", sep="\t", index=False)
if os.environ.get("MF_PRINT"):
    for sp in ("Fringilla_montifringilla", RSP):
        C = sim_change(sp); print(sp); print(C.groupby("method")[[a for a, _, _ in ARMS]].mean().round(1).to_string())
    print(R.groupby(["method", "period"]).h.mean().unstack().assign(
        pct=lambda d: (100 * (d.modern / d.historical - 1)).round(1)).to_string())

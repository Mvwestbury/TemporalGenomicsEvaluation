"""Composite figure: Liu et al. 2026 birds + Chatham Island black robin, empirical and simulated.

  a  deamination per bird (mapDamage 5' C->T, peak positions 1-3 minus interior 10-25)
  b  transversion excess per bird (TvCount combined excess G->T/C->A over the G->C/C->G control)
  c-f per-individual heterozygosity (log scale) with % change modern vs historical:
      c ANGSD native depth, d ANGSD depth-matched + mapDamage-rescaled, e HEARTY (5 reads, >=2 minor), f ATLAS 5x
  g  relative change (modern vs historical, %) per species and method
  h  absolute change |modern - historical| per species and method (log; filled = decline, open = increase)
  i  SCAR: modern reads degraded with a historical profile, % change vs the same unmodified reads
  j  SCAR artefact (fragmentation + all damage, absolute) vs the empirical historical-modern difference

Damage: Liu = pilot/pilot_damage_table.tsv (bwa mem pilot, 1M fragments); robin = robin/damage (rb_damage.sbatch,
mapDamage/TvCount on 1M reads). Liu heterozygosity = results/figures/fig1_per_individual.tsv (figures.py).
Writes results/figures/Fig_liu_robin.{pdf,png} and results/figures/fig_liu_robin_data.tsv (all plotted values).
Run: /path/to/software/miniconda3/bin/python scripts/figure_liu_robin.py
"""
import csv, glob, gzip, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

B = "/path/to/TemporalGenomicsEvaluation"
OUT = f"{B}/results/figures"
SP = ["Fringilla_montifringilla", "Actinodura_cyanouroptera", "Trochalopteron_erythrocephalum",
      "Yuhina_occipitalis", "Petroica_traversi"]
SPL = {s: s.split("_")[0] for s in SP}
HIST, MOD = "#eb6834", "#2a78d6"                       # validated pair (CVD dE 24.7)
ARMC = ["#b3b2ab", "#77766f", "#2f2f2c"]               # validated ordinal ramp: frag, +deam, +all damage
INK, INK2, GRID = "#1f1f1e", "#5f5e58", "#e4e3df"
METH = [("b", "ANGSD, native depth", "o"), ("c", "ANGSD, depth-matched + rescaled", "s"), ("e", "HEARTY", "^"), ("g", "ATLAS, 5×", "D")]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})
pat = re.compile(r": (\d+) sites, ([0-9.]+) sites with")
rows = []   # long table of everything plotted


def put(**k): rows.append(k)


def hearty_H(f):
    n = h = 0.0
    for l in gzip.open(f, "rt"):
        m = pat.search(l); n += int(m.group(1)); h += float(m.group(2))
    return h / n


def het_file(f):
    """.het files end in (sites, H) or (H, sites) depending on the script; H is the one below 1."""
    a, b = (float(x) for x in open(f).read().split("\t")[-2:])
    return a if a < 1 else b


def atlas(f): return float(pd.read_csv(f, sep="\t").HKY85_het.iloc[0])


# ---------------- damage
dm = pd.read_csv(f"{B}/pilot/pilot_damage_table.tsv", sep="\t")
for r in dm[dm["mode"] == "modern"].itertuples():
    put(panel="damage", species=r.species.replace(" ", "_"), sample=r.voucher, period=r.period, metric="CtoT_excess", value=r.md_CtoT_excess)
    put(panel="damage", species=r.species.replace(" ", "_"), sample=r.voucher, period=r.period, metric="tv_excess", value=r.tvcount_excess_rate)
RS = pd.read_csv(f"{B}/robin/robin_samples.tsv", sep="\t")
for r in RS.itertuples():
    s = f"{r.voucher}_{r.CRR}"; q = f"{B}/robin/damage/{s}"
    if not os.path.exists(f"{q}/tvcount.summary.tsv"): print("no damage yet", s); continue
    ct = pd.read_csv(f"{q}/mapdamage/5pCtoT_freq.txt", sep="\t").set_index("pos").iloc[:, 0]
    od = next(csv.DictReader(open(f"{q}/tvcount.summary.tsv"), delimiter="\t"))
    put(panel="damage", species="Petroica_traversi", sample=r.voucher, period=r.period, metric="CtoT_excess",
        value=ct.loc[1:3].max() - ct.loc[10:25].mean())
    put(panel="damage", species="Petroica_traversi", sample=r.voucher, period=r.period, metric="tv_excess",
        value=float(od["overall_combined_excess_rate"]))

# ---------------- empirical heterozygosity
D = pd.read_csv(f"{OUT}/fig1_per_individual.tsv", sep="\t")
for r in D[D.method.isin([m[0] for m in METH])].itertuples():
    put(panel="het", species=r.species, sample=r.sample, period=r.period, metric=r.method, value=r.h)
sp = "Petroica_traversi"
for r in RS.itertuples():
    s = f"{r.voucher}_{r.CRR}"
    put(panel="het", species=sp, sample=s, period=r.period, metric="b", value=het_file(f"{B}/het/all/{sp}/{s}.het"))
    put(panel="het", species=sp, sample=s, period=r.period, metric="c", value=het_file(f"{B}/het_matched/rescaled/all/{sp}/{s}.het"))
    put(panel="het", species=sp, sample=s, period=r.period, metric="e",
        value=hearty_H(glob.glob(f"{B}/hearty/raw/{sp}/ds5_exact/{s}_het0.25*windows.txt.gz")[0]))
    put(panel="het", species=sp, sample=s, period=r.period, metric="g", value=atlas(f"{B}/atlas/results/{sp}/{s}/5/heterozygosity_stats.txt.gz"))
    put(panel="het", species=sp, sample=s, period=r.period, metric="r", value=het_file(f"{B}/het_atlasrecal/all/{sp}/{s}.het"))

# ---------------- SCAR
SIMM = [("angsd", "ANGSD"), ("recal", "ANGSD\nrecal."), ("hearty", "HEARTY"), ("atlas", "ATLAS")]
ARMS = [("frag", "fragmentation"), ("frag_deam", "+ deamination"), ("frag_deam_tvex", "+ all damage")]
for spp in ("Fringilla_montifringilla", "Petroica_traversi"):
    for d in sorted(glob.glob(f"{B}/scar_sim/{spp}/*/")):
        s = d.rstrip("/").split("/")[-1]
        for arm in ["original"] + [a for a, _ in ARMS]:
            v = {}
            for k, f in (("angsd", f"{d}het/all/{s}_{arm}.het"), ("recal", f"{d}het_atlasrecal/all/{s}_{arm}.het")):
                if os.path.exists(f): v[k] = het_file(f)
            f = glob.glob(f"{d}hearty/{arm}/*_het0.25*exact.windows.txt.gz")
            if f: v["hearty"] = hearty_H(f[0])
            f = f"{d}atlas/{arm}/heterozygosity_stats.txt.gz"
            if os.path.exists(f): v["atlas"] = atlas(f)
            for k, x in v.items(): put(panel="sim", species=spp, sample=s, period=arm, metric=k, value=x)

T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/fig_liu_robin_data.tsv", sep="\t", index=False, float_format="%.6g")


# ---------------- plotting helpers
def letter(ax, l):
    ax.annotate(l, xy=(0, 1), xycoords="axes fraction", xytext=(-26, 6), textcoords="offset points",
                fontsize=9, fontweight="bold", ha="left", va="bottom")


def style(ax, grid="y"):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, lw=0.5, which="major"); ax.set_axisbelow(True)
    ax.tick_params(length=2, pad=1.5)


def species_axis(ax, n=5):
    ax.set_xticks(range(n)); ax.set_xticklabels([SPL[s] for s in SP[:n]], fontstyle="italic", fontsize=5.7)
    ax.set_xlim(-0.55, n - 0.45)
    ax.axvline(3.5, color=GRID, lw=0.8, zorder=0)   # Liu | robin


rng = np.random.default_rng(3)


def strip(ax, q, log=False):
    for j, sp in enumerate(SP):
        for k, (per, col) in enumerate((("historical", HIST), ("modern", MOD))):
            y = q[(q.species == sp) & (q.period == per)].value.values; x0 = j + (k - 0.5) * 0.38
            ax.scatter(x0 + rng.uniform(-0.07, 0.07, len(y)), y, s=13, color=col, edgecolor="white", linewidth=0.5, zorder=3)
            if len(y): ax.hlines(y.mean(), x0 - 0.14, x0 + 0.14, color=INK, lw=1.0, zorder=4)
    if log: ax.set_yscale("log")
    species_axis(ax)


def change(q):
    h, m = q[q.period == "historical"].value.mean(), q[q.period == "modern"].value.mean()
    return 100 * (m - h) / h, m - h


fig = plt.figure(figsize=(7.2, 10.2))
fig.subplots_adjust(left=0.09, right=0.985, top=0.975, bottom=0.07)
gs = fig.add_gridspec(5, 2, height_ratios=[0.85, 1, 1, 0.95, 1.05], hspace=0.5, wspace=0.2)

# a, b damage
dmg = T[T.panel == "damage"]
for i, (met, lab, lt) in enumerate((("CtoT_excess", "5′ C→T excess (deamination)", "a"),
                                     ("tv_excess", "G→T/C→A excess (transversion excess)", "b"))):
    ax = fig.add_subplot(gs[0, i]); strip(ax, dmg[dmg.metric == met])
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_ylabel(lab, fontsize=6.5, color=INK2); style(ax); letter(ax, lt)

# c-f heterozygosity
het = T[T.panel == "het"]
lo, hi = het[het.metric.isin(["b", "c", "e", "g"])].value.min() / 1.6, het.value.max() * 4
for i, (m, lab, _) in enumerate(METH):
    ax = fig.add_subplot(gs[1 + i // 2, i % 2]); q = het[het.metric == m]; strip(ax, q, log=True)
    ax.set_ylim(lo, hi)
    for j, sp in enumerate(SP):
        c, _ = change(q[q.species == sp])
        ax.text(j, hi * 0.82, f"{c:+.0f}%", ha="center", va="top", fontsize=6.5, fontweight="bold")
    ax.set_title(lab, loc="left", fontsize=6.8, color=INK, pad=3)
    if i % 2 == 0: ax.set_ylabel("heterozygosity (per site)", fontsize=6.5, color=INK2)
    style(ax); letter(ax, "cdef"[i])

# g, h change per species x method
axg, axh = fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1])
off = np.linspace(-0.27, 0.27, len(METH))
for o, (m, lab, mk) in zip(off, METH):
    for j, sp in enumerate(SP):
        c, a = change(het[(het.metric == m) & (het.species == sp)])
        axg.scatter(j + o, c, marker=mk, s=18, color=INK, edgecolor="white", linewidth=0.4, zorder=3)
        axh.scatter(j + o, abs(a), marker=mk, s=18, zorder=3, linewidth=0.8,
                    facecolor=INK if a < 0 else "white", edgecolor=INK)
        put(panel="change", species=sp, sample="", period="", metric=f"{m}_pct", value=c)
        put(panel="change", species=sp, sample="", period="", metric=f"{m}_abs", value=a)
axg.axhline(0, color=INK2, lw=0.6); axg.set_ylabel("change, modern vs historical (%)", fontsize=6.5, color=INK2)
axh.set_yscale("log"); axh.set_ylabel("|modern − historical| (per site)", fontsize=6.5, color=INK2)
for ax, lt in ((axg, "g"), (axh, "h")): species_axis(ax); style(ax); letter(ax, lt)
axh.legend(handles=[Line2D([], [], marker="o", ls="", mfc=INK, mec=INK, ms=4, label="decline"),
                    Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=4, label="increase")],
           frameon=False, fontsize=6, loc="lower left", handletextpad=0.2, borderaxespad=0.1)

# i, j simulations
sim = T[T.panel == "sim"]
orig = sim[sim.period == "original"].set_index(["species", "sample", "metric"]).value
sim = sim[sim.period != "original"].copy()
sim["orig"] = [orig.get((a, b, c), np.nan) for a, b, c in zip(sim.species, sim["sample"], sim.metric)]
sim["pct"], sim["abs"] = 100 * (sim.value / sim.orig - 1), sim.value - sim.orig
axi, axj = fig.add_subplot(gs[4, 0]), fig.add_subplot(gs[4, 1])
w = 0.2; xt, xl = [], []
for si, spp in enumerate(("Fringilla_montifringilla", "Petroica_traversi")):
    for mi, (m, mlab) in enumerate(SIMM):
        x = si * (len(SIMM) + 0.8) + mi
        for ai, (arm, _) in enumerate(ARMS):
            q = sim[(sim.species == spp) & (sim.metric == m) & (sim.period == arm)]
            axi.bar(x + (ai - 1) * w, q.pct.mean(), width=w * 0.9, color=ARMC[ai], zorder=3)
            put(panel="sim_mean", species=spp, sample=f"n={len(q)}", period=arm, metric=f"{m}_pct", value=q.pct.mean())
            put(panel="sim_mean", species=spp, sample=f"n={len(q)}", period=arm, metric=f"{m}_abs", value=q["abs"].mean())
        # j: artefact (all damage) vs empirical difference, absolute
        q = sim[(sim.species == spp) & (sim.metric == m) & (sim.period == "frag_deam_tvex")]
        axj.bar(x - 0.17, abs(q["abs"].mean()), width=0.3, color=ARMC[2], zorder=3)
        em = {"angsd": "b", "recal": "r", "hearty": "e", "atlas": "g"}[m]
        qe = het[(het.metric == em) & (het.species == spp)]
        if len(qe):
            _, a = change(qe); axj.bar(x + 0.17, abs(a), width=0.3, color="white", edgecolor=INK, linewidth=0.8, zorder=3)
        xt.append(x); xl.append(mlab)
    for ax in (axi, axj):
        ax.annotate(SPL[spp], xy=(si * (len(SIMM) + 0.8) + 1.5, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -20), textcoords="offset points", ha="center", va="top", fontsize=6.5, fontstyle="italic")
for ax in (axi, axj):
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5.6)
    ax.tick_params(axis="x", length=0, pad=2)
axi.axhline(0, color=INK2, lw=0.6); axi.set_ylabel("degraded vs unmodified reads (%)", fontsize=6.5, color=INK2)
axj.set_yscale("log"); axj.set_ylabel("absolute difference (per site)", fontsize=6.5, color=INK2)
for ax, lt in ((axi, "i"), (axj, "j")): style(ax); letter(ax, lt)
axi.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, label=l) for c, (_, l) in zip(ARMC, ARMS)],
           frameon=False, fontsize=5.8, loc="lower right", bbox_to_anchor=(1, 1.0), ncol=3, handlelength=1, borderaxespad=0.1, columnspacing=0.8)
axj.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=ARMC[2], label="simulated artefact (frag. + all damage)"),
                    plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=INK, label="empirical difference")],
           frameon=False, fontsize=5.8, loc="lower right", bbox_to_anchor=(1, 1.0), ncol=2, handlelength=1, borderaxespad=0.1, columnspacing=0.8)

fig.legend(handles=[Line2D([], [], marker="o", ls="", color=HIST, label="historical (museum)", markersize=4.5),
                    Line2D([], [], marker="o", ls="", color=MOD, label="modern", markersize=4.5),
                    Line2D([], [], color=INK, lw=1, label="mean")] +
                   [Line2D([], [], marker=mk, ls="", color=INK, label=lab, markersize=4) for _, lab, mk in METH],
           loc="lower center", ncol=7, frameon=False, fontsize=6, bbox_to_anchor=(0.5, 0.0), handlelength=1.1, columnspacing=1.0)
pd.DataFrame(rows).to_csv(f"{OUT}/fig_liu_robin_data.tsv", sep="\t", index=False, float_format="%.6g")
for ext in ("pdf", "png"): fig.savefig(f"{OUT}/Fig_liu_robin.{ext}", dpi=300)
print("written", f"{OUT}/Fig_liu_robin.png")

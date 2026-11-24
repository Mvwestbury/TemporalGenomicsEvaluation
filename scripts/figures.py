"""Final figures for the eLetter.

Fig 1  per-individual heterozygosity, historical vs modern, by method (rows a-g) and species (columns).
Fig 2  % change modern vs historical for every analysis incl. population statistics; Tajima's D; SCAR.
Fig 3  ROH in modern birds at native depth: bcftools roh and HEARTY windows (two cutoffs), FROH by segment length.
Writes results/figures/Fig{1,2,3}.{pdf,png} and results/figures/fig1_per_individual.tsv.
Run: /path/to/software/miniconda3/bin/python scripts/figures.py  (after summary_table.py)
"""
import glob, gzip, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/path/to/TemporalGenomicsEvaluation"
OUT = f"{B}/results/figures"; os.makedirs(OUT, exist_ok=True)
SP = ["Fringilla_montifringilla", "Actinodura_cyanouroptera", "Trochalopteron_erythrocephalum", "Yuhina_occipitalis"]
SPL = {s: f"{s.split('_')[0][0]}. {s.split('_')[1]}" for s in SP}
HIST, MOD = "#eb6834", "#2a78d6"          # validated pair (dataviz validator: all checks pass)
INK, INK2, GRID = "#1f1f1e", "#5f5e58", "#e4e3df"
ROBUST, ANGSDC = "#1baf7a", "#8a8980"     # Fig 2: damage-robust methods vs ANGSD-GL-based
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})

man = pd.read_csv(f"{B}/metadata/full_map_samples.tsv", sep="\t")
man["sample"] = man.voucher + "_" + man.CRR
period = dict(zip(man["sample"], man.period))
recs = []  # long table: method, species, sample, period, h


def put(method, sp, sample, per, h):
    recs.append(dict(method=method, species=sp, sample=sample, period=per, h=h))


# a. paper
p = pd.read_csv(f"{B}/metadata/Liu2026_TableS1_h.tsv", sep="\t")
for r in p.itertuples():
    put("a", r.species.replace(" ", "_"), r.voucher, r.period, r.h_paper)
# b. ANGSD native (all sites)
for f in glob.glob(f"{B}/het/all/*/*.het"):
    x = open(f).read().split("\t"); put("b", x[0].replace(" ", "_"), x[1], x[2], float(x[-1]))
# c/d. ANGSD depth-matched + rescaled (both periods from the rescaled, matched runs)
for m, v in (("c", "all"), ("d", "tv")):
    for f in glob.glob(f"{B}/het_matched/rescaled/{v}/*/*.het"):
        x = open(f).read().split("\t"); put(m, x[0].replace(" ", "_"), x[1], x[2], float(x[9]))
# e/f. HEARTY exact, plain BAMs
pat = re.compile(r": (\d+) sites, ([0-9.]+) sites with")
for f in glob.glob(f"{B}/hearty/raw/*/ds5_*/*_het0.25.mincov5_*downsampled_depth5_exact.windows.txt.gz"):
    s, sp = os.path.basename(f).split("_het0.25")[0], f.split("/")[-3]
    if any(r["sample"] == s and r["method"] == "e" for r in recs): continue
    n = h = 0.0
    for l in gzip.open(f, "rt"):
        mm = pat.search(l); n += int(mm.group(1)); h += float(mm.group(2))
    hv = sum(float(pat.search(l).group(2)) for l in gzip.open(f.replace("windows.txt", "windows_trv.txt"), "rt"))
    put("e", sp, s, period[s], h / n); put("f", sp, s, period[s], hv / n)
# g. ATLAS 5x
at = pd.read_csv(f"{B}/atlas/atlas_summary.tsv", sep="\t")
for r in at[at.level.astype(str) == "5"].itertuples():
    put("g", next(s for s in SP if s.startswith(r.species)), r.sample, r.period, r.het)
D = pd.DataFrame(recs)
D.to_csv(f"{OUT}/fig1_per_individual.tsv", sep="\t", index=False)

ROWS = [("a", "Liu et al. 2026, as published (ANGSD)"),
        ("b", "ANGSD, native depth"),
        ("c", "ANGSD, depth-matched + mapDamage-rescaled"),
        ("d", "ANGSD, depth-matched, transversions only"),
        ("e", "HEARTY, 5 reads/site, ≥2 reads for the minor allele"),
        ("f", "HEARTY, transversions only"),
        ("g", "ATLAS (per-sample damage + error model), 5×")]


def pct(q):
    h, m = q[q.period == "historical"].h.mean(), q[q.period == "modern"].h.mean()
    return 100 * (m - h) / h


# ---------------- Fig 1
rng = np.random.default_rng(3)
fig, axs = plt.subplots(len(ROWS), 4, figsize=(7.2, 10.6), gridspec_kw=dict(hspace=1.0, wspace=0.42))
fig.subplots_adjust(left=0.1, right=0.98, top=0.92, bottom=0.05)
for i, (m, lab) in enumerate(ROWS):
    for j, sp in enumerate(SP):
        ax = axs[i, j]; q = D[(D.method == m) & (D.species == sp)]
        for k, (per, col) in enumerate((("historical", HIST), ("modern", MOD))):
            y = q[q.period == per].h.values
            ax.scatter(k + rng.uniform(-0.12, 0.12, len(y)), y, s=16, color=col, edgecolor="white", linewidth=0.6, zorder=3)
            if len(y): ax.hlines(y.mean(), k - 0.28, k + 0.28, color=INK, lw=1.2, zorder=4)
        top = q.h.max() * 1.32
        ax.set_ylim(0, top); ax.set_xlim(-0.6, 1.6)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["hist.", "mod."] if i == len(ROWS) - 1 else ["", ""])
        ax.yaxis.set_major_locator(plt.MaxNLocator(3)); ax.ticklabel_format(axis="y", style="plain")
        ax.tick_params(length=2, pad=1.5)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
        c = pct(q); c = 0.0 if abs(c) < 0.5 else c
        ax.text(0.5, 0.97, f"{c:+.0f}%", transform=ax.transAxes, ha="center", va="top", fontsize=7,
                fontweight="bold", color=INK)
        if i == 0: fig.text((ax.get_position().x0 + ax.get_position().x1) / 2, 0.965, SPL[sp], ha="center", fontstyle="italic", fontsize=8.5)
    axs[i, 0].set_ylabel("heterozygosity", fontsize=6.5, color=INK2)
    pos = axs[i, 0].get_position()
    fig.text(0.02, pos.y1 + 0.022, f"{m}", fontweight="bold", fontsize=9, va="bottom")
    fig.text(0.05, pos.y1 + 0.022, lab, fontsize=7.5, va="bottom")
h1 = plt.Line2D([], [], marker="o", ls="", color=HIST, label="historical (museum)", markersize=5)
h2 = plt.Line2D([], [], marker="o", ls="", color=MOD, label="modern", markersize=5)
h3 = plt.Line2D([], [], color=INK, lw=1.2, label="mean")
fig.legend(handles=[h1, h2, h3], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
for ext in ("pdf", "png"): fig.savefig(f"{OUT}/Fig1_heterozygosity.{ext}", dpi=300)
plt.close(fig)

# ---------------- Fig 2
S = pd.read_csv(f"{B}/results/summary_table.tsv", sep="\t")
SC = {s.split("_")[0]: s for s in SP}


def ch(analysis):
    q = S[S.analysis == analysis]; return {r.species: r.change for r in q.itertuples()}


F2 = [("Per individual, ANGSD", None),
      ("Liu et al. 2026, as published", ch("Liu et al. 2026 (Table S1; ANGSD, as published)"), ANGSDC),
      ("ANGSD, native depth", ch("ANGSD GL, native depth, all sites"), ANGSDC),
      ("ANGSD, depth-matched + rescaled", ch("ANGSD GL, mapDamage-rescaled + depth-matched, all sites"), ANGSDC),
      ("ANGSD, depth-matched, transversions", ch("ANGSD GL, mapDamage-rescaled + depth-matched, transversions only"), ANGSDC),
      ("Population, authors' ANGSD", None),
      ("Authors' H (SFS singletons)", ch("Authors' 'heterozygosity' (singleton class of unfolded SFS)"), ANGSDC),
      ("Watterson's θ", ch("Watterson's theta per site"), ANGSDC),
      ("Nucleotide diversity π", ch("Nucleotide diversity (pi) per site"), ANGSDC),
      ("Per individual, damage-robust", None),
      ("HEARTY, all sites", ch("HEARTY, 5 reads/site, >=2 minor reads, plain BAMs, all sites"), ROBUST),
      ("HEARTY, transversions", ch("HEARTY, 5 reads/site, >=2 minor reads, plain BAMs, transversions only"), ROBUST),
      ("ATLAS, 5×", ch("ATLAS (per-sample PMD + recalibration, GERP sites), all birds downsampled to 5x"), ROBUST)]
fig = plt.figure(figsize=(7.2, 8.2))
fig.subplots_adjust(left=0.3, right=0.985, top=0.94, bottom=0.135)
gs = fig.add_gridspec(2, 4, height_ratios=[2.3, 1.65], hspace=0.3, wspace=0.12)
axa = []
ypos, ylab, y = [], [], 0
for row in F2:
    ylab.append(row[0]); ypos.append(y); y -= 1
lo = min(v for r in F2 if r[1] for v in r[1].values())
for j, sp in enumerate(SP):
    ax = fig.add_subplot(gs[0, j]); axa.append(ax)
    for yy, row in zip(ypos, F2):
        if row[1] is None: continue
        v = row[1][sp]; v = 0.0 if abs(v) < 0.5 else v
        ax.barh(yy, v, height=0.62, color=row[2], zorder=3)
        ax.text(v + (1.5 if v >= 0 else -1.5), yy, f"{v:+.0f}", va="center", ha="left" if v >= 0 else "right", fontsize=6, color=INK)
    ax.axvline(0, color=INK, lw=0.7, zorder=4)
    ax.set_xlim(-72, 30); ax.set_ylim(ypos[-1] - 0.7, 0.7)
    ax.set_yticks(ypos); ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.set_title(SPL[sp], fontstyle="italic", fontsize=8)
    if j == 0:
        for yy, row in zip(ypos, F2):
            fig.text(0, 0, "")  # placeholder
            ax.text(-75, yy, row[0], ha="right", va="center", fontsize=7 if row[1] else 7.2,
                    fontweight="normal" if row[1] else "bold", color=INK if row[1] else INK)
pa = axa[0].get_position()
fig.text(0.01, pa.y1 + 0.012, "a", fontweight="bold", fontsize=10)
fig.text((pa.x0 + axa[-1].get_position().x1) / 2, pa.y0 - 0.03, "% change, modern vs historical",
         ha="center", fontsize=6.5, color=INK2)
# b: Tajima's D
gb = gs[1, :].subgridspec(1, 2, wspace=0.75, width_ratios=[1, 1.15])
axb = fig.add_subplot(gb[0])
for k, sp in enumerate(SP):
    r = S[(S.species == sp) & S.analysis.str.startswith("Tajima")].iloc[0]
    axb.plot([r.historical, r.modern], [k, k], color=GRID, lw=2.5, zorder=2, solid_capstyle="round")
    axb.scatter(r.historical, k, color=HIST, s=30, zorder=3, edgecolor="white", lw=0.6)
    axb.scatter(r.modern, k, color=MOD, s=30, zorder=3, edgecolor="white", lw=0.6)
axb.axvline(0, color=INK, lw=0.7)
axb.set_yticks(range(4)); axb.set_yticklabels([SPL[s] for s in SP], fontstyle="italic"); axb.invert_yaxis()
axb.set_xlabel("Tajima's D (mean of 50 kb windows)", fontsize=6.5, color=INK2); axb.set_xlim(-1.05, 0.1)
for s in ("top", "right"): axb.spines[s].set_visible(False)
axb.tick_params(axis="y", length=0)
axb.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=HIST, label="historical", markersize=5),
                    plt.Line2D([], [], marker="o", ls="", color=MOD, label="modern", markersize=5)],
           frameon=False, loc="lower left", fontsize=6.5)
# c: SCAR - same three degraded modern Fringilla, ANGSD vs the damage-robust methods
axc = fig.add_subplot(gb[1])
steps = [("fragmentation", "fragmentation (historical length distribution)"), ("+ deamination", "+ deamination"),
         ("+ transversion excess", "+ deamination + transversion excess")]
meths = [("ANGSD", "", ANGSDC, ("all", "tv")), ("HEARTY", "HEARTY: ", ROBUST, ("all", "tv")),
         ("ATLAS", "ATLAS: ", ROBUST, ("all",))]
obs = {v: S[(S.species == SC["Fringilla"]) & (S.analysis == f"ANGSD GL, native depth, {lab}")].change.iloc[0]
       for v, lab in (("all", "all sites"), ("tv", "transversions only"))}
LABV = {"all": "all sites", "tv": "transversions only"}
cy, cl, cb, y = [], [], [], 0
for k, (nm, key) in enumerate(steps):
    cb.append((y + 0.75, nm)); 
    for mname, pre, col, vs in meths:
        for o, v in enumerate(vs):
            r = S[S.analysis == f"{pre}{key}, {LABV[v]} (change vs same birds unmodified)"].iloc[0]
            yy = y + (o - (len(vs) - 1) / 2) * 0.44
            axc.barh(yy, r.change, height=0.4, color=col, zorder=3, hatch="//////" if v == "tv" else None,
                     edgecolor="white", linewidth=0)
            axc.text(r.change + (2 if r.change >= 0 else -2), yy, f"{r.change:+.0f}", va="center",
                     ha="left" if r.change >= 0 else "right", fontsize=5.8)
        cy.append(y); cl.append(mname); y -= 1.15
    y -= 0.8
for yy, nm in cb:
    axc.text(-0.02, yy, nm, transform=axc.get_yaxis_transform(), ha="right", va="center", fontsize=6.5, fontweight="bold")
for v, ls in (("all", "--"), ("tv", ":")):
    axc.axvline(-obs[v] / (1 + obs[v] / 100), color=INK, lw=0.8, ls=ls, zorder=4)
axc.axvline(0, color=INK, lw=0.7, zorder=4)
axc.set_yticks(cy); axc.set_yticklabels(cl, fontsize=6); axc.set_ylim(y + 0.5, 1.2)
axc.tick_params(axis="y", length=0); axc.set_xlim(-12, 100)
axc.grid(axis="x", color=GRID, lw=0.5); axc.set_axisbelow(True)
for s in ("top", "right", "left"): axc.spines[s].set_visible(False)
axc.set_xlabel("% change, degraded vs unmodified\n(same 3 modern F. montifringilla)", fontsize=6.5, color=INK2)
hatch = plt.Rectangle((0, 0), 1, 1, facecolor=INK2, hatch="//////", edgecolor="white", linewidth=0, label="transversions")
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=ANGSDC, label="ANGSD"),
                    plt.Rectangle((0, 0), 1, 1, color=ROBUST, label="HEARTY / ATLAS"),
                    plt.Rectangle((0, 0), 1, 1, color=INK2, label="all sites"), hatch,
                    plt.Line2D([], [], color=INK, ls="--", lw=0.8, label="observed hist. excess, all sites"),
                    plt.Line2D([], [], color=INK, ls=":", lw=0.8, label="observed hist. excess, transversions")],
           frameon=False, loc="lower center", bbox_to_anchor=(0.62, 0.0), ncol=3, fontsize=6)
pb, pc = axb.get_position(), axc.get_position()
fig.text(0.01, pc.y1 + 0.012, "b", fontweight="bold", fontsize=10)
fig.text(pc.x0 - 0.17, pc.y1 + 0.012, "c", fontweight="bold", fontsize=10)
for ext in ("pdf", "png"): fig.savefig(f"{OUT}/Fig2_changes.{ext}", dpi=300)
plt.close(fig)

# ---------------- Fig 3 ROH
R = pd.read_csv(f"{B}/hearty/native_roh_summary.tsv", sep="\t")
# bcftools roh (joint calls of the 5 modern birds per species, AF estimated from them); file is in %, R in fractions
bc = pd.read_csv(f"{B}/bcftools_roh/bcftools_froh_summary.tsv", sep="\t")
bc = bc[bc.AF == "AFest"].assign(species=lambda x: x.species.map(SC), threshold_mode="bcftools")
bc = bc.rename(columns={k: "FROH_" + k for k in ("gt100kb", "gt1Mb", "gt5Mb")})
for k in ("FROH_gt100kb", "FROH_gt1Mb", "FROH_gt5Mb"): bc[k] /= 100
R = pd.concat([R, bc[R.columns.intersection(bc.columns)]], ignore_index=True)
dep = pd.read_csv(f"{B}/hearty/manifest.tsv", sep="\t")
dep = dict(zip(dep["sample"], dep.mean_depth))
CL = [("0.1–1 Mb", "#9ec5f4"), ("1–5 Mb", "#3987e5"), ("> 5 Mb", "#104281")]  # one-hue sequential steps
fig, axs = plt.subplots(1, 3, figsize=(7.2, 4.5), sharey=True, gridspec_kw=dict(wspace=0.1))
for a, (mode, title) in enumerate((("bcftools", "a   bcftools roh\n     (genotype calls, AF from 5 birds)"),
                                   ("abs1e-4", "b   HEARTY windows,\n     heterozygosity ≤ 10⁻⁴"),
                                   ("median0.3", "c   HEARTY windows,\n     heterozygosity ≤ 0.3 × median"))):
    ax = axs[a]; yt, yl, y = [], [], 0
    for sp in SP:
        q = R[(R.species == sp) & (R.threshold_mode == mode)].copy()
        q["d"] = q["sample"].map(lambda s: next(v for k, v in dep.items() if k.startswith(s)))
        for r in q.sort_values("d").itertuples():
            seg = [r.FROH_gt100kb - r.FROH_gt1Mb, r.FROH_gt1Mb - r.FROH_gt5Mb, r.FROH_gt5Mb]
            left = 0
            for (nm, col), v in zip(CL, seg):
                if v > 0: ax.barh(y, 100 * v, left=left, height=0.72, color=col, edgecolor="white", linewidth=0.8, zorder=3)
                left += 100 * v
            if sum(seg) < 0.001: ax.text(0.1, y, "< 0.1%", va="center", fontsize=5.5, color=INK2)
            yt.append(y); yl.append(f"{r.sample.split('_CRR')[0]} ({r.d:.0f}×)"); y -= 1
        if a == 0: ax.text(-0.35, yt[-len(q)] + 1, SPL[sp], transform=ax.get_yaxis_transform(), ha="left", va="center",
                           fontstyle="italic", fontsize=7, fontweight="bold")
        y -= 1.2
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_title(title, loc="left", fontsize=7.5)
    ax.set_xlabel("FROH (% of autosomes in ROH)", fontsize=6.5, color=INK2)
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
for ax in axs: ax.set_xlim(0, 12)
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, label=n) for n, c in CL], title="ROH segment length",
           loc="lower center", ncol=3, frameon=False, fontsize=6.5, title_fontsize=6.5, bbox_to_anchor=(0.55, -0.01))
fig.subplots_adjust(left=0.2, right=0.98, top=0.93, bottom=0.17)
for ext in ("pdf", "png"): fig.savefig(f"{OUT}/Fig3_ROH_modern.{ext}", dpi=300)
plt.close(fig)
print("written to", OUT)

"""Composite of mapDamage and TvCount results for the pilot (1M fragments per sample).

Writes pilot/pilot_damage_table.tsv (both mapping modes) and
pilot/pilot_damage_composite.{pdf,png} (modern = bwa mem mode only).
"""
import csv
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/path/to/TemporalGenomicsEvaluation"
NPOS = 25
INTERIOR = range(10, 26)          # mapDamage positions used as the interior baseline
COL = {"historical": "#eb6834", "modern": "#2a78d6"}
LABEL = {"historical": "historical", "modern": "modern control"}
INK, INK2 = "#0b0b0b", "#52514e"


def first_row(path):
    with open(path) as f:
        return next(csv.DictReader(f, delimiter="\t"))


def md_freq(path):
    return pd.read_csv(path, sep="\t").set_index("pos").iloc[:, 0]


samples = pd.read_csv(f"{B}/metadata/pilot_samples.tsv", sep="\t")
rows, curves = [], []
for r in samples.itertuples():
    sp = r.species.replace(" ", "_")
    pre = f"{r.voucher}_{r.CRR}"
    d = f"{B}/pilot/{r.mode}/{sp}"
    if not os.path.exists(f"{d}/tvcount/{pre}.summary.tsv"):     # e.g. Actinodura ancient-all, not rerun on GLGS5798
        print(f"skip {r.mode} {sp} {pre}: no output", file=sys.stderr)
        continue
    st = first_row(f"{d}/{pre}.plainmap.stats.tsv")
    od = first_row(f"{d}/tvcount/{pre}.summary.tsv")
    ct = md_freq(f"{d}/{pre}_mapdamage/5pCtoT_freq.txt")
    ga = md_freq(f"{d}/{pre}_mapdamage/3pGtoA_freq.txt")
    pos = pd.read_csv(f"{d}/tvcount/{pre}.positions.tsv", sep="\t")
    ct_int, ga_int = ct.loc[list(INTERIOR)].mean(), ga.loc[list(INTERIOR)].mean()
    ct_pk, ga_pk = ct.loc[1:3].max(), ga.loc[1:3].max()      # the profiles peak at position 2, not 1
    rows.append({
        "species": r.species, "voucher": r.voucher, "CRR": r.CRR, "period": r.period, "mode": r.mode,
        "endogenous": float(st["endog"]), "dup_rate": float(st["duprate"]), "mean_read_len": float(st["avg_readlen"]),
        "md_CtoT_5p_peak": ct_pk, "md_CtoT_interior": ct_int, "md_CtoT_excess": ct_pk - ct_int,
        "md_GtoA_3p_peak": ga_pk, "md_GtoA_interior": ga_int, "md_GtoA_excess": ga_pk - ga_int,
        "tvcount_p_GtoT": float(od["overall_p_gt_over_g"]), "tvcount_p_GtoC": float(od["overall_p_gc_over_g"]),
        "tvcount_p_CtoA": float(od["overall_p_ca_over_c"]), "tvcount_p_CtoG": float(od["overall_p_cg_over_c"]),
        "tvcount_GT_over_GC": float(od["overall_gt_over_gc"]), "tvcount_CA_over_CG": float(od["overall_ca_over_cg"]),
        "tvcount_excess_rate": float(od["overall_combined_excess_rate"]),
        "tvcount_5p_terminal_interior_fold": float(od["5p_terminal_interior_fold"]),
    })
    if r.mode == "modern":
        p5 = pos[pos.End == "5p"].set_index("Pos")
        curves.append((r.species, r.period, ct.loc[1:NPOS], ga.loc[1:NPOS],
                       p5.loc[1:NPOS, "excess_gt_rate"], rows[-1]["tvcount_excess_rate"]))

tab = pd.DataFrame(rows).sort_values(["species", "mode", "period", "voucher"])
tab.to_csv(f"{B}/pilot/pilot_damage_table.tsv", sep="\t", index=False, float_format="%.5g")

species = sorted(samples.species.unique())
titles = ["mapDamage 5′ C→T", "mapDamage 3′ G→A", "TvCount excess G→T (G→T − G→C), 5′ end",
          "TvCount excess rate"]
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(len(species), 4, figsize=(13, 2.4 * len(species)),
                       gridspec_kw={"width_ratios": [1, 1, 1, 0.55]})
for i, sp in enumerate(species):
    sc = [c for c in curves if c[0] == sp]
    for j in range(3):
        a = ax[i, j]
        for period in ("modern", "historical"):        # historical drawn on top
            for c in (c for c in sc if c[1] == period):
                y = c[2 + j]
                a.plot(y.index, y.values, color=COL[period], lw=1.2, alpha=0.9)
        if j == 2:
            a.axhline(0, color=INK2, lw=0.6, ls=":")
        a.grid(axis="y", color="#e4e3df", lw=0.5)
        a.set_xlim(1, NPOS)
        if i == len(species) - 1:
            a.set_xlabel("position from read end")
    a = ax[i, 3]
    for k, period in enumerate(("historical", "modern")):
        v = [c[5] for c in sc if c[1] == period]
        a.scatter([k] * len(v), v, s=28, color=COL[period], edgecolor="white", lw=1, zorder=3)
    a.axhline(0, color=INK2, lw=0.6, ls=":")
    a.set_xticks([0, 1], ["historical", "modern"])
    a.set_xlim(-0.6, 1.6)
    a.grid(axis="y", color="#e4e3df", lw=0.5)
    ax[i, 0].set_ylabel(sp, color=INK, fontstyle="italic", fontsize=9)
for j, t in enumerate(titles):
    ax[0, j].set_title(t, color=INK, fontsize=9, loc="left")
handles = [plt.Line2D([], [], color=COL[p], lw=2, label=LABEL[p]) for p in ("historical", "modern")]
fig.legend(handles=handles, loc="upper right", frameon=False, ncol=2)
fig.suptitle("Pilot damage profiles, bwa mem mapping, 1M fragments per sample", x=0.01, ha="left",
             color=INK, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
for ext in ("pdf", "png"):
    fig.savefig(f"{B}/pilot/pilot_damage_composite.{ext}", dpi=150)

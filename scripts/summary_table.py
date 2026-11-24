"""Master summary: modern vs historical diversity for every analysis, all four species.

%change = (mean modern - mean historical) / mean historical * 100 (as in Liu et al. 2026).
Writes results/summary_table.tsv (long format: one row per analysis x species) and
results/summary_table.md (wide: analyses x species, '%change (n hist/n mod)').
Rerun any time; each cell uses whatever outputs exist and reports its sample sizes.
"""
import glob, gzip, os, re
import pandas as pd

B = "/path/to/TemporalGenomicsEvaluation"
SP = ["Actinodura_cyanouroptera", "Fringilla_montifringilla", "Trochalopteron_erythrocephalum", "Yuhina_occipitalis"]
short = {s: s.split("_")[0] for s in SP}
rows = []


def add(section, analysis, sp, hist, mod, note="", kind="pct"):
    """hist/mod: lists of per-sample values (or a single population value)."""
    hist, mod = list(hist), list(mod)
    if not hist or not mod:
        rows.append(dict(section=section, analysis=analysis, species=sp, n_hist=len(hist), n_mod=len(mod),
                         historical=None, modern=None, change=None, note=note, kind=kind)); return
    h, m = sum(hist) / len(hist), sum(mod) / len(mod)
    ch = 100 * (m - h) / h if kind == "pct" else m - h
    rows.append(dict(section=section, analysis=analysis, species=sp, n_hist=len(hist), n_mod=len(mod),
                     historical=h, modern=m, change=ch, note=note, kind=kind))


man = pd.read_csv(f"{B}/metadata/full_map_samples.tsv", sep="\t")
man["sp"] = man.species.str.replace(" ", "_"); man["sample"] = man.voucher + "_" + man.CRR
period = dict(zip(man["sample"], man.period))

# 0. paper (Table S1)
p = pd.read_csv(f"{B}/metadata/Liu2026_TableS1_h.tsv", sep="\t"); p["sp"] = p.species.str.replace(" ", "_")
for sp in SP:
    q = p[p.sp == sp]
    add("A. Per-individual heterozygosity", "Liu et al. 2026 (Table S1; ANGSD, as published)", sp,
        q[q.period == "historical"].h_paper, q[q.period == "modern"].h_paper)


# 1. ANGSD per-individual (native depth; depth-matched; rescaled)
def het_native(v):
    d = {}
    for f in glob.glob(f"{B}/het/{v}/*/*.het"):
        x = open(f).read().split("\t"); d[x[1]] = (x[0].replace(" ", "_"), x[2], float(x[-1]))
    return d


def het_matched(inp, v):
    d = {}
    for f in glob.glob(f"{B}/het_matched/{inp}/{v}/*/*.het"):
        x = open(f).read().split("\t"); d[x[1]] = (x[0].replace(" ", "_"), x[2], float(x[9]))
    return d


lab = {"all": "all sites", "tv": "transversions only"}
for v in ("all", "tv"):
    nat, mr, rs = het_native(v), het_matched("raw", v), het_matched("rescaled", v)
    for sp in SP:
        g = lambda d, per: [h for s, (spp, pp, h) in d.items() if spp == sp and pp == per]
        add("A. Per-individual heterozygosity", f"ANGSD GL, native depth, {lab[v]}", sp, g(nat, "historical"), g(nat, "modern"))
        add("A. Per-individual heterozygosity", f"ANGSD GL, modern downsampled to paired historical depth, {lab[v]}", sp,
            g(nat, "historical"), g(mr, "modern"))
        add("A. Per-individual heterozygosity", f"ANGSD GL, mapDamage-rescaled + depth-matched, {lab[v]}", sp,
            g(rs, "historical"), g(rs, "modern"))

# 2. HEARTY (5 reads per site, het if minor allele on >=2 reads; t = 0.25)
pat = re.compile(r": (\d+) sites, ([0-9.]+) sites with")


def hearty_H(f):
    n = h = 0.0
    for l in gzip.open(f, "rt"):
        m = pat.search(l); n += int(m.group(1)); h += float(m.group(2))
    hv = sum(float(pat.search(l).group(2)) for l in gzip.open(f.replace("windows.txt", "windows_trv.txt"), "rt"))
    return h / n, hv / n


def hearty(inp, tag_glob, folder="ds5_*"):
    d = {}
    for f in glob.glob(f"{B}/hearty/{inp}/*/{folder}/*_het0.25.mincov5_*{tag_glob}.windows.txt.gz"):
        s = os.path.basename(f).split("_het0.25")[0]
        if s in d or not os.path.exists(f.replace("windows.txt", "windows_trv.txt")): continue
        d[s] = (f.split("/")[-3], period[s], *hearty_H(f))
    return d


hx = {inp: hearty(inp, "downsampled_depth5_exact") for inp in ("raw", "rescaled")}
h1 = {inp: hearty(inp, "downsampled_depth5_reps1", "ds5_r1_s1") for inp in ("raw", "rescaled")}  # old single draw, seed 1
h10 = hearty("raw", "downsampled_depth5_reps10")
for inp in ("raw", "rescaled"):
    comb = dict(h1[inp]); comb.update(hx[inp])  # exact where available, else one draw (same expectation)
    n_fill = sum(1 for s in comb if s not in hx[inp])
    for i, v in ((2, "all"), (3, "tv")):
        for sp in SP:
            g = lambda per: [x[i] for s, x in comb.items() if x[0] == sp and x[1] == per]
            nf = sum(1 for s, x in comb.items() if x[0] == sp and s not in hx[inp])
            add("A. Per-individual heterozygosity",
                f"HEARTY, 5 reads/site, >=2 minor reads, {'plain' if inp == 'raw' else 'rescaled'} BAMs, {lab[v]}", sp,
                g("historical"), g("modern"), note=f"{nf} birds from single-draw run (exact otherwise)" if nf else "")
for i, v in ((2, "all"), (3, "tv")):
    for sp in SP:
        g = lambda per: [x[i] for s, x in h10.items() if x[0] == sp and x[1] == per]
        add("A. Per-individual heterozygosity", f"HEARTY old default (mean allele frequency of 10 draws), {lab[v]}", sp,
            g("historical"), g("modern"), note="superseded; shown for comparison")

# 3. ATLAS
at = pd.read_csv(f"{B}/atlas/atlas_summary.tsv", sep="\t")
for lev, name in (("5", "all birds downsampled to 5x"), ("native", "native depth")):
    for sp in SP:
        q = at[(at.species == short[sp]) & (at.level.astype(str) == lev)]
        add("A. Per-individual heterozygosity", f"ATLAS (per-sample PMD + recalibration, GERP sites), {name}", sp,
            q[q.period == "historical"].het, q[q.period == "modern"].het)

# 4. population-level (authors' pipeline, autosomes)
for sp in SP:
    o = f"{B}/pop/{sp}/{sp}"
    H = {l.split()[1]: float(l.split()[3]) for per in ("historical", "modern") for l in open(f"{o}.{per}.H")}
    add("B. Population-level (authors' ANGSD pipeline, autosomes)", "Authors' 'heterozygosity' (singleton class of unfolded SFS)", sp,
        [H["historical"]], [H["modern"]])
    st = {}
    for per in ("historical", "modern"):
        w = t = n = dd = 0.0
        for l in list(open(f"{o}.{per}.Win_50K.pestPG"))[1:]:
            x = l.split("\t"); w += float(x[3]); t += float(x[4]); n += float(x[13]); dd += float(x[8]) * float(x[13])
        st[per] = (w / n, t / n, dd / n)
    add("B. Population-level (authors' ANGSD pipeline, autosomes)", "Watterson's theta per site", sp, [st["historical"][0]], [st["modern"][0]])
    add("B. Population-level (authors' ANGSD pipeline, autosomes)", "Nucleotide diversity (pi) per site", sp, [st["historical"][1]], [st["modern"][1]])
    add("B. Population-level (authors' ANGSD pipeline, autosomes)", "Tajima's D (mean of 50 kb windows; change = difference)", sp,
        [st["historical"][2]], [st["modern"][2]], kind="diff")
    fst = open(f"{o}.fst.txt").read().split()
    rows.append(dict(section="B. Population-level (authors' ANGSD pipeline, autosomes)", analysis="Fst historical vs modern (weighted)",
                     species=sp, n_hist=5, n_mod=5, historical=None, modern=None, change=float(fst[1]), note="value, not a change", kind="value"))

# 5. SCAR (modern Fringilla degraded; change relative to the same birds' unmodified reads, equal depth)
sd = f"{B}/scar_sim/Fringilla_montifringilla"
for v in ("all", "tv"):
    vals = {}
    for f in glob.glob(f"{sd}/*/het/{v}/*.het"):
        b = os.path.basename(f)[:-4]; s, arm = b.split("_CRR")[0], b.split("_", 3)[3]
        vals.setdefault(arm, {})[s] = float(open(f).read().split("\t")[-2])
    for arm, name in (("frag", "fragmentation (historical length distribution)"), ("frag_deam", "+ deamination"),
                      ("frag_deam_tvex", "+ deamination + transversion excess")):
        both = [s for s in vals.get(arm, {}) if s in vals.get("original", {})]
        add("C. SCAR: modern Fringilla degraded to look historical (ANGSD GL, equal depth ~4.7x)",
            f"{name}, {lab[v]} (change vs same birds unmodified)", "Fringilla_montifringilla",
            [vals["original"][s] for s in both], [vals[arm][s] for s in both], note="'historical' = unmodified, 'modern' = degraded")

# 6. SCAR with the damage-robust methods (same arms; HEARTY exact 5 reads/site on full arm BAMs,
#    ATLAS with all arms of a bird downsampled to a common depth <= 5x)
SCAR_ARMS = (("frag", "fragmentation (historical length distribution)"), ("frag_deam", "+ deamination"),
             ("frag_deam_tvex", "+ deamination + transversion excess"))
SECD = "D. SCAR: same degraded Fringilla, damage-robust methods"
hv = {}
for f in glob.glob(f"{sd}/*/hearty/*/*_het0.25*exact.windows.txt.gz"):
    s, arm = f.split("/")[-4].split("_CRR")[0], f.split("/")[-2]
    hv.setdefault(arm, {})[s] = hearty_H(f)
av = {}
for f in glob.glob(f"{sd}/*/atlas/*/heterozygosity_stats.txt.gz"):
    s, arm = f.split("/")[-4].split("_CRR")[0], f.split("/")[-2]
    t = pd.read_csv(f, sep="\t"); av.setdefault(arm, {})[s] = float(t.HKY85_het.iloc[0])
for arm, name in SCAR_ARMS:
    for k, v in ((0, "all"), (1, "tv")):
        both = [s for s in hv.get(arm, {}) if s in hv.get("original", {})]
        add(SECD, f"HEARTY: {name}, {lab[v]} (change vs same birds unmodified)", "Fringilla_montifringilla",
            [hv["original"][s][k] for s in both], [hv[arm][s][k] for s in both], note="'historical' = unmodified, 'modern' = degraded")
    both = [s for s in av.get(arm, {}) if s in av.get("original", {})]
    add(SECD, f"ATLAS: {name}, all sites (change vs same birds unmodified)", "Fringilla_montifringilla",
        [av["original"][s] for s in both], [av[arm][s] for s in both], note="'historical' = unmodified, 'modern' = degraded")

d = pd.DataFrame(rows)
os.makedirs(f"{B}/results", exist_ok=True)
d.to_csv(f"{B}/results/summary_table.tsv", sep="\t", index=False, float_format="%.6g")


def cell(r):
    if r.change is None or pd.isna(r.change): return "–"
    if r.kind == "value": return f"{r.change:.3f}"
    if r.kind == "diff": return f"{r.historical:.2f} → {r.modern:.2f}"
    s = f"{r.change:+.1f}%"
    return s if r.section.startswith("B") else f"{s} ({r.n_hist}/{r.n_mod})"


d["cell"] = d.apply(cell, axis=1)
with open(f"{B}/results/summary_table.md", "w") as o:
    for sec, g in d.groupby("section", sort=False):
        cols = [s for s in SP if s in set(g.species)]
        o.write(f"\n### {sec}\n\n| analysis | " + " | ".join(short[s] for s in cols) + " |\n|---|" + "---|" * len(cols) + "\n")
        for a, gg in g.groupby("analysis", sort=False):
            c = dict(zip(gg.species, gg.cell))
            o.write(f"| {a} | " + " | ".join(c.get(s, "–") for s in cols) + " |\n")
print(open(f"{B}/results/summary_table.md").read())

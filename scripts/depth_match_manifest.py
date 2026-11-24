"""Depth-matching manifest for ANGSD heterozygosity.

Within each species, historical and modern samples are sorted by mean autosomal depth (PlainMap
coverage table) and paired by rank; each modern BAM is subsampled to its historical partner's depth
(fraction 1 if the modern sample is already shallower). Historical samples keep their native depth.

Writes one row per BAM to analyse, for two inputs:
  raw      : PlainMap BAMs; only the downsampled modern samples (historical raw H already exists in het/)
  rescaled : mapDamage-rescaled BAMs; historical native + modern downsampled
columns: species sample period input bam ref native_depth fraction target_depth partner seed

usage: depth_match_manifest.py [species ...] > manifest.tsv
"""
import os
import sys

B = "/path/to/TemporalGenomicsEvaluation"


def autosomal_depth(cov, ref):
    rf = f"{B}/references/regions/{os.path.basename(ref)[:-6]}.autosomes.txt"
    auto = {l.strip().rstrip(":") for l in open(rf) if l.strip()}
    L = D = 0.0
    for line in open(cov):
        if line.startswith("#"):
            continue
        x = line.split("\t")
        if x[0] in auto:
            L += float(x[2]); D += float(x[2]) * float(x[6])
    return D / L


want = set(sys.argv[1:])
rows = [l.rstrip("\n").split("\t") for l in open(f"{B}/metadata/full_map_samples.tsv")][1:]
by = {}
for n, (species, voucher, crr, period, ref) in enumerate(rows, 1):
    if want and species.replace(" ", "_") not in want:
        continue
    sp = species.replace(" ", "_")
    prefix = f"{voucher}_{crr}"
    cov = f"{B}/mapped/modern/{sp}/{prefix}.plainmap.coverage.tsv"
    if not os.path.exists(cov):
        sys.exit(f"not mapped yet: {prefix}")
    by.setdefault(species, {}).setdefault(period, []).append((autosomal_depth(cov, ref), prefix, ref, n))

print("\t".join("species sample period input bam ref native_depth fraction target_depth partner seed".split()))
for species, d in by.items():
    sp = species.replace(" ", "_")
    hist = sorted(d["historical"]); mod = sorted(d["modern"])
    if len(hist) != len(mod):
        print(f"# {species}: {len(hist)} historical vs {len(mod)} modern; pairing the first min(n)", file=sys.stderr)
    for (hd, hp, ref, hn), (md, mp, _, mn) in zip(hist, mod):
        f = min(1.0, hd / md)
        for inp in ("raw", "rescaled"):
            bam = (f"{B}/mapped/modern/{sp}/{mp}.final.RG.bam" if inp == "raw"
                   else f"{B}/mapped/rescaled/{sp}/{mp}.rescaled.bam")
            print(f"{species}\t{mp}\tmodern\t{inp}\t{bam}\t{ref}\t{md:.4f}\t{f:.6f}\t{md * f:.4f}\t{hp}\t{2000 + mn}")
    for hd, hp, ref, hn in hist:
        print(f"{species}\t{hp}\thistorical\trescaled\t{B}/mapped/rescaled/{sp}/{hp}.rescaled.bam\t{ref}\t"
              f"{hd:.4f}\t1\t{hd:.4f}\t-\t0")

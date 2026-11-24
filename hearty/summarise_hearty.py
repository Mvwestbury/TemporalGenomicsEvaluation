"""Genome-wide HEARTY heterozygosity per BAM and setting (sum of 100 kb window counts) -> hearty/hearty_summary.tsv"""
import glob, gzip, re, os
B = "/path/to/TemporalGenomicsEvaluation"
man = {}
for l in list(open(f"{B}/hearty/manifest.tsv"))[1:]:
    sp, s, per, inp, *_ = l.rstrip("\n").split("\t"); man[(inp, s)] = (sp, per)
pat = re.compile(r": (\d+) sites, (\d+) sites with")
rows = []
for f in glob.glob(f"{B}/hearty/*/*/ds5_*/*_het0.25.*windows.txt.gz"):
    inp, sp, setting = f.split("/")[-4:-1]
    s = os.path.basename(f).split("_het0.25")[0]
    fv = f.replace("windows.txt.gz", "windows_trv.txt.gz")
    if not os.path.exists(fv): continue
    n = h = 0
    for l in gzip.open(f, "rt"):
        m = pat.search(l); n += int(m.group(1)); h += int(m.group(2))
    hv = sum(int(pat.search(l).group(2)) for l in gzip.open(fv, "rt"))
    rows.append((sp, s, man[(inp, s)][1], inp, setting, n, h / n, hv / n))
with open(f"{B}/hearty/hearty_summary.tsv", "w") as o:
    o.write("species\tsample\tperiod\tinput\tsetting\tsites\tH_all\tH_tv\n")
    for r in sorted(rows):
        o.write("\t".join(map(str, r[:6])) + f"\t{r[6]:.6g}\t{r[7]:.6g}\n")
print(len(rows), "rows")

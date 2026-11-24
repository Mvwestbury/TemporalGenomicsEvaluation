"""Build a SCAR damage profile holding only the excess damage over the read-interior background.

Deamination (from a mapDamage misincorporation.txt, summed over Chr and strand):
    5p C->T and 3p G->A at positions 1..N, minus the mean rate over interior positions
    (default 30-60 from that end), clamped at 0; position N+1 written as 0.
Transversion excess (optional, from an TvCount positions.tsv):
    G->T and C->A excess (excess_gt_rate / excess_ca_rate) from the 5p end at positions 1..N,
    clamped at 0; position N+1 = mean of positions 10..N (the flat interior plateau).

SCAR caps positions at the profile's maximum and tries both ends for every base, so
  - deamination rows end with an explicit 0 -> no damage in the read interior;
  - transversion excess is written on 5p rows only -> applied once per base, plateau along the rest of the read.

usage: scar_profile.py --misinc misincorporation.txt [--tvcount positions.tsv] [-n 25] > profile.txt
"""
import argparse
import sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--misinc", required=True)
ap.add_argument("--tvcount")
ap.add_argument("-n", type=int, default=25)
ap.add_argument("--interior", default="30-60")
a = ap.parse_args()
lo, hi = map(int, a.interior.split("-"))

num = defaultdict(float)    # (end, pos, sub) -> count
den = defaultdict(float)    # (end, pos, base) -> count
with open(a.misinc) as f:
    header = None
    for line in f:
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if header is None:
            header = x
            continue
        r = dict(zip(header, x))
        end, pos = r["End"], int(r["Pos"])
        for b in "CG":
            den[end, pos, b] += float(r[b])
        for s in ("C>T", "G>A"):
            num[end, pos, s] += float(r[s])


def rate(end, pos, sub):
    d = den[end, pos, sub[0]]
    return num[end, pos, sub] / d if d else 0.0


out = ["# SCAR profile: excess over interior background", "# end\tposition\tfrom\tto\tfrequency"]
for end, sub in (("5p", "C>T"), ("3p", "G>A")):
    base = sum(rate(end, p, sub) for p in range(lo, hi + 1)) / (hi - lo + 1)
    print(f"{end} {sub}: interior baseline {base:.5f}, pos1 {rate(end, 1, sub):.5f}, pos2 {rate(end, 2, sub):.5f}",
          file=sys.stderr)
    for p in range(1, a.n + 1):
        out.append(f"{end}\t{p}\t{sub[0]}\t{sub[2]}\t{max(0.0, rate(end, p, sub) - base):.6g}")
    out.append(f"{end}\t{a.n + 1}\t{sub[0]}\t{sub[2]}\t0")

if a.tvcount:
    import csv
    rows = [r for r in csv.DictReader(open(a.tvcount), delimiter="\t") if r["End"] == "5p"]
    by = {int(r["Pos"]): r for r in rows}
    for col, (fr, to) in (("excess_gt_rate", "GT"), ("excess_ca_rate", "CA")):
        vals = {p: max(0.0, float(by[p][col])) for p in range(1, a.n + 1) if p in by}
        plateau = sum(vals[p] for p in range(10, a.n + 1)) / (a.n - 9)
        print(f"transversion excess {fr}->{to}: pos1 {vals.get(1, 0):.5f}, plateau {plateau:.5f}", file=sys.stderr)
        for p in range(1, a.n + 1):
            out.append(f"5p\t{p}\t{fr}\t{to}\t{vals.get(p, 0.0):.6g}")
        out.append(f"5p\t{a.n + 1}\t{fr}\t{to}\t{plateau:.6g}")

print("\n".join(out))

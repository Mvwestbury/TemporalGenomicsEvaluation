"""Expected 1-draw HEARTY het call per site (exact, no RNG): P(call = HET) when DEPTH reads are drawn
without replacement from the observed counts, using HEARTY's classify_counts_for_threshold.
Also: mean of per-draw calls over R draws (HEARTY RNG), for comparison.
usage: expected_het.py basecall.txt.gz min_depth max_depth [chrom]"""
import sys, gzip, itertools
from functools import lru_cache
from math import comb
sys.path.insert(0, "/path/to/HEARTY"); import HEARTY as H
DEPTH, T = 5, 0.25
f, mind, maxd = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); chrom = sys.argv[4] if len(sys.argv) > 4 else None
draws = [k for k in itertools.product(range(DEPTH + 1), repeat=4) if sum(k) == DEPTH]

@lru_cache(maxsize=None)
def expected(counts):
    n = sum(counts); tot = comb(n, DEPTH); p_all = p_tv = 0.0
    for k in draws:
        if any(k[i] > counts[i] for i in range(4)): continue
        _, st = H.classify_counts_for_threshold(H.counts_to_freqs(list(k)), T)
        if st.startswith("HET"):
            p = 1.0
            for i in range(4): p *= comb(counts[i], k[i])
            p /= tot; p_all += p; p_tv += p if st == "HET_trv" else 0
    return p_all, p_tv

def call(counts, ch, pos, seed):
    s = H.sample_counts_without_replacement(counts, DEPTH, H.site_rng(seed, ch, pos))
    return H.classify_counts_for_threshold(H.counts_to_freqs(s), T)[1]

n = 0; e = [0.0, 0.0]; one = {s: [0, 0] for s in (1, 2, 3)}; avg10 = [0.0, 0.0]
fh = gzip.open(f, "rt"); fh.readline()
for l in fh:
    x = l.split("\t")
    if chrom and x[0] != chrom:
        if n: break
        continue
    d = int(x[2]); c = tuple(int(v) for v in x[3:7])
    if not (mind <= d <= maxd) or sum(c) < DEPTH: continue
    n += 1; pa, pt = expected(c); e[0] += pa; e[1] += pt
    if pa == 0: continue  # no draw can be HET -> all methods give 0
    for s in (1, 2, 3):
        st = call(list(c), x[0], int(x[1]), s); one[s][0] += st.startswith("HET"); one[s][1] += st == "HET_trv"
    rng = H.site_rng(1, x[0], int(x[1]))
    for _ in range(10):
        st = H.classify_counts_for_threshold(H.counts_to_freqs(H.sample_counts_without_replacement(list(c), DEPTH, rng)), T)[1]
        avg10[0] += st.startswith("HET") / 10; avg10[1] += (st == "HET_trv") / 10
print(f"sites {n}")
for s in (1, 2, 3): print(f"1 draw seed {s:<11} H_all {one[s][0]/n:.6f}  H_tv {one[s][1]/n:.6f}")
print(f"10 draws, mean of calls H_all {avg10[0]/n:.6f}  H_tv {avg10[1]/n:.6f}")
print(f"exact expectation     H_all {e[0]/n:.6f}  H_tv {e[1]/n:.6f}")

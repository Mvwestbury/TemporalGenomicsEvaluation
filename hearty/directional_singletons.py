"""Reference-directional singleton/doubleton spectrum from a HEARTY basecall table.

Each site passing HEARTY's filters (min..max depth, >= DEPTH counted bases) is downsampled to DEPTH
reads with HEARTY's own per-site RNG (seed, one draw), so the draws are the ones HEARTY used for the
ds5_r1_s<seed> summaries. Non-reference reads among the DEPTH are counted:
  k = 1 -> singleton of type Ref>Alt;  k = 2 with both reads the same base -> doubleton Ref>Alt.
Output (tsv): sites (passing, Ref in ACGT), then one row per k and substitution type.

usage: directional_singletons.py basecall.txt.gz out.tsv min_depth max_depth [depth=5] [seed=1] [workers=6]
"""
import gzip
import sys
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, "/path/to/HEARTY")
import HEARTY as H  # noqa: E402

BASES = "ACGT"


def work(args):
    rows, mind, maxd, depth, seed = args
    c = Counter(); n = 0
    for chrom, site, d, counts, ref in rows:
        if not (mind <= d <= maxd) or sum(counts) < depth or ref not in BASES:
            continue
        n += 1
        s = H.sample_counts_without_replacement(counts, depth, H.site_rng(seed, chrom, site))
        r = BASES.index(ref)
        k = depth - s[r]
        if k == 1 or k == 2:
            alts = [i for i in range(4) if i != r and s[i]]
            if len(alts) == 1:
                c[(k, f"{ref}>{BASES[alts[0]]}")] += 1
    return n, c


def chunks(path, mind, maxd, depth, seed, size=200_000):
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        ix = [h.index(x) for x in ("chr", "pos", "totDepth", "count_A", "count_C", "count_G", "count_T", "Ref")]
        rows = []
        for line in f:
            x = line.rstrip("\n").split("\t")
            rows.append((x[ix[0]], int(x[ix[1]]), int(x[ix[2]]), [int(x[i]) for i in ix[3:7]], x[ix[7]]))
            if len(rows) == size:
                yield rows, mind, maxd, depth, seed
                rows = []
        if rows:
            yield rows, mind, maxd, depth, seed


def main():
    a = sys.argv[1:]
    path, out, mind, maxd = a[0], a[1], int(a[2]), int(a[3])
    depth = int(a[4]) if len(a) > 4 else 5
    seed = int(a[5]) if len(a) > 5 else 1
    workers = int(a[6]) if len(a) > 6 else 6
    total = 0; spec = Counter()
    with Pool(workers) as p:
        for n, c in p.imap_unordered(work, chunks(path, mind, maxd, depth, seed), chunksize=1):
            total += n; spec.update(c)
    with open(out + ".tmp", "w") as o:
        o.write(f"#sites\t{total}\n")
        o.write("k\ttype\tcount\trate\n")
        for k in (1, 2):
            for r in BASES:
                for b in BASES:
                    if r != b:
                        n = spec[(k, f"{r}>{b}")]
                        o.write(f"{k}\t{r}>{b}\t{n}\t{n / total:.6g}\n")
    import os
    os.replace(out + ".tmp", out)


if __name__ == "__main__":
    main()

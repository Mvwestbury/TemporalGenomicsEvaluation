"""Where do the mismatches against the reference come from?

For the first N_READS MAPQ>=20 primary reads on autosomes of each BAM, tabulates mismatches
(from the MD tag) by base quality, by sequencing cycle and by substitution type, plus indels.
Usage: mismatch_diagnostics.py OUT_PREFIX BAM [BAM ...]
"""
import collections, re, sys
import pysam

N_READS = 200000
AUTO = re.compile(r"^(Pseudochr\d+|.*_chr\d+)$")
QBINS = [(0, 20), (20, 30), (30, 36), (36, 40), (40, 99)]
COMP = str.maketrans("ACGT", "TGCA")

out = sys.argv[1]
summ, byq, bycyc, spec = [], [], [], []
for path in sys.argv[2:]:
    name = path.split("/")[-1].replace(".final.RG.bam", "")
    bam = pysam.AlignmentFile(path)
    qb = collections.Counter(); qm = collections.Counter()
    cb = collections.Counter(); cm = collections.Counter()
    sub = collections.Counter(); n = bases = mm = ins = dele = clipped = 0
    for chrom in (c for c in bam.references if AUTO.match(c)):
        for r in bam.fetch(chrom):
            if r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
                continue
            q = r.query_qualities; seq = r.query_sequence; L = r.query_length
            for op, l in r.cigartuples:
                ins += l if op == 1 else 0; dele += l if op == 2 else 0; clipped += l if op == 4 else 0
            for qp, rp, ref in r.get_aligned_pairs(matches_only=True, with_seq=True):
                cyc = L - qp if r.is_reverse else qp + 1
                k = next(i for i, (a, b) in enumerate(QBINS) if a <= q[qp] < b)
                bases += 1; qb[k] += 1; cb[cyc] += 1
                if ref.islower():
                    mm += 1; qm[k] += 1; cm[cyc] += 1
                    s = f"{ref.upper()}>{seq[qp]}"
                    if r.is_reverse:     # report on the read strand
                        s = s.translate(COMP)
                    sub[s] += 1
            n += 1
            if n >= N_READS: break
        if n >= N_READS: break
    summ.append([name, n, bases, mm / bases, ins / bases, dele / bases, clipped / (bases + clipped),
                 sum(qb[k] for k in (0, 1)) / bases, sum(qm[k] for k in (0, 1)) / mm,
                 (mm - qm[0] - qm[1]) / (bases - qb[0] - qb[1])])
    for k, (a, b) in enumerate(QBINS):
        byq.append([name, f"{a}-{b - 1}", qb[k], qm[k] / qb[k] if qb[k] else "NA"])
    for c in sorted(cb):
        bycyc.append([name, c, cb[c], cm[c] / cb[c]])
    for s in sorted(sub):
        spec.append([name, s, sub[s], sub[s] / mm])

def write(fn, header, rows):
    with open(f"{out}.{fn}.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(f"{x:.5g}" if isinstance(x, float) else str(x) for x in r) + "\n")

write("summary", ["sample", "reads", "aligned_bases", "mismatch_rate", "ins_rate", "del_rate", "softclip_frac",
                  "frac_bases_Qlt30", "frac_mismatches_at_Qlt30", "mismatch_rate_Q30plus"], summ)
write("by_quality", ["sample", "qual_bin", "bases", "mismatch_rate"], byq)
write("by_cycle", ["sample", "cycle", "bases", "mismatch_rate"], bycyc)
write("spectrum", ["sample", "substitution", "count", "fraction_of_mismatches"], spec)

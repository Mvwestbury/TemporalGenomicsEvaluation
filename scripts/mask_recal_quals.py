"""Keep ATLAS-recalibrated base qualities only where the ORIGINAL quality passed a cutoff.

Bases whose original quality was below --minq get quality 0 in the output, so that `angsd -minQ 1`
uses exactly the same bases as the original `-minQ 30` runs, but with ATLAS-recalibrated qualities.

usage: mask_recal_quals.py original.bam recal.bam out.bam [--minq 30]
Both BAMs must be coordinate-sorted; recal.bam (from `atlas filterBAM`) is a subset of original.bam.
"""
import argparse
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("original"); ap.add_argument("recal"); ap.add_argument("out")
ap.add_argument("--minq", type=int, default=30)
a = ap.parse_args()

orig = pysam.AlignmentFile(a.original)
rec = pysam.AlignmentFile(a.recal)
out = pysam.AlignmentFile(a.out, "wb", template=rec)
oit = orig.fetch(until_eof=True)
buf, cur = {}, None   # originals at positions <= current recal position, keyed by (tid, pos, name, flag)
n = masked = total = 0


def key(r): return (r.reference_id, r.reference_start, r.query_name, r.flag & 0x10)


for r in rec.fetch(until_eof=True):
    k = key(r)
    while k not in buf:
        o = next(oit)  # StopIteration here means the recal read is missing from the original: fail loudly
        if o.is_secondary or o.is_supplementary: continue
        buf[key(o)] = o
        if len(buf) > 200000:  # drop originals that are behind the current position (filtered by ATLAS)
            buf = {kk: v for kk, v in buf.items() if (kk[0], kk[1]) >= (k[0], k[1])}
    o = buf.pop(k)
    if o.query_sequence != r.query_sequence:
        raise SystemExit(f"sequence mismatch for {r.query_name}")
    oq, rq = o.query_qualities, list(r.query_qualities)
    for i, q in enumerate(oq):
        if q < a.minq: rq[i] = 0; masked += 1
    total += len(oq)
    r.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in rq))
    out.write(r); n += 1
out.close()
print(f"reads {n}, bases {total}, masked (original Q<{a.minq}) {masked} ({masked / max(total, 1):.4f})")

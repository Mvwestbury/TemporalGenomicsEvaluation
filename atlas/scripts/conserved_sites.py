"""(Copied from Black_bears ATLAS pipeline; lifted mode caps --count at the eligible sites.)
Rank published GERP scores before and after ATLAS position-based liftover."""
import argparse
import collections
import json
from pathlib import Path

import numpy as np


def cutoff(histogram, count):
    above = 0
    for score in sorted(histogram, reverse=True):
        if above + histogram[score] >= count:
            return score, count - above
        above += histogram[score]
    raise ValueError(f"Only {above} eligible sites, require {count}")


def select_bw(args):
    import pyBigWig
    bw = pyBigWig.open(str(args.input))
    histogram = collections.Counter()
    chromosomes = bw.chroms()

    def chunks():
        for chrom, length in chromosomes.items():
            for start in range(0, length, 2_000_000):
                values = bw.values(chrom, start, min(start + 2_000_000, length), numpy=True)
                yield chrom, start, values

    for chrom, start, values in chunks():
        scores, counts = np.unique(values[np.isfinite(values) & (values > 0)], return_counts=True)
        histogram.update(dict(zip(scores.tolist(), counts.tolist())))
    threshold, ties_left = cutoff(histogram, args.count)
    selected = 0
    with args.output.open("w") as out:
        for chrom, start, values in chunks():
            mask = np.isfinite(values) & (values > threshold)
            ties = np.flatnonzero(values == threshold)[:ties_left]
            mask[ties] = True
            ties_left -= len(ties)
            positions = np.flatnonzero(mask)
            if not len(positions):
                continue
            # Equal-score adjacent bases may be one BED interval; ATLAS expands them.
            boundaries = np.r_[0, np.flatnonzero((np.diff(positions) != 1) |
                               (np.diff(values[positions]) != 0)) + 1, len(positions)]
            for left, right in zip(boundaries[:-1], boundaries[1:]):
                first, last = positions[left], positions[right - 1]
                out.write(f"{chrom}\t{start+first}\t{start+last+1}\t{float(values[first]):.9g}\n")
            selected += len(positions)
    assert selected == args.count, (selected, args.count)
    args.output.with_suffix(".json").write_text(json.dumps({"sites": selected, "min_GERP": threshold}, indent=2) + "\n")
    bw.close()


def select_lifted(args):
    allowed = set(args.autosomes.read_text().split())
    histogram = collections.Counter()

    def rows():
        with args.input.open() as inp:
            for line in inp:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end, score = line.split()
                if chrom in allowed:
                    yield chrom, int(start), int(end), float(score)

    for _, start, end, score in rows():
        histogram[score] += end - start
    # birds: fewer sites may survive liftover from a distant source genome; take what exists
    args.count = min(args.count, sum(histogram.values()))
    threshold, ties_left = cutoff(histogram, args.count)
    selected = 0
    previous = None
    with args.output.open("w") as out:
        for chrom, start, end, score in rows():
            if score < threshold:
                continue
            if score == threshold:
                end = min(end, start + ties_left)
                ties_left -= end - start
            if start == end:
                continue
            selected += end - start
            if previous and previous[0] == chrom and start <= previous[2]:
                if start < previous[2]:
                    raise ValueError("Overlapping liftover positions")
                previous = (chrom, previous[1], end)
            else:
                if previous:
                    out.write("\t".join(map(str, previous)) + "\n")
                previous = (chrom, start, end)
        if previous:
            out.write("\t".join(map(str, previous)) + "\n")
    assert selected == args.count, (selected, args.count)
    args.output.with_suffix(".json").write_text(json.dumps({"sites": selected, "min_GERP": threshold,
        "eligible_lifted_sites": sum(histogram.values())}, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["source", "lifted"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--count", type=int, required=True)
    p.add_argument("--autosomes", type=Path)
    args = p.parse_args()
    if args.mode == "source":
        select_bw(args)
    else:
        select_lifted(args)


if __name__ == "__main__":
    main()

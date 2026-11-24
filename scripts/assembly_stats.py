"""Assembly size, scaffold N50 and contig N50 (scaffolds split at runs of >=10 N) for FASTA files."""
import re, sys

def n50(lengths):
    lengths = sorted(lengths, reverse=True); half = sum(lengths) / 2; s = 0
    for l in lengths:
        s += l
        if s >= half: return l

print("assembly\tsequences\ttotal_bp\tN_bp\tscaffold_N50\tcontig_N50\tlargest_contig\tbp_in_Pseudochr_autosomes")
for path in sys.argv[1:]:
    scaf, ctg, nbp, auto = [], [], 0, 0
    def add(name, seq):
        global nbp, auto
        seq = "".join(seq).upper(); scaf.append(len(seq)); nbp += seq.count("N")
        if re.fullmatch(r"(Pseudochr\d+|.*_chr\d+)", name): auto += len(seq)
        ctg.extend(len(c) for c in re.split(r"N{10,}", seq) if c)
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name: add(name, seq)
            name, seq = line[1:].split()[0], []
        else: seq.append(line.strip())
    add(name, seq)
    print(f"{path.split('/')[-1]}\t{len(scaf)}\t{sum(scaf)}\t{nbp}\t{n50(scaf)}\t{n50(ctg)}\t{max(ctg)}\t{auto}")

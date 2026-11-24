# TemporalGenomicsEvaluation

Historical-heterozygosity bias: reanalysis of Liu et al. (2026) and the Chatham Island black robin.

Code to reproduce the analyses of the paper. The repository contains:
- every cluster script (SLURM `sbatch`) and analysis script used for the paper;
- the per-sample damage tables;
- the exact inputs given to the read simulator SCAR;
- the plotted values of every figure.

The scripts are the versions that were run. Three things were changed: absolute cluster paths are
replaced by placeholders (see **Paths**), the third SCAR arm and its profiles are named `tvex`
(transversion excess), and the in-house damage counter is called TvCount, consistently throughout. Reads, BAMs and other large intermediate files are not included.

## Data

| Dataset | Accessions | Sample tables |
|---|---|---|
| Liu et al. 2026: *Fringilla montifringilla*, *Actinodura cyanouroptera*, *Trochalopteron erythrocephalum*, *Yuhina occipitalis* (5 historical + 5 modern each) | GSA CRA020079, CRA032748 | `metadata/full_map_samples.tsv`, `metadata/pilot_samples.tsv`, `metadata/four_species_*.tsv` |
| Black robin, *Petroica traversi* (5 historical + 5 modern) | ENA, runs in `robin/robin_files.tsv` | `robin/robin_samples.tsv` |
| References | bFriMon1 (GCA_988331165); assemblies of Liu et al. (MT275, GLG_GLGS1218, GLG_GLGS5798); *P. traversi* NRM v1 | `prep_refs.sbatch`, `robin/rb_ref.sbatch` |

## Layout (the cluster project layout, kept so that the scripts' relative paths are valid)

| Path | What it does |
|---|---|
| `download_reads.*`, `fetch_chunked.*`, `wget_files.*`, `robin/rb_reads.sbatch` | Download and MD5-check the reads |
| `prep_refs.sbatch`, `actinodura_*.sbatch`, `robin/rb_ref.sbatch` | Index the references; check the *Actinodura* reference (read mismatch rates against candidate assemblies) |
| `pilot_map.sbatch` | 1 M-fragment pilot mapping, mapDamage and TvCount (both bwa mem and bwa aln modes) |
| `full_map.sbatch`, `rescale_md.sbatch`, `robin/rb_map.sbatch`, `robin/rb_rescale_md.sbatch` | Full mapping with PlainMap (bwa mem) and mapDamage quality rescaling |
| `het_ind.sbatch`, `het_matched.sbatch`, `scripts/depth_match_manifest.py`, `robin/rb_het_ind.sbatch` | ANGSD per-individual heterozygosity (all sites and transversions; native and depth-matched) |
| `pop_angsd.sbatch` | Replication of the population-level ANGSD analysis of Liu et al. |
| `hearty/` | HEARTY heterozygosity; `directional_singletons.py` gives singleton/doubleton counts by substitution class |
| `atlas/` | ATLAS heterozygosity (calibration sites from GERP, depth series native/15/10/7.5/5×) |
| `bcftools_roh/` | bcftools ROH in modern birds |
| `scar_*.sbatch`, `robin/rb_profile.sbatch`, `robin/rb_scar_*.sbatch`, `scripts/scar_profile.py` | SCAR degradation simulations (see below) |
| `robin/submit_all.sh` | Full black-robin chain with SLURM dependencies |
| `scripts/` | Summary tables and figures (`summary_table.py`, `figures.py`, `figure_liu_robin.py`, `main_figures.py`, `singleton_doubleton.py`, `pilot_damage_composite.py`) |
| `results/figures/*.tsv` | Plotted values of every figure |
| `damage/` | Reduced TvCount tables (see below) |
| `scar_sim/profiles/` | SCAR inputs and their sources (see below) |
| `*/manifest*.tsv`, `atlas/atlas_summary.tsv`, `hearty/*summary.tsv`, `bcftools_roh/bcftools_froh_summary.tsv`, `pilot/*.tsv` | Per-sample manifests and per-method summaries |

## Paths

Cluster paths are replaced by placeholders. To run the scripts, point them at your installation:

```bash
grep -rl '/path/to' . | xargs sed -i -e 's#/path/to/TemporalGenomicsEvaluation#'"$PWD"'#g' \
    -e 's#/path/to/software#/your/software#g' -e 's#/path/to/HEARTY#/your/HEARTY#g'
```

`/path/to/GT_damage` (in `robin/rb_ref.sbatch` and `robin/rb_reads.sbatch`) is where the robin reference and
historical reads were stored locally. `scripts/safe_mkdir.sh` works around a cluster-specific
`mkdir -p` race and can be removed elsewhere.

## Software

| Tool | Version |
|---|---|
| PlainMap (mapping wrapper: fastp, bwa, samtools, mapDamage) | commit 5d0126f |
| bwa | 0.7.17-r1188 |
| samtools / bcftools | 1.24 / 1.22 |
| fastp | 1.0.1 |
| mapDamage | 2.2.3 |
| ANGSD | 0.935 |
| HEARTY | commit d8305fe |
| ATLAS | 2.0.4, commit a37654b9 |
| SCAR | commit 38719e8 (unpublished, in prep.; see below) |
| TvCount | unpublished (see below) |
| Python | 3.13 (pandas 3.0, numpy 2.4, matplotlib 3.10) |

## Damage tables (TvCount)

TvCount is not yet published. The calls to it (`tools/TvCount/TvCount.py` in `pilot_map.sbatch`, `scar_map.sbatch`,
`robin/rb_damage.sbatch` and `robin/rb_profile.sbatch`) are left in so the settings are visible. Its outputs are
released as count tables, so every rate can be recomputed without the software.

For each aligned read (1 M reads per sample; MAPQ and duplicate filters as in the scripts), TvCount counts
reference G and C bases and the substitutions G>T, G>C, C>A and C>G by distance from the 5′ and 3′ read ends.
G>C and C>G are strand-symmetric transversions that damage does not produce, so they serve as the
background control for G>T and C>A.

`damage/tvcount_positions.tsv`: one row per sample × read end (`End`, 5p or 3p) × position (`Pos`, 1 = terminal base).

| Column | Definition |
|---|---|
| `G`, `C` | reference G / C bases observed at that position |
| `G>T`, `G>C`, `C>A`, `C>G` | read bases differing from the reference |
| `excess_gt_rate` | `G>T/G − G>C/G` |
| `excess_ca_rate` | `C>A/C − C>G/C` |

`damage/tvcount_summary.tsv`: one row per sample, with genome-wide counts over all positions, the four rates
(`overall_p_gt_over_g` = gt_count / g_bases_observed, and so on), and
`overall_combined_excess_rate` = ((gt − gc)/G + (ca − cg)/C) / 2.

Both tables share the metadata columns `dataset`, `species`, `sample`, `period`, `mapping` and `arm`:

| `dataset` | What it holds |
|---|---|
| `liu_pilot_1M` | Liu birds from the pilot; `mapping` = bwa_mem or bwa_aln |
| `robin_1M` | Robins |
| `scar_1M` | Simulated SCAR arms (`arm` = original, frag, frag_deam, frag_deam_tvex); this checks that the simulated damage reached the reads |

A second line of evidence for the transversion excess needs no TvCount at all: the singleton/doubleton spectrum
from the HEARTY 5-read draws (`hearty/directional_singletons.py`, `scripts/singleton_doubleton.py`, plotted values in
`results/figures/figS_singleton_doubleton_data.tsv`).

## SCAR simulations

SCAR (commit 38719e8) is not yet public. A simulation is fully specified by its inputs, so these are released
exactly as given to SCAR:

**Inputs (`scar_sim/profiles/`)**

| File | Content |
|---|---|
| `Fringilla_worst_deam.txt`, `Fringilla_worst_deam_tvex.txt` | Damage profiles used for the three modern *Fringilla* in `scar_sim/fringilla_samples.tsv` |
| `Robin_BRH08_deam.txt`, `Robin_BRH08_deam_tvex.txt` | Damage profiles used for the five modern robins in `robin/robin_scar_samples.tsv` |
| `KIZ_013086_CRR2269149.fraglen.txt.gz`, `BRH08_ERR6335811.fraglen.txt.gz` | Fragment lengths (one per line) sampled from a historical bird; gunzip before use |

Profile format: `end  position  from  to  frequency`. The frequency is the per-base probability of that
substitution in excess of the read-interior background. Positions run 1–25 from the read end; the last row
(26) applies to the rest of the read. The profiles contain the following rows:
- 5′ C→T and 3′ G→A: deamination, from mapDamage `misincorporation.txt` minus the interior mean (positions 30–60), clamped at 0, interior 0;
- `_tvex` files add 5′ G→T and C→A: the TvCount `excess_gt_rate` / `excess_ca_rate`, clamped at 0, with the interior row set to the mean of positions 10–25.

**Sources (`scar_sim/profiles/sources/`).** These are the exact tables the profiles were built from with `scripts/scar_profile.py`:

| Profile | Deamination from | Transversion excess from |
|---|---|---|
| Fringilla_worst | mapDamage, KIZ_013087 (pilot, bwa mem) | TvCount, KIZ_013083 (pilot, bwa mem; highest excess of the historical *Fringilla*) |
| Robin_BRH08 | mapDamage, BRH08 (1 M reads) | TvCount, BRH08 (1 M reads) |

To regenerate the four profiles byte-identically:

```bash
S=scar_sim/profiles/sources
python scripts/scar_profile.py --misinc $S/KIZ_013087_CRR2269153.misincorporation.txt > Fringilla_worst_deam.txt
python scripts/scar_profile.py --misinc $S/KIZ_013087_CRR2269153.misincorporation.txt \
    --tvcount $S/KIZ_013083_CRR2269150.tvcount_positions.tsv > Fringilla_worst_deam_tvex.txt
python scripts/scar_profile.py --misinc $S/BRH08_ERR6335811.misincorporation.txt > Robin_BRH08_deam.txt
python scripts/scar_profile.py --misinc $S/BRH08_ERR6335811.misincorporation.txt \
    --tvcount $S/BRH08_ERR6335811.tvcount_positions.tsv > Robin_BRH08_deam_tvex.txt
```

**Fragment lengths.** These are MAPQ ≥ 30 primary reads from a 0.5% sample of the historical BAM: single-end read length,
or |TLEN| of the first mate for proper pairs. The command is in `robin/rb_profile.sbatch`; the *Fringilla* file
was made with the same command on KIZ_013086.

**Commands** (`scar_prep.sbatch`). The modern reads are first trimmed and merged with fastp (`-l 30`) into one
single-end file (`original`). Then:

```bash
scar -i original.fq.gz -o frag         --fd empirical --fdf fraglen.txt --mfl 30 --mutation-rate 0 --seed 11
scar -i frag.fastq.gz  -o frag_deam    -d profile_deam.txt    --seed 12
scar -i frag.fastq.gz  -o frag_deam_tvex -d profile_deam_tvex.txt --seed 12
```

Both damage arms start from the same fragmented reads, so the arms differ only in the damage added.
Every arm, including `original`, is then mapped and analysed exactly like the empirical data (`scar_map.sbatch`,
`scar_het.sbatch`, `scar_hearty.sbatch`, `scar_atlas.sbatch`; robin copies in `robin/`).


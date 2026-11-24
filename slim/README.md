# Metapopulation simulations of genomic erosion

This directory contains the forward-time simulations used to investigate the expected loss of genetic diversity following demographic collapse and reduced connectivity in a metapopulation. The model is implemented in SLiM 4 using a diploid Wright–Fisher model with five demes, non-overlapping generations and symmetric migration.

Following an ancestral burn-in, all five demes simultaneously change in population size and migration rate. Genomic metrics are recorded at the end of the burn-in (historical sample), during the bottleneck and at the end of the bottleneck (modern sample). The model uses the genomic architecture of chromosome 12 of the collared flycatcher (*Ficedula albicollis*), including variation in recombination rates along the chromosome.

Question can be send to Hernán E. Morales Herman [@] sund.ku.dk

## Files

| File | Description |
| --- | --- |
| `metapop_5demes_WF_with_FST.slim` | Main SLiM simulation script, including demographic scenarios, genomic architecture, sampling and calculation of genomic metrics. |
| `run_metapop_grid_FST.sbatch` | SLURM job-array script for running combinations of demographic parameters with independent replicates. |
| `coords_SLiM_GCF_000247815.1_FicAlb1.5_chr12.txt` | Genomic coordinates of exons, introns and non-coding regions for chromosome 12 of the collared flycatcher reference assembly (FicAlb1.5). |
| `out_rho4Ne_win2e+05_col.hun.snpPair_runALL_chr12.txt` | Chromosome 12 recombination map based on 200-kb windows. |

## Running the simulations

Requires **SLiM 4** and, for parallel execution on an HPC cluster, **SLURM**. Place all four files in the same directory. Simulation parameters can be supplied using the `-d` option in SLiM; additional settings are documented in the main script.

Example of a short neutral test run, using the default demographic parameters and a 100-generation burn-in (not intended to establish equilibrium):

```bash
mkdir -p out
slim -s 12345 \
    -d 'basePath="."' \
    -d 'outDir="./out"' \
    -d 'runType=1' \
    -d 'burnLen=100' \
    metapop_5demes_WF_with_FST.slim
```

To run the full parameter grid, first edit `SLIM_SCRIPT` in the SLURM script to point to `metapop_5demes_WF_with_FST.slim` and set `SLIM_BIN` to the local SLiM executable (or `slim` if it is on the PATH). Adjust the requested wall time and memory to your cluster and simulation settings.

```bash
mkdir -p logs
sbatch run_metapop_grid_FST.sbatch
```

## Simulation parameters

The supplied SLURM script defines the following neutral simulation grid:

| Parameter | Values |
| --- | --- |
| Ancestral deme size (`N_anc1`–`N_anc5`) | 1,000 individuals per deme |
| Ancestral migration (`M_anc`) | 0.1 |
| Burn-in (`burnLen`) | 50,000 generations |
| Bottleneck deme size (`N_btl1`–`N_btl5`) | 10, 20, 50, 100, 1,000 |
| Bottleneck migration (`M_btl`) | 0, 0.001, 0.01, 0.1 |
| Bottleneck duration (`btl_gens`) | 10, 20 generations |
| Replicates | 20 per scenario |
| Sample size (`sampleN`) | 5 diploid individuals per deme per sampling point |
| Minimum ROH length (`rohCutoff`) | 500,000 bp (0.5 Mb) |

This gives 40 demographic scenarios and 800 independent simulations. Migration is expressed as the total probability of offspring being derived from parents in another deme, distributed equally among the four other demes (`M/4` per source deme). The ancestral size, migration and burn-in can be overridden in the SLURM script or through its environment variables.

## Outputs

Each simulation produces tab-separated files with a filename prefix identifying the demographic scenario and random seed. The SLURM script writes results to `tmp_out_mjolnir_fst/`.

| Output suffix | Content |
| --- | --- |
| `_parameters.tsv` | Simulation settings and random seed. |
| `_populations.tsv` | Deme sizes, migration and sampling generations. |
| `_diversity.tsv` | Sample-based and census-wide mean individual heterozygosity, Watterson's theta and Tajima's D (the latter two from the sample). |
| `_pi.tsv`, `_piTotal.tsv` | Sample mean individual heterozygosity (neutral and all mutations, respectively). |
| `_individual.tsv` | Individual-level genomic statistics. |
| `_roh.tsv` | Individual ROH counts, lengths, coordinates and FROH. |
| `_comparison.tsv` | Historical-to-modern heterozygosity loss and mean FROH. |
| `_fst.tsv` | Temporal Hudson's FST within each deme, calculated for sampled individuals and the entire deme. |

ROH are estimated from contiguous regions without heterozygous mutation markers on the original simulated chromosome, excluding the appended synthetic exome. These are proxies for homozygous tracts rather than direct identity-by-descent calls.

# Non-B DNA analysis pipeline

This Snakemake workflow intersects human transcript annotations with nine
Non-B DNA motif sets and compares protein-coding and lncRNA transcripts. The
default `all` target produces per-motif summaries and plots, an extended feature
matrix, contingency tests, statistical tests, and an analysis completion marker.

The configured publication dataset, `gencode.v47.transcripts`, measures motifs
across **genomic transcript spans, including introns**. Its `transcript_length`
is the span length, not the length of a spliced RNA molecule. Interpret results
and manuscript claims accordingly.

## Inputs

Run commands below from the `nonb-pipeline/` directory. The default dataset
requires these external files, which are absent from a fresh checkout:

| Path | Format and role |
|---|---|
| `resources/gencode.v47.annotation.gtf` | GENCODE v47 GRCh38 GTF; transcript and exon records provide the annotation intervals. |
| `resources/gencode.v47.pc_transcripts.fa` | Protein-coding transcript FASTA; versioned transcript IDs in headers provide coding labels. |
| `resources/gencode.v47.lncRNA_transcripts.fa` | lncRNA transcript FASTA; versioned transcript IDs in headers provide lncRNA labels. |
| `resources/GRCh38.p14.genome.fa.gz` | GRCh38 assembly FASTA archive used to generate motif BED files. |

The GTF and both transcript FASTAs must describe the same GENCODE release.
The assembly must match their GRCh38 coordinates. Obtain these reference files
from GENCODE and the GRCh38 assembly provider; the pipeline checks for them
before constructing a production DAG. The workflow builds the motif resources
under `resources/GRCh38_NonBDNA/` using pinned `non-B_gfa` and `g4Discovery`
revisions. The nine clean BED inputs are named
`GCA_000001405.15_GRCh38_no_alt_analysis_set.{motif}_clean.bed`, where
`{motif}` can be:

- `gfa.APR`
- `gfa.DR`
- `gfa.IR`
- `gfa.MR`
- `gfa.STR`
- `gfa.TRI`
- `gfa.Z`
- `g4Discovery_plus`
- `g4Discovery_minus`

For the `toy` dataset, Snakemake downloads a checksum-verified GENCODE v47
basic GTF and hg38 chromosome 22 FASTA. It extracts
`resources/toy/toy_chr22.gtf`, runs the pinned `non-B_gfa` and `g4Discovery`
tools on chromosome 22, and creates all nine
`resources/toy/toy_{motif}_chr22_clean.bed` files. A fresh checkout needs
network access, but no toy BED fixtures. The exact filenames are listed in
`workflow/rules/toy.smk`. Toy transcript labels come from the GTF's
`transcript_type` attribute. Production labels come from the two transcript
FASTAs listed above.

## Requirements

- Snakemake 9 with Conda or Mamba.
- Network access on first use to fetch the pinned Snakemake bedtools wrapper,
  create rule environments, download the toy GTF and chr22 FASTA, or clone
  upstream motif tools.
- For motif regeneration: the rule environments provide the C compiler, R and
  Python packages; the workflow also uses standard shell tools. No GPU is
  required.
- The default profile uses SLURM. Use `--executor local` for a local run.

## Run

From `nonb-pipeline/`, after providing the production inputs:

```bash
# Local production run
snakemake --executor local --cores 8 --use-conda all

# SLURM production run (uses profiles/default/config.yaml)
snakemake --profile profiles/default all
```

For the smaller chromosome 22 check, Snakemake builds the toy inputs:

```bash
snakemake --executor local --cores 1 --use-conda --config datasets='[toy]' all
```

To build only the genome-wide motif BED resources, run `all_upstream` with the
same local or SLURM options. To make the optional combined text report, target
`results/gencode.v47.transcripts/complete_analysis_report.txt`.

## Configuration

`Snakefile` loads `config/config.yaml`, `config/samples.yaml`, and
`config/toy.yaml`, in that order. Command-line `--config` and `--configfile`
values can override loaded values. Copy a config before editing it if you want
to keep the default publication setup intact.

### Functioning parameters

| Parameter | Default | Meaning |
|---|---|---|
| `datasets` | `[gencode.v47.transcripts]` | Dataset names to build under `results/{dataset}/`; each name needs a `samples` entry. |
| `samples.gencode.v47.transcripts.gtf` | `resources/gencode.v47.annotation.gtf` | Production GTF input. |
| `samples.gencode.v47.transcripts.feature_source` | `transcripts` | Uses full genomic transcript BED intervals, including introns. `exons` instead uses exon intervals and summed exon length, but requires a matching dataset/configuration and regenerated results. |
| `upstream.assembly_gz` | `resources/GRCh38.p14.genome.fa.gz` | External assembly archive for motif generation. |
| `upstream.g4_contigs` | `chr1`–`chr22`, `chrX`, `chrY`, `chrM` | Contigs processed separately by g4Discovery. |
| `gfa_motifs` | Nine motifs | Declared in config but not used; the active motif list is fixed in the Snakefile. |
| `analysis.*` | `run_basic: true`, `run_extended: true`, `keep_other_transcripts: false`, `alpha: 0.05`, `correction_method: fdr` | Declared but not read by the active workflow; changing them does not change analysis behavior. |

For a one-off toy run, use `--config datasets='[toy]'` as shown above. The
`upstream` section also pins tool repositories and Git revisions; see
`config/config.yaml` before changing resource generation.

### Environment

| Parameter | Default | Meaning |
|---|---|---|
| `profiles/default/config.yaml:executor` | `slurm` | Default executor; override with `--executor local` for local work. |
| `profiles/default/config.yaml:jobs` | `100` | Maximum concurrent SLURM jobs. |
| `profiles/default/config.yaml:use-conda` | `true` | Activate rule-specific Conda environments. |
| `profiles/default/config.yaml:latency-wait` | `40` seconds | Wait for outputs on shared storage. |
| `--cores` | `8` local / scheduler allocation on SLURM | Local CPU limit; choose a value supported by the machine. |

## Outputs

Outputs are under `results/{dataset}/`. The feature matrix for the configured
publication dataset is the principal result for downstream analysis:

| File | Meaning |
|---|---|
| `extended_analysis/features_nonb_features.csv` | One row per annotated transcript, with motif features, span length and coding class. Transcripts without motif hits retain zero-valued motif features; IDs absent from both production FASTAs can have class `other`. |
| `extended_analysis/contingency_motif_type_chi_square.csv` | Coding versus lncRNA motif-presence tests. |
| `extended_analysis/statistics_univariate_tests.csv` | Per-feature coding versus lncRNA tests and adjusted p values. |
| `extended_analysis/statistics_feature_importance.csv` | Random Forest feature importance. |
| `transcript_gfa.{motif}_summary.tsv`, `gfa.{motif}_distributions.png` | Basic per-motif transcript summaries and KDE plots. |
| `extended_analysis/analysis_complete.txt` | Completion marker, not an analysis table. |
| `complete_analysis_report.txt` | Optional combined text report; not part of `all`. |

## Reproducibility and limits

The motif tool Git revisions, toy annotation and FASTA URLs and MD5s, and rule
Conda environments are declared in config and `workflow/envs/`. The bedtools
wrapper is pinned at `v7.3.0`. Toy analysis passed twice with the earlier local
BED files on 2026-09-29; the new toy BED generation rules have passed a forced
resource DAG dry run but still need a clean-checkout integration run. For
resource boundaries and the rule inventory, see `nbd_pipeline_manifest.yaml`.

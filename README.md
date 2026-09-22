# Non-B DNA Analysis Pipeline

Snakemake workflow that intersects human transcript annotations (GTF) with
Non-B DNA structural-motif BED files and produces per-transcript feature
summaries, enrichment tests, and statistical analyses comparing protein-coding
and lncRNA transcripts.

## Requirements

- Conda / Mamba (environments: `workflow/envs/gfa.yaml`, `g4discovery.yaml`)
- Snakemake ≥ 7
- Internet access for the pinned upstream tool clones when regenerating motifs

## Non-B Motif Generation

`all_upstream` regenerates the GRCh38 motif BED files used by non-toy datasets.
It runs the default `non-B_gfa` predictor for APR, DR, IR, MR, STR, and Z-DNA;
G4 prediction is skipped there because it is produced separately by the default
`g4Discovery` workflow. TRI is the subset of MR records whose gfa `Subset`
field is flagged as triplex-prone. G4Discovery is run once per configured GRCh38
contig because it accepts a single FASTA record per invocation; the merged BED
retains the pqsfinder and G4Hunter score columns and is then split by strand.

The tool URLs, pinned revisions, FASTA paths, and contigs are declared under
`upstream` in `config/config.yaml`. Generate the resource set with:

```bash
snakemake --use-conda --cores 1 all_upstream
```

## Quick start

```bash
# copy and edit the default config
cp config/config.yaml config/config_mine.yaml

snakemake --snakefile Snakefile_integrated \
          --configfile config/config_mine.yaml \
          --use-conda -j 8

# Extended analysis only
snakemake --snakefile Snakefile_integrated \
          --configfile config/config_mine.yaml \
          --use-conda -j 8 extended_analysis_all
```

## Analysis layers

| Layer | Entry rule | Description |
|-------|-----------|-------------|
| **Basic** | `all` | Per-motif overlap summaries + KDE distribution plots |
| **Extended** | `extended_analysis_all` | 100+ feature extraction, chi-square contingency tests, univariate statistical tests, and Random Forest feature importance |

## Pipeline steps

| Step | Rule | Output |
|------|------|--------|
| 1 | `create_transcripts_bed` | `transcripts.bed` (BED6 from GTF) |
| 1 | `create_biotypes_from_fasta` | `biotypes.tsv` |
| 1 | `prepare_transcript_ids` | `annotation/pc_transcript_ids.txt`, `lncrna_transcript_ids.txt` |
| 2 | `basic_motif_analysis` | `transcript_gfa.{motif}_summary.tsv`, `gfa.{motif}_distributions.png` (×8 motifs) |
| 3.1 | `extended_feature_extraction` | `extended_analysis/features_nonb_features.csv`, `features_nonb_summary.txt` |
| 3.2 | `extended_contingency_analysis` | `extended_analysis/contingency_motif_type_chi_square.csv`, `contingency_contingency_report.txt` |
| 3.3 | `extended_statistical_analysis` | `extended_analysis/statistics_univariate_tests.csv`, `statistics_feature_importance.csv` |
| 4 | `create_summary_report` | `complete_analysis_report.txt` |

All outputs are namespaced under `results/{dataset}/`.

## Configs

| File | Purpose |
|------|---------|
| `config/config.yaml` | Default config (toy + gencode.v47 datasets) |

## Dev info

See `nbd_pipeline_manifest.yaml` for the full rule catalogue, known bugs, and
script inventory.

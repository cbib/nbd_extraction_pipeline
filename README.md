# Non-B DNA Analysis Pipeline

Snakemake workflow that intersects human transcript annotations (GTF) with
Non-B DNA structural-motif BED files and produces per-transcript feature
summaries, enrichment tests, and statistical analyses comparing protein-coding
and lncRNA transcripts.

## Requirements

- Conda / Mamba (environments: `workflow/envs/gfa.yaml`,
    `workflow/envs/g4discovery.yaml`, and `workflow/envs/nonb_analysis.yaml`)
- Snakemake ≥ 7
- Internet access when a pinned upstream tool, full reference resource, or toy
    annotation archive has not already been acquired

## Resource Policy

The workflow checks configured external prerequisites before expanding the DAG.
Non-toy runs require the configured Gencode annotation and transcript FASTAs,
plus the compressed GRCh38 assembly. Missing upstream tools and derived motif
BED resources are built by pinned workflow rules; missing external resources
produce an actionable preflight error rather than silently changing datasets.

Resource provenance and the declared external/buildable boundary are recorded
in `nbd_pipeline_manifest.yaml`.

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

The active motif set contains nine entries: `APR`, `DR`, `g4Discovery_plus`,
`g4Discovery_minus`, `IR`, `MR`, `STR`, `TRI`, and `Z`. The G4 outputs are
separate plus- and minus-strand BED files.

## Quick start

```bash
# copy and edit the default config; it enables gencode.v47.transcripts
cp config/config.yaml config/config_mine.yaml

snakemake --snakefile Snakefile \
          --configfile config/config_mine.yaml \
          --use-conda --cores 8

# Extended analysis for every configured dataset
snakemake --snakefile Snakefile \
          --configfile config/config_mine.yaml \
          --use-conda --cores 8 extended_analysis_all_datasets
```

## Toy Dataset

Run the toy workflow explicitly; the default configuration does not enable it:

```bash
snakemake --snakefile Snakefile --configfile config/config.yaml \
          --use-conda --cores 1 --config datasets='[toy]'
```

On its first run, the workflow downloads the pinned Gencode v47 basic
annotation, verifies its MD5, extracts `chr22`, and writes
`resources/toy/toy_provenance.yaml`. The provenance file records checksums and
record counts for the extracted annotation and the checked-in motif BED
fixtures. Source URL, checksum, assembly, chromosome, and output paths are
declared in `config/toy.yaml`.

## Analysis layers

| Layer | Entry rule | Description |
|-------|-----------|-------------|
| **Basic** | `all` | Per-motif overlap summaries + KDE distribution plots |
| **Extended** | `extended_analysis_all_datasets` | 100+ feature extraction, chi-square contingency tests, univariate statistical tests, and Random Forest feature importance |

## Pipeline steps

| Step | Rule | Output |
|------|------|--------|
| 1 | `create_transcripts_bed` | `transcripts.bed` (BED6 from GTF) |
| 1 | `create_biotypes_from_fasta` | `biotypes.tsv` |
| 1 | `prepare_transcript_ids` | `annotation/pc_transcript_ids.txt`, `lncrna_transcript_ids.txt` |
| 2 | `basic_motif_analysis` | `transcript_gfa.{motif}_summary.tsv`, `gfa.{motif}_distributions.png` (×9 motifs) |
| 3.1 | `extended_feature_extraction` | `extended_analysis/features_nonb_features.csv`, `features_nonb_summary.txt` |
| 3.2 | `extended_contingency_analysis` | `extended_analysis/contingency_motif_type_chi_square.csv`, `contingency_contingency_report.txt` |
| 3.3 | `extended_statistical_analysis` | `extended_analysis/statistics_univariate_tests.csv`, `statistics_feature_importance.csv` |
| 4 | `create_summary_report` | `complete_analysis_report.txt` |

All outputs are namespaced under `results/{dataset}/`.

## Configs

| File | Purpose |
|------|---------|
| `config/config.yaml` | Default configuration (`gencode.v47.transcripts`) |
| `config/samples.yaml` | Dataset annotation, motif BED, and feature-source settings |
| `config/toy.yaml` | Pinned toy annotation source and reproducible-output paths |

## Toy Dataset Labels

The reproducible toy dataset extracts `chr22` from the pinned Gencode v47 GTF
and derives its coding/lncRNA labels from the GTF `transcript_type` attribute.
This is a temporary workaround for a deterministic toy workflow; it is not
equivalent to the transcript classification required to reproduce the article's
findings. Production datasets must import coding/lncRNA labels derived from the
presence of transcripts in the GENCODE fasta files.

## Dev info

See `nbd_pipeline_manifest.yaml` for the full rule catalogue, known bugs, and
script inventory.

# Testing Framework Design

**Date:** 2026-06-19
**Branch:** feat/strand-specific-g4

## Goals

1. **Regression detection** — catch unintended output changes via content snapshots committed to git.
2. **Fast iteration** — unit tests run in seconds without invoking Snakemake.

## Approach: Two-layer pytest

### Layer 1 — Unit tests (`tests/unit/`)

Test the four core Python scripts directly by importing their classes/functions. No Snakemake involved.

**Fixtures (conftest.py):**
- *Synthetic*: small in-memory DataFrames (~10 rows) for fast logic assertions. Self-contained, no filesystem dependency.
- *Toy file paths*: `pytest.fixture` returning paths into `resources/toy/` for snapshot comparison runs.

**Covered scripts:**
| Test file | Target script |
|-----------|--------------|
| `test_feature_extractor.py` | `workflow/scripts/nonb_feature_extractor.py` |
| `test_parse_overlaps.py` | `workflow/scripts/parse_overlaps.py` |
| `test_contingency_analyzer.py` | `workflow/scripts/nonb_contingency_analyzer.py` |
| `test_statistical_analyzer.py` | `workflow/scripts/nonb_statistical_analyzer.py` |

Wrapper scripts (`extended_*_wrapper.py`) are Snakemake glue — covered by the integration test only.

**Test pattern per file:**
```
test_<script>.py
├── test_<logic_1>_synthetic   # fast, in-memory
├── test_<logic_2>_synthetic
└── test_output_matches_snapshot  # uses toy files, diffs against tests/snapshots/unit/
```

**Snapshot comparison:**
```python
pd.testing.assert_frame_equal(actual, expected, check_exact=False, atol=1e-6)
```
Float columns get `atol=1e-6` tolerance; string/integer columns are exact.

### Layer 2 — Integration test (`tests/integration/`)

One test function that runs the full pipeline on the toy dataset and checks outputs.

**Steps:**
1. `subprocess.run(["snakemake", "--config", "datasets=[toy]", "-j1", "--use-conda"], cwd=<pipeline_root>)`
2. Assert exit code 0.
3. For each expected output in `results/toy/`, assert it exists and matches `tests/snapshots/integration/`.
4. `yield`-based fixture cleans up `results/toy/` after the test.

Marked `@pytest.mark.integration` — excluded from default `pytest tests/unit/` runs, included in CI.

## Directory Layout

```
tests/
├── conftest.py
├── snapshots/
│   ├── unit/
│   │   ├── features_nonb_features.csv
│   │   ├── transcript_gfa.APR_summary.tsv
│   │   └── ...
│   └── integration/
│       ├── transcript_gfa.APR_summary.tsv
│       └── ...
├── unit/
│   ├── test_feature_extractor.py
│   ├── test_parse_overlaps.py
│   ├── test_contingency_analyzer.py
│   └── test_statistical_analyzer.py
└── integration/
    └── test_pipeline_toy.py
```

Snapshot files are TSV/CSV (same format as pipeline output) — human-readable and git-diffable.

## Snapshot Updates

`conftest.py` registers `--update-snapshots` via `pytest_addoption`. When the flag is set, snapshot tests write actual output to `tests/snapshots/` instead of comparing. When absent, they read and diff.

```bash
# Update unit snapshots after intentional logic change
pytest tests/unit/ --update-snapshots

# Update integration snapshots after pipeline output changes
pytest tests/integration/ --update-snapshots
```

## Running Tests

```bash
# Fast — unit tests only (~seconds)
pytest tests/unit/

# Full — unit + integration (~minutes, needs conda + toy data)
pytest tests/

# Integration only
pytest tests/integration/ -m integration

# Regenerate all snapshots
pytest tests/ --update-snapshots
```

## What Is Not Tested

- Shell/awk rules (`create_transcripts_bed`, `create_exons_bed`, etc.) — covered implicitly by the integration test.
- Snakemake DAG structure (wildcard expansion, rule ordering) — Snakemake validates this at dry-run time; the integration test exercises it end-to-end.
- Upstream rules (`workflow/rules/upstream.smk`) — require large external resources, out of scope.

## Constraints

- No new testing libraries beyond `pytest`. Snapshot logic is implemented directly in `conftest.py`.
- Toy dataset (`resources/toy/`) must be present for both unit snapshot tests and the integration test.
- Integration test requires a working `snakemake` + `conda` install in the environment.

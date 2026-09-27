"""Integration test: full pipeline on the toy dataset.

Requires:
  - snakemake + conda available in PATH
  - resources/toy/ data present
  - run with: pytest tests/integration/ -m integration
"""

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

PIPELINE_ROOT = Path(__file__).parent.parent.parent
TOY_RESULTS = PIPELINE_ROOT / "results" / "toy"
SNAPSHOTS_DIR = Path(__file__).parent.parent / "snapshots" / "integration"

# Files produced by the pipeline that we compare against snapshots
EXPECTED_OUTPUTS = [
    TOY_RESULTS / "extended_analysis" / "features_nonb_features.csv",
    TOY_RESULTS / "extended_analysis" / "contingency_motif_type_chi_square.csv",
    TOY_RESULTS / "extended_analysis" / "statistics_univariate_tests.csv",
]


@pytest.fixture(scope="module")
def run_pipeline(request):
    """Run snakemake on the toy dataset; clean results/toy/ after the test module."""
    update = request.config.getoption("--update-snapshots")

    result = subprocess.run(
        [
            "snakemake",
            "--config",
            "datasets=[toy]",
            "-j1",
            "--use-conda",
            "--rerun-incomplete",
        ],
        cwd=str(PIPELINE_ROOT),
        capture_output=True,
        text=True,
    )
    yield result, update

    # Cleanup: remove toy results after the test module finishes
    if TOY_RESULTS.exists():
        shutil.rmtree(TOY_RESULTS)


@pytest.mark.integration
def test_pipeline_exits_successfully(run_pipeline):
    result, _ = run_pipeline
    assert (
        result.returncode == 0
    ), f"Snakemake failed (exit {result.returncode}):\n{result.stderr[-3000:]}"


@pytest.mark.integration
def test_expected_output_files_exist(run_pipeline):
    _, _ = run_pipeline
    for path in EXPECTED_OUTPUTS:
        assert path.exists(), f"Missing pipeline output: {path}"


@pytest.mark.integration
def test_outputs_match_snapshots(run_pipeline):
    _, update = run_pipeline
    SNAPSHOTS_DIR.mkdir(parents=True, exist_ok=True)

    for output_path in EXPECTED_OUTPUTS:
        if not output_path.exists():
            pytest.skip(f"Output missing: {output_path}")

        snapshot = SNAPSHOTS_DIR / output_path.name
        sep = "\t" if output_path.suffix == ".tsv" else ","

        actual = pd.read_csv(output_path, sep=sep)
        if update:
            actual.to_csv(snapshot, sep=sep, index=False)
            continue

        assert snapshot.exists(), (
            f"Snapshot missing for {output_path.name}. "
            "Run with --update-snapshots to create it."
        )
        expected = pd.read_csv(snapshot, sep=sep)
        float_cols = actual.select_dtypes("float").columns.tolist()
        other_cols = [c for c in actual.columns if c not in float_cols]
        pd.testing.assert_frame_equal(
            actual[other_cols].reset_index(drop=True),
            expected[other_cols].reset_index(drop=True),
        )
        pd.testing.assert_frame_equal(
            actual[float_cols].reset_index(drop=True),
            expected[float_cols].reset_index(drop=True),
            check_exact=False,
            atol=1e-6,
        )

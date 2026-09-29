"""Unit tests for nonb_statistical_analyzer.NonBStatisticalAnalyzer."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))
from nonb_statistical_analyzer import NonBStatisticalAnalyzer

PIPELINE_ROOT = Path(__file__).parent.parent.parent
TOY_RESULTS = PIPELINE_ROOT / "results" / "toy"


def _make_analyzer(tmp_path, features_df):
    csv = tmp_path / "features.csv"
    features_df.to_csv(csv, index=False)
    return NonBStatisticalAnalyzer(str(csv), str(tmp_path / "out"))


# ---------------------------------------------------------------------------
# Synthetic tests
# ---------------------------------------------------------------------------


def test_perform_univariate_tests_returns_expected_columns(
    tmp_path, synthetic_features_df
):
    analyzer = _make_analyzer(tmp_path, synthetic_features_df)
    analyzer.load_features()
    results = analyzer.perform_univariate_tests()
    for col in (
        "feature",
        "mean_coding",
        "mean_lncrna",
        "t_statistic",
        "p_value_ttest",
        "u_statistic",
        "p_value_mannwhitney",
        "cohens_d",
    ):
        assert col in results.columns, f"Missing column: {col}"


def test_perform_univariate_tests_p_values_in_range(tmp_path, synthetic_features_df):
    analyzer = _make_analyzer(tmp_path, synthetic_features_df)
    analyzer.load_features()
    results = analyzer.perform_univariate_tests()
    assert (results["p_value_ttest"].dropna().between(0, 1)).all()
    assert (results["p_value_mannwhitney"].dropna().between(0, 1)).all()


def test_load_features_raises_without_coding_class(tmp_path):
    bad = pd.DataFrame({"transcript_id": ["T1"], "val": [1.0]})
    csv = tmp_path / "bad.csv"
    bad.to_csv(csv, index=False)
    analyzer = NonBStatisticalAnalyzer(str(csv), str(tmp_path / "out"))
    with pytest.raises(ValueError, match="coding_class"):
        analyzer.load_features()


# ---------------------------------------------------------------------------
# Snapshot test
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (TOY_RESULTS / "extended_analysis" / "features_nonb_features.csv").exists(),
    reason="toy data not present",
)
def test_output_matches_snapshot(
    tmp_path, toy_paths, snapshot_dir_unit, update_snapshots
):
    """Reproduce univariate test results on toy features; diff against snapshot."""
    analyzer = NonBStatisticalAnalyzer(
        str(toy_paths["features_csv"]),
        str(tmp_path / "stats"),
    )
    analyzer.load_features()
    actual = analyzer.perform_univariate_tests()

    snapshot = snapshot_dir_unit / "statistics_univariate_tests.csv"
    if update_snapshots:
        actual.to_csv(snapshot, index=False)
        return

    assert (
        snapshot.exists()
    ), "Snapshot missing. Run with --update-snapshots to create it."
    expected = pd.read_csv(snapshot)
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

"""Unit tests for nonb_contingency_analyzer.NonBContingencyAnalyzer."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))
from nonb_contingency_analyzer import NonBContingencyAnalyzer

PIPELINE_ROOT = Path(__file__).parent.parent.parent
TOY_RESULTS = PIPELINE_ROOT / "results" / "toy"


def _make_analyzer(tmp_path, features_df):
    csv = tmp_path / "features.csv"
    features_df.to_csv(csv, index=False)
    return NonBContingencyAnalyzer(str(csv), str(tmp_path / "out"))


# ---------------------------------------------------------------------------
# Synthetic tests
# ---------------------------------------------------------------------------


def test_perform_overall_test_returns_tuple(tmp_path, synthetic_features_df):
    analyzer = _make_analyzer(tmp_path, synthetic_features_df)
    analyzer.load_features()
    chi2, p_value, significant = analyzer.perform_overall_test()
    assert chi2 >= 0
    assert 0.0 <= p_value <= 1.0
    assert significant in (True, False)


def test_perform_motif_type_tests_synthetic(tmp_path, synthetic_features_df):
    analyzer = _make_analyzer(tmp_path, synthetic_features_df)
    analyzer.load_features()
    results = analyzer.perform_motif_type_tests()
    # Only motifs whose presence column exists in the fixture
    assert {"APR", "DR", "GQ_PLUS", "GQ_MINUS"} <= set(
        results["motif_type"].str.upper()
    )
    assert "p_value" in results.columns
    assert "cramers_v" in results.columns
    assert (results["p_value"] >= 0).all()
    assert (results["p_value"] <= 1).all()


def test_load_features_raises_without_coding_class(tmp_path):
    bad = pd.DataFrame({"transcript_id": ["T1"], "some_feature": [1.0]})
    csv = tmp_path / "bad.csv"
    bad.to_csv(csv, index=False)
    analyzer = NonBContingencyAnalyzer(str(csv), str(tmp_path / "out"))
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
    """Reproduce contingency results on toy features; diff against snapshot."""
    analyzer = NonBContingencyAnalyzer(
        str(toy_paths["features_csv"]),
        str(tmp_path / "contingency"),
    )
    analyzer.load_features()
    motif_results = analyzer.perform_motif_type_tests()

    snapshot = snapshot_dir_unit / "contingency_motif_type_chi_square.csv"
    if update_snapshots:
        motif_results.to_csv(snapshot, index=False)
        return

    assert (
        snapshot.exists()
    ), "Snapshot missing. Run with --update-snapshots to create it."
    expected = pd.read_csv(snapshot)
    float_cols = motif_results.select_dtypes("float").columns.tolist()
    other_cols = [c for c in motif_results.columns if c not in float_cols]
    pd.testing.assert_frame_equal(
        motif_results[other_cols].reset_index(drop=True),
        expected[other_cols].reset_index(drop=True),
    )
    pd.testing.assert_frame_equal(
        motif_results[float_cols].reset_index(drop=True),
        expected[float_cols].reset_index(drop=True),
        check_exact=False,
        atol=1e-6,
    )

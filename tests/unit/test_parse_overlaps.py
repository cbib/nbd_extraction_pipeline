"""Unit tests for parse_overlaps.py (legacy GTF-based path).

parse_overlaps.py uses `snakemake` as a global at module level; we inject a
mock before import so the module-level logging setup succeeds.
"""

import builtins
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

# Inject mock snakemake before the module is imported.
# setup_snakemake_logging checks `snakemake_obj.log`; falsy → basicConfig only.
_snakemake_mock = MagicMock()
_snakemake_mock.log = []
builtins.snakemake = _snakemake_mock

PIPELINE_ROOT = Path(__file__).parent.parent.parent
SCRIPTS_DIR = PIPELINE_ROOT / "workflow" / "scripts"
WORKFLOW_DIR = PIPELINE_ROOT / "workflow"
TOY_RESULTS = PIPELINE_ROOT / "results" / "toy"

for p in (str(SCRIPTS_DIR), str(WORKFLOW_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

import parse_overlaps  # noqa: E402 — must follow mock injection

# ---------------------------------------------------------------------------
# Synthetic tests
# ---------------------------------------------------------------------------


def test_filter_features_keeps_only_target_feature():
    df = pd.DataFrame(
        {
            "feature": ["gene", "transcript", "exon", "transcript"],
            "value": [1, 2, 3, 4],
        }
    )
    result = parse_overlaps.filter_features(df, feature="transcript")
    assert list(result["feature"]) == ["transcript", "transcript"]
    assert list(result["value"]) == [2, 4]


def test_extract_attributes_transcript_id():
    df = pd.DataFrame(
        {
            "feature": ["transcript", "transcript"],
            "attributes": [
                'gene_id "ENSG1.1"; transcript_id "ENST1.1"; gene_name "FOO";',
                'gene_id "ENSG2.1"; transcript_id "ENST2.1"; gene_name "BAR";',
            ],
        }
    )
    result = parse_overlaps.extract_attributes(df, attr="transcript_id")
    assert list(result["transcript_id"]) == ["ENST1.1", "ENST2.1"]


def test_extract_attributes_rejects_unknown():
    df = pd.DataFrame({"attributes": ['transcript_id "X";']})
    with pytest.raises(ValueError, match="not recognized"):
        parse_overlaps.extract_attributes(df, attr="bad_attr")


def test_classify_transcripts_assigns_correct_classes():
    df = pd.DataFrame(
        {
            "transcript_id": ["ENST1.1", "ENST2.1", "ENST3.1"],
        }
    )
    pc_ids = {"ENST1.1"}
    lnc_ids = {"ENST2.1"}
    result = parse_overlaps.classify_transcripts(df, pc_ids, lnc_ids, keep_other=True)
    classes = dict(zip(result["transcript_id"], result["coding_class"]))
    assert classes["ENST1.1"] == "coding"
    assert classes["ENST2.1"] == "lncRNA"
    assert classes["ENST3.1"] == "other"


def test_classify_transcripts_filters_other_by_default():
    df = pd.DataFrame({"transcript_id": ["ENST1.1", "ENST_OTHER.1"]})
    result = parse_overlaps.classify_transcripts(df, {"ENST1.1"}, set())
    assert "ENST_OTHER.1" not in result["transcript_id"].values


def test_calculate_lengths():
    df = pd.DataFrame(
        {
            "start": [100, 500],
            "end": [200, 600],
            "chrom_start": [110, 510],
            "chrom_end": [130, 540],
        }
    )
    result = parse_overlaps.calculate_lengths(df)
    assert list(result["transcript_length"]) == [101, 101]
    assert list(result["nbd_length"]) == [21, 31]


def test_aggregate_transcript_metrics_counts_hits():
    df = pd.DataFrame(
        {
            "transcript_id": ["T1", "T1", "T2"],
            "transcript_length": [1000, 1000, 2000],
            "coding_class": ["coding", "coding", "lncRNA"],
            "chrom": ["chr1", "chr1", "chr1"],
            "nbd_length": [50, 30, 100],
        }
    )
    result = parse_overlaps.aggregate_transcript_metrics(df, "g4")
    assert result.loc["T1", "g4_hit_count"] == 2
    assert result.loc["T2", "g4_hit_count"] == 1
    # length-normalized: hits per kb
    assert abs(result.loc["T1", "g4_hit_count_per_kb"] - 2.0) < 1e-9


# ---------------------------------------------------------------------------
# Snapshot test
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (TOY_RESULTS / "isect_g4Discovery.bed").exists(),
    reason="toy data not present",
)
def test_output_matches_snapshot(toy_paths, snapshot_dir_unit, update_snapshots):
    """Reproduce transcript_g4Discovery_summary.tsv and diff against snapshot."""
    pc_ids, lnc_ids = parse_overlaps.load_transcript_ids(
        str(toy_paths["pc_ids"]), str(toy_paths["lnc_ids"])
    )
    df = parse_overlaps.read_overlap_file(str(toy_paths["isect_g4"]))
    df = parse_overlaps.filter_features(df)
    df = parse_overlaps.extract_attributes(df)
    df = parse_overlaps.classify_transcripts(df, pc_ids, lnc_ids, version=False)
    df = parse_overlaps.calculate_lengths(df)
    actual = parse_overlaps.aggregate_transcript_metrics(df, "g4Discovery")

    snapshot = snapshot_dir_unit / "transcript_g4Discovery_summary.tsv"
    if update_snapshots:
        actual.to_csv(snapshot, sep="\t", index=True)
        return

    assert (
        snapshot.exists()
    ), "Snapshot missing. Run with --update-snapshots to create it."
    expected = pd.read_csv(snapshot, sep="\t", index_col=0)
    float_cols = actual.select_dtypes("float").columns.tolist()
    other_cols = [c for c in actual.columns if c not in float_cols]
    pd.testing.assert_frame_equal(actual[other_cols], expected[other_cols])
    pd.testing.assert_frame_equal(
        actual[float_cols], expected[float_cols], check_exact=False, atol=1e-6
    )

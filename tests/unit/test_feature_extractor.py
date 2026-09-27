"""Unit tests for nonb_feature_extractor.NonBFeatureExtractor."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "workflow" / "scripts"))
from nonb_feature_extractor import NonBFeatureExtractor

PIPELINE_ROOT = Path(__file__).parent.parent.parent
TOY_RESULTS = PIPELINE_ROOT / "results" / "toy"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_extractor(tmp_path, transcripts_df, isect_dfs: dict):
    """Write temp BED files and return an initialised (not yet loaded) extractor."""
    bed = tmp_path / "transcripts.bed"
    transcripts_df.to_csv(bed, sep="\t", index=False, header=False)

    intersections = {}
    for motif, df in isect_dfs.items():
        f = tmp_path / f"isect_{motif}.bed"
        df.to_csv(f, sep="\t", index=False, header=False)
        intersections[motif] = str(f)

    return NonBFeatureExtractor(
        nonb_intersections=intersections,
        transcripts_bed=str(bed),
        output_prefix=str(tmp_path / "out"),
    )


# ---------------------------------------------------------------------------
# Synthetic tests
# ---------------------------------------------------------------------------


def test_load_data_aggregates_exon_lengths(tmp_path, synthetic_transcripts_bed_df):
    """load_data() should sum exon segment lengths per transcript."""
    extractor = _make_extractor(tmp_path, synthetic_transcripts_bed_df, {})
    extractor.load_data()

    result = extractor.transcripts.set_index("transcript_id")
    # ENST1.1 has two exons: [100,200) + [500,600) = 100 + 100 = 200
    assert result.loc["ENST1.1", "length"] == 200
    # ENST2.1 has one exon: [200,400) = 200
    assert result.loc["ENST2.1", "length"] == 200
    assert len(result) == 2


def test_load_data_bed6_schema_detection(
    tmp_path, synthetic_transcripts_bed_df, synthetic_isect_bed6_df
):
    """load_data() detects the BED6-based intersection schema (13 cols)."""
    extractor = _make_extractor(
        tmp_path,
        synthetic_transcripts_bed_df,
        {"APR": synthetic_isect_bed6_df},
    )
    extractor.load_data()
    assert "APR" in extractor.nonb_elements
    assert len(extractor.nonb_elements["APR"]) > 0


def test_resolve_overlaps_merges_overlapping_motifs(
    tmp_path, synthetic_transcripts_bed_df
):
    """resolve_overlaps() must collapse overlapping intervals into one segment."""
    # Two overlapping motifs on ENST1.1: [120,160) and [140,180)
    isect = pd.DataFrame(
        {
            0: ["chr1", "chr1"],
            1: [100, 100],  # exon start
            2: [200, 200],  # exon end
            3: ["ENST1.1", "ENST1.1"],
            4: [0, 0],
            5: ["+", "+"],
            6: ["chr1", "chr1"],
            7: [120, 140],  # motif start
            8: [160, 180],  # motif end
            9: [".", "."],
            10: [0, 0],
            11: ["+", "+"],
            12: ["0,0,0", "0,0,0"],
        }
    )
    extractor = _make_extractor(tmp_path, synthetic_transcripts_bed_df, {"APR": isect})
    extractor.load_data()
    extractor.resolve_overlaps()

    resolved = extractor.resolved_elements["APR"]
    assert len(resolved) == 1
    assert resolved["motif_start"].iloc[0] == 120
    assert resolved["motif_end"].iloc[0] == 180


def test_resolve_overlaps_keeps_non_overlapping_separate(
    tmp_path, synthetic_transcripts_bed_df
):
    """Two non-overlapping motifs on the same transcript stay as two segments."""
    isect = pd.DataFrame(
        {
            0: ["chr1", "chr1"],
            1: [100, 500],
            2: [200, 600],
            3: ["ENST1.1", "ENST1.1"],
            4: [0, 0],
            5: ["+", "+"],
            6: ["chr1", "chr1"],
            7: [110, 520],
            8: [130, 560],
            9: [".", "."],
            10: [0, 0],
            11: ["+", "+"],
            12: ["0,0,0", "0,0,0"],
        }
    )
    extractor = _make_extractor(tmp_path, synthetic_transcripts_bed_df, {"APR": isect})
    extractor.load_data()
    extractor.resolve_overlaps()

    assert len(extractor.resolved_elements["APR"]) == 2


# ---------------------------------------------------------------------------
# Snapshot test
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (TOY_RESULTS / "transcripts.bed").exists(),
    reason="toy data not present",
)
def test_output_matches_snapshot(
    tmp_path, toy_paths, snapshot_dir_unit, update_snapshots
):
    """Run feature extractor on toy data; compare CSV to stored snapshot."""
    ext_dir = toy_paths["extended_dir"]
    intersections = {}
    motif_map = {
        "APR": "isect_gfa.APR.bed",
        "DR": "isect_gfa.DR.bed",
        "IR": "isect_gfa.IR.bed",
        "MR": "isect_gfa.MR.bed",
        "STR": "isect_gfa.STR.bed",
        "TRI": "isect_gfa.TRI.bed",
        "Z": "isect_gfa.Z.bed",
        "GQ_PLUS": "isect_g4Discovery.bed",
    }
    for motif, fname in motif_map.items():
        p = ext_dir / fname
        if p.exists():
            intersections[motif] = str(p)

    extractor = NonBFeatureExtractor(
        nonb_intersections=intersections,
        transcripts_bed=str(toy_paths["transcripts_bed"]),
        biotypes_file=str(toy_paths["biotypes_tsv"]),
        output_prefix=str(tmp_path / "out"),
        motif_name_mapping={"GQ_PLUS": "g4discovery"},
    )
    extractor.load_data()
    extractor.resolve_overlaps()
    actual = extractor.extract_features()

    snapshot = snapshot_dir_unit / "features_nonb_features.csv"
    if update_snapshots:
        actual.to_csv(snapshot, index=False)
        return

    assert (
        snapshot.exists()
    ), f"Snapshot missing. Run with --update-snapshots to create it."
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

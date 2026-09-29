import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PIPELINE_ROOT = Path(__file__).parent.parent
SCRIPTS_DIR = PIPELINE_ROOT / "workflow" / "scripts"
WORKFLOW_DIR = PIPELINE_ROOT / "workflow"
TOY_RESULTS_DIR = PIPELINE_ROOT / "results" / "toy"
SNAPSHOTS_DIR = Path(__file__).parent / "snapshots"

# Add workflow paths so scripts can import each other and utils
for p in (str(SCRIPTS_DIR), str(WORKFLOW_DIR), str(PIPELINE_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)


def pytest_addoption(parser):
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        default=False,
        help="Overwrite snapshot files with actual output instead of comparing",
    )


@pytest.fixture
def update_snapshots(request):
    return request.config.getoption("--update-snapshots")


@pytest.fixture
def snapshot_dir_unit():
    d = SNAPSHOTS_DIR / "unit"
    d.mkdir(parents=True, exist_ok=True)
    return d


@pytest.fixture
def toy_paths():
    return {
        "transcripts_bed": TOY_RESULTS_DIR / "transcripts.bed",
        "extended_dir": TOY_RESULTS_DIR / "extended_analysis",
        "annotation_dir": TOY_RESULTS_DIR / "annotation",
        "features_csv": TOY_RESULTS_DIR
        / "extended_analysis"
        / "features_nonb_features.csv",
        "biotypes_tsv": TOY_RESULTS_DIR / "biotypes.tsv",
        "isect_g4": TOY_RESULTS_DIR / "isect_g4Discovery.bed",
        "pc_ids": TOY_RESULTS_DIR / "annotation" / "pc_transcript_ids.txt",
        "lnc_ids": TOY_RESULTS_DIR / "annotation" / "lncrna_transcript_ids.txt",
    }


# --- Synthetic fixtures (in-memory, no filesystem) ---


@pytest.fixture
def synthetic_transcripts_bed_df():
    """Two-transcript BED: ENST1.1 spans two exons, ENST2.1 is single."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [100, 500, 200],
            "end": [200, 600, 400],
            "transcript_id": ["ENST1.1", "ENST1.1", "ENST2.1"],
            "score": [0, 0, 0],
            "strand": ["+", "+", "-"],
        }
    )


@pytest.fixture
def synthetic_isect_bed6_df():
    """BED6+B format intersection: two overlapping motifs on ENST1.1."""
    return pd.DataFrame(
        {
            0: ["chr1", "chr1"],
            1: [100, 500],
            2: [200, 600],
            3: ["ENST1.1", "ENST1.1"],
            4: [0, 0],
            5: ["+", "+"],
            6: ["chr1", "chr1"],
            7: [120, 510],
            8: [160, 580],
            9: [".", "."],
            10: [0, 0],
            11: ["+", "+"],
            12: ["0,0,0", "0,0,0"],
        }
    )


@pytest.fixture
def synthetic_features_df():
    """Minimal feature DataFrame for contingency/statistical tests."""
    rng = np.random.default_rng(42)
    n = 40
    return pd.DataFrame(
        {
            "transcript_id": [f"ENST{i:05d}.1" for i in range(n)],
            "transcript_id_base": [f"ENST{i:05d}" for i in range(n)],
            "coding_class": ["coding"] * (n // 2) + ["lncRNA"] * (n // 2),
            "transcript_length": rng.integers(1000, 50000, n).tolist(),
            "apr_hit_count": rng.poisson(3, n).tolist(),
            "apr_total_length": rng.integers(50, 500, n).tolist(),
            "apr_present": (rng.random(n) > 0.4).astype(int).tolist(),
            "dr_hit_count": rng.poisson(2, n).tolist(),
            "dr_present": (rng.random(n) > 0.5).astype(int).tolist(),
            "gq_plus_hit_count": rng.poisson(5, n).tolist(),
            "gq_plus_present": (rng.random(n) > 0.3).astype(int).tolist(),
            "gq_minus_hit_count": rng.poisson(4, n).tolist(),
            "gq_minus_present": (rng.random(n) > 0.35).astype(int).tolist(),
            "any_nonb_present": (rng.random(n) > 0.25).astype(int).tolist(),
            "apr_hit_count_per_kb": rng.exponential(1.5, n).tolist(),
            "dr_hit_count_per_kb": rng.exponential(1.0, n).tolist(),
        }
    )

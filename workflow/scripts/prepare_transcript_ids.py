#!/usr/bin/env python3
"""
Extract transcript IDs for protein-coding and lncRNA transcripts.

Usage in Snakemake:
    script: "scripts/prepare_transcript_ids.py"

Input (via snakemake object):
    input.biotypes: Path to biotypes TSV file

Output (via snakemake object):
    output.pc_ids: Path to protein-coding IDs output
    output.lnc_ids: Path to lncRNA IDs output
"""

import sys
from pathlib import Path

import pandas as pd

# Add workflow directory to path
workflow_dir = Path(__file__).parent.parent.parent
sys.path.insert(0, str(workflow_dir))

try:
    from workflow.utils.logging_utils import setup_snakemake_logging

    logger = setup_snakemake_logging(snakemake, script_name=__file__)
except ImportError:
    import logging

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

# Read biotypes
logger.info(f"Reading biotypes from: {snakemake.input.biotypes}")
df = pd.read_csv(snakemake.input.biotypes, sep="\t")
logger.info(f"Loaded {len(df)} transcript biotypes")

# Extract protein-coding IDs
pc_df = df[df["transcript_type"] == "protein_coding"]
logger.info(f"Found {len(pc_df)} protein-coding transcripts")
pc_df["transcript_id_base"].to_csv(snakemake.output.pc_ids, index=False, header=False)

# Extract lncRNA IDs
lnc_df = df[df["transcript_type"].str.contains("lncRNA", na=False)]
logger.info(f"Found {len(lnc_df)} lncRNA transcripts")
lnc_df["transcript_id_base"].to_csv(snakemake.output.lnc_ids, index=False, header=False)

logger.info("Transcript ID extraction complete")

#!/usr/bin/env python3
"""
Snakemake wrapper for Non-B DNA Feature Extractor
Interfaces the feature extractor with Snakemake's input/output/params
"""

import sys
from pathlib import Path

# Add workflow directory to path
script_dir = Path(__file__).parent
workflow_dir = script_dir.parent.parent
sys.path.insert(0, str(workflow_dir))

from workflow.utils.logging_utils import log_job_completion, setup_snakemake_logging

# Setup logging
logger = setup_snakemake_logging(snakemake, script_name=__file__)

# Import the feature extractor
from nonb_feature_extractor import NonBFeatureExtractor

# Get parameters from Snakemake
input_dir = snakemake.params.input_dir
sample = snakemake.wildcards.dataset
transcripts_bed = str(snakemake.input.transcripts_bed)
biotypes_file = (
    str(snakemake.input.biotypes) if hasattr(snakemake.input, "biotypes") else None
)
output_prefix = snakemake.params.output_prefix
exon_mode = bool(snakemake.params.get("exon_mode", True))

logger.info("=" * 80)
logger.info("Non-B DNA Extended Feature Extraction")
logger.info("=" * 80)
logger.info(f"Sample: {sample}")
logger.info(f"Input directory: {input_dir}")
logger.info(f"Transcripts BED: {transcripts_bed}")
logger.info(f"Biotypes file: {biotypes_file}")
logger.info(f"Output prefix: {output_prefix}")
logger.info(f"Exon mode: {exon_mode}")

# Build intersections dictionary from explicit Snakemake inputs
# This ensures all files are tracked by Snakemake's dependency system
intersections = {
    "APR": str(snakemake.input.apr),
    "DR": str(snakemake.input.dr),
    "GQ_PLUS": str(snakemake.input.gq_plus),
    "GQ_MINUS": str(snakemake.input.gq_minus),
    "IR": str(snakemake.input.ir),
    "MR": str(snakemake.input.mr),
    "STR": str(snakemake.input._str),
    "TRI": str(snakemake.input.tri),
    "Z": str(snakemake.input.z),
}

# Map internal motif names to output file names
motif_name_mapping = {
    "GQ_PLUS": "g4Discovery_plus",
    "GQ_MINUS": "g4Discovery_minus",
}

logger.info(
    f"Processing {len(intersections)} intersection files: {list(intersections.keys())}"
)

# Initialize extractor
extractor = NonBFeatureExtractor(
    nonb_intersections=intersections,
    transcripts_bed=transcripts_bed,
    biotypes_file=biotypes_file,
    output_prefix=output_prefix,
    motif_name_mapping=motif_name_mapping,
    exon_mode=exon_mode,
)

# Run pipeline
try:
    extractor.load_data()
    extractor.resolve_overlaps()
    extractor.extract_features()
    extractor.save_features()
    extractor.save_resolved_elements()

    logger.info("=" * 80)
    logger.info("Feature extraction completed successfully!")
    logger.info("=" * 80)

    log_job_completion(logger)

except Exception as e:
    logger.error(f"Feature extraction failed: {e}", exc_info=True)
    sys.exit(1)

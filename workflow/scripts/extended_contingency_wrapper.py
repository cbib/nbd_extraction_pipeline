#!/usr/bin/env python3
"""
Snakemake wrapper for Non-B DNA Contingency Analyzer
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

# Import the analyzer
from nonb_contingency_analyzer import NonBContingencyAnalyzer

# Get parameters from Snakemake
features_file = str(snakemake.input.features)
output_prefix = snakemake.params.output_prefix

logger.info("=" * 80)
logger.info("Non-B DNA Contingency Analysis")
logger.info("=" * 80)
logger.info(f"Features file: {features_file}")
logger.info(f"Output prefix: {output_prefix}")

# Initialize analyzer
analyzer = NonBContingencyAnalyzer(
    features_file=features_file, output_prefix=output_prefix
)

# Run analysis
try:
    analyzer.load_features()
    overall_stats = analyzer.perform_overall_test()
    motif_results = analyzer.perform_motif_type_tests()
    analyzer.save_results(overall_stats, motif_results)

    logger.info("=" * 80)
    logger.info("Contingency analysis completed successfully!")
    logger.info("=" * 80)

    log_job_completion(logger)

except Exception as e:
    logger.error(f"Contingency analysis failed: {e}", exc_info=True)
    sys.exit(1)

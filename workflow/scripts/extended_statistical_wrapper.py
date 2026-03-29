#!/usr/bin/env python3
"""
Snakemake wrapper for Non-B DNA Statistical Analyzer
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
from nonb_statistical_analyzer import NonBStatisticalAnalyzer

# Get parameters from Snakemake
features_file = str(snakemake.input.features)
output_prefix = snakemake.params.output_prefix

logger.info("=" * 80)
logger.info("Non-B DNA Statistical Analysis")
logger.info("=" * 80)
logger.info(f"Features file: {features_file}")
logger.info(f"Output prefix: {output_prefix}")

# Initialize analyzer
analyzer = NonBStatisticalAnalyzer(
    features_file=features_file, output_prefix=output_prefix
)

# Run analysis
try:
    analyzer.load_features()
    univariate_results = analyzer.perform_univariate_tests()
    importance_results, performance_results = analyzer.perform_feature_importance()
    analyzer.save_results(univariate_results, importance_results, performance_results)

    logger.info("=" * 80)
    logger.info("Statistical analysis completed successfully!")
    logger.info("=" * 80)

    log_job_completion(logger)

except Exception as e:
    logger.error(f"Statistical analysis failed: {e}", exc_info=True)
    sys.exit(1)

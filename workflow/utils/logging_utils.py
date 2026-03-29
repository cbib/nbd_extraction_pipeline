import logging
import os
import sys
from datetime import datetime
from pathlib import Path


def setup_snakemake_logging(snakemake_obj, script_name=None):
    """
    Configure logging for Snakemake scripts with enhanced metadata.

    Parameters:
    -----------
    snakemake_obj : snakemake object
        The snakemake object available in script directive scripts
    script_name : str, optional
        Name of the calling script. If None, will attempt to determine automatically.

    Returns:
    --------
    logger : logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger()

    # Determine script name if not provided
    if script_name is None:
        # Try to get from snakemake object's scriptdir attribute
        if hasattr(snakemake_obj, "scriptdir"):
            # scriptdir gives the directory, we need to combine with the actual script
            script_name = "Unknown"
        else:
            script_name = "Unknown"

    if snakemake_obj.log:
        logging.basicConfig(
            filename=snakemake_obj.log[0],
            level=logging.DEBUG,
            format="%(asctime)s [%(levelname)s] %(message)s",
        )
        sys.stderr = open(snakemake_obj.log[0], "a")

        # Log header with metadata
        logger.info("=" * 80)
        logger.info("SNAKEMAKE JOB METADATA")
        logger.info("=" * 80)
        logger.info(f"Rule: {snakemake_obj.rule}")
        logger.info(f"Script: {script_name}")
        logger.info(f"Log file: {snakemake_obj.log[0]}")

        # SLURM job information
        slurm_job_id = os.environ.get("SLURM_JOB_ID", "N/A")
        slurm_job_name = os.environ.get("SLURM_JOB_NAME", "N/A")
        slurm_array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID", "N/A")
        logger.info(f"SLURM Job ID: {slurm_job_id}")
        logger.info(f"SLURM Job Name: {slurm_job_name}")
        if slurm_array_task_id != "N/A":
            logger.info(f"SLURM Array Task ID: {slurm_array_task_id}")

        # Wildcards
        if snakemake_obj.wildcards:
            logger.info(f"Wildcards: {dict(snakemake_obj.wildcards)}")

        # Resources
        logger.info(f"Threads: {snakemake_obj.threads}")
        if hasattr(snakemake_obj, "resources"):
            logger.info(f"Resources: {dict(snakemake_obj.resources)}")

        # Input/Output files
        logger.info(f"Input files: {list(snakemake_obj.input)}")
        logger.info(f"Output files: {list(snakemake_obj.output)}")

        # Timestamp
        logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 80)
        logger.info("")

        # Exception handler
        def handle_unhandled_exception(exc_type, exc_value, exc_traceback):
            if issubclass(exc_type, KeyboardInterrupt):
                sys.__excepthook__(exc_type, exc_value, exc_traceback)
                return
            logger.critical(
                "Unhandled exception", exc_info=(exc_type, exc_value, exc_traceback)
            )

        sys.excepthook = handle_unhandled_exception
    else:
        logging.basicConfig(level=logging.DEBUG)

    return logger


def log_job_completion(logger):
    """Log job completion metadata."""
    logger.info("")
    logger.info("=" * 80)
    logger.info(f"Job completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 80)

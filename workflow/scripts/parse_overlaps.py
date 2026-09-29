"""DEPRECATED: legacy basic-analysis parser for GTF-based intersection files.

The extended pipeline now uses BED-based transcript intersections together with
nonb_feature_extractor.py. Keep this script only for the older basic
analysis path until that workflow is retired or migrated.
"""

import sys
from pathlib import Path
from typing import Set, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Add workflow directory to sys.path for imports
script_dir = Path(__file__).parent
workflow_dir = script_dir.parent
sys.path.insert(0, str(workflow_dir))
print(f"Added {workflow_dir} to sys.path for imports")

from utils.logging_utils import log_job_completion, setup_snakemake_logging

logger = setup_snakemake_logging(snakemake, script_name=__file__)
logger.warning(
    "parse_overlaps.py is deprecated and only supports the legacy GTF-based "
    "intersection schema used by the basic analysis path."
)


def load_transcript_ids(pc_id_file: str, lnc_id_file: str) -> Tuple[Set[str], Set[str]]:
    """
    Load protein-coding and lncRNA transcript IDs from files.

    Args:
        pc_id_file: Path to protein-coding transcript IDs file
        lnc_id_file: Path to lncRNA transcript IDs file

    Returns:
        Tuple of (protein-coding IDs set, lncRNA IDs set)
    """
    logger.info(f"Loading protein-coding transcript IDs from: {pc_id_file}")
    pc_ids = set(pd.read_csv(pc_id_file, header=None)[0])
    logger.info(f"Loaded {len(pc_ids)} protein-coding transcript IDs")

    logger.info(f"Loading lncRNA transcript IDs from: {lnc_id_file}")
    lnc_ids = set(pd.read_csv(lnc_id_file, header=None)[0])
    logger.info(f"Loaded {len(lnc_ids)} lncRNA transcript IDs")

    return pc_ids, lnc_ids


def read_overlap_file(overlap_file: str) -> pd.DataFrame:
    """
    Read and parse the overlap BED file with transcript and NBD features.

    Args:
        overlap_file: Path to the overlap file

    Returns:
        DataFrame with parsed columns
    """
    logger.info(f"Reading overlap file: {overlap_file}")

    # Define expected columns
    transcript_cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]
    nbd_cols = [
        "chrom",
        "chrom_start",
        "chrom_end",
        "nbd_id",
        "nbd_score",
        "nbd_strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
    ]
    combined_cols = transcript_cols + nbd_cols

    # Read file
    df = pd.read_csv(overlap_file, sep="\t", header=None, dtype=str)
    logger.info(f"Read {len(df)} rows from overlap file")

    # Assign columns dynamically
    if len(df.columns) < len(combined_cols):
        logger.warning(
            f"Overlap file has {len(df.columns)} columns, expected {len(combined_cols)}. Truncating column names."
        )
        combined_cols = combined_cols[: len(df.columns)]

    df.columns = combined_cols
    logger.info(f"Assigned column names: {', '.join(df.columns)}")

    numeric_cols = ["start", "end", "chrom_start", "chrom_end"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def filter_features(df: pd.DataFrame, feature="transcript") -> pd.DataFrame:
    """
    Filter dataframe to only include transcript features.

    Args:
        df: Input DataFrame
        feature: Feature type to filter on (default: 'transcript')

    Returns:
        Filtered DataFrame containing only a specific feature type
    """
    logger.info(f"Filtering annotation DataFrame for {feature} features only")
    initial_count = len(df)
    df_filtered = df[df["feature"] == feature].copy()
    logger.info(
        f"Filtered from {initial_count} to {len(df_filtered)} {feature} records"
    )

    return df_filtered


def extract_attributes(
    df: pd.DataFrame, attr: Union[str, list] = "transcript_id"
) -> pd.DataFrame:
    """
    Extract transcript IDs from GTF attributes field.

    Args:
        df: Input DataFrame with 'attributes' column
        attr: Attribute to extract (default: 'transcript_id')

    Returns:
        DataFrame with a new column for the extracted attribute
    """
    available_attrs = ["transcript_id", "gene_id", "gene_name"]
    if isinstance(attr, str):
        attr = [attr]
    for a in attr:
        if a not in available_attrs:
            raise ValueError(
                f"Attribute '{a}' not recognized. Available attributes: {available_attrs}"
            )

        logger.info(f"Extracting {a} from attributes field")

        # Use regex extraction for Gencode format
        df[a] = df["attributes"].str.extract(rf'{a} "([^"]+)"')

        n_extracted = df[a].notna().sum()
        n_total = len(df)
        logger.info(f"Successfully extracted {n_extracted}/{n_total} {a} attributes")

        if n_extracted < n_total:
            logger.warning(f"{n_total - n_extracted} records missing {a} attribute")

    return df


def classify_transcripts(
    df: pd.DataFrame,
    pc_ids: Set[str],
    lnc_ids: Set[str],
    version=True,
    keep_other=False,
) -> pd.DataFrame:
    """
    Classify transcripts as coding, lncRNA, or other based on the provided ID sets

    Args:
        df: DataFrame with transcript_id column
        pc_ids: Set of protein-coding transcript IDs
        lnc_ids: Set of lncRNA transcript IDs
        version: Whether transcript IDs include version numbers (default: True)

    Returns:
        DataFrame with added 'coding_class' column
    """
    logger.info("Classifying transcripts")

    # Initialize all as 'other'
    df["coding_class"] = "other"

    # Extract base transcript ID
    if version:
        base_ids = df["transcript_id"]
    else:
        base_ids = df["transcript_id"].str.split(".").str[0]

    # Classify coding transcripts
    coding_mask = base_ids.isin(pc_ids)
    df.loc[coding_mask, "coding_class"] = "coding"
    n_coding = coding_mask.sum()

    # Classify lncRNA transcripts
    lnc_mask = base_ids.isin(lnc_ids)
    df.loc[lnc_mask, "coding_class"] = "lncRNA"
    n_lnc = lnc_mask.sum()

    # Count other
    n_other = (df["coding_class"] == "other").sum()

    logger.info(
        f"Classification results: {n_coding} coding, {n_lnc} lncRNA, {n_other} other"
    )

    if not keep_other:
        df = df[df["coding_class"] != "other"].copy()
        logger.info(f"Filtered out 'other' transcripts, remaining records: {len(df)}")

    return df


def calculate_lengths(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate transcript and NBD feature lengths.

    Args:
        df: DataFrame with coordinate columns

    Returns:
        DataFrame with added length columns
    """
    logger.info("Calculating transcript and NBD lengths")

    df["transcript_length"] = df["end"] - df["start"] + 1
    df["nbd_length"] = df["chrom_end"] - df["chrom_start"] + 1

    logger.info(
        f"Transcript length range: {df['transcript_length'].min():.0f} - {df['transcript_length'].max():.0f} bp"
    )
    logger.info(
        f"NBD length range: {df['nbd_length'].min():.0f} - {df['nbd_length'].max():.0f} bp"
    )

    return df


def aggregate_transcript_metrics(df: pd.DataFrame, nbd_type: str) -> pd.DataFrame:
    """
    Aggregate NBD metrics per transcript, normalized by transcript length.

    Args:
        df: DataFrame with per-overlap records
        nbd_type: Type of NBD feature being analyzed (e.g., 'g4', 'z_dna', 'triplex')

    Returns:
        DataFrame with per-transcript aggregated metrics (normalized by transcript length)
    """
    logger.info(f"Aggregating {nbd_type.upper()} metrics per transcript")

    # Create dynamic column names based on NBD type
    hit_count_col = f"{nbd_type}_hit_count"
    total_length_col = f"total_{nbd_type}_length"
    max_length_col = f"max_{nbd_type}_length"
    coverage_col = f"{nbd_type}_coverage"

    hits = df.groupby("transcript_id").agg(
        transcript_length=("transcript_length", "first"),
        coding_class=("coding_class", "first"),
        **{
            total_length_col: ("nbd_length", "sum"),
            max_length_col: ("nbd_length", "max"),
            hit_count_col: ("chrom", "size"),
        },
    )
    # NOTE: The **{} syntax is used to dynamically create aggregation columns.

    ############
    # Calculate length-normalized metrics (per kilobase)
    ############
    transcript_length_kb = hits["transcript_length"] / 1000

    # Hits per kb
    hits[f"{hit_count_col}_per_kb"] = hits[hit_count_col] / transcript_length_kb
    # Calculate coverage in percentage
    hits[coverage_col] = hits[total_length_col] / hits["transcript_length"] * 100
    # Length of largest NBD feature as percentage of transcript length
    hits[f"{max_length_col}_as_%"] = (
        hits[max_length_col] / hits["transcript_length"] * 100
    )

    # Drop raw columns (keep them if you want to save both)
    # hits = hits.drop(columns=['raw_hit_count', 'raw_total_length', 'raw_max_length'])

    logger.info(
        f"Aggregated metrics for {len(hits)} unique transcripts (normalized per kb)"
    )

    logger.info(f"Aggregated metrics for {len(hits)} unique transcripts")

    return hits


def print_summary_statistics(hits: pd.DataFrame, nbd_type: str) -> None:
    """
    Print summary statistics (length-normalized).

    Args:
        hits: DataFrame with aggregated transcript metrics
        nbd_type: Type of NBD feature being analyzed
    """
    nbd_label = nbd_type.upper()

    logger.info("=" * 60)
    logger.info(f"SUMMARY STATISTICS - {nbd_label} (Length-Normalized)")
    logger.info("=" * 60)

    # Dynamic column names (now normalized per kb)
    hit_count_col = f"{nbd_type}_hit_count_per_kb"
    # total_length_col = f'total_{nbd_type}_length_per_kb'
    max_length_col = f"max_{nbd_type}_length_as_%"
    coverage_col = f"{nbd_type}_coverage"

    # Hit Count Summary (per kb)
    logger.info(f"\n{nbd_label} Hit Count per kb Summary:")
    hit_summary = hits.groupby("coding_class")[hit_count_col].describe()
    logger.info(f"\n{hit_summary}")

    # Max NBD Length Summary
    logger.info(f"\n{nbd_label} Max Length as % of transcript Summary:")
    max_length_summary = hits.groupby("coding_class")[max_length_col].describe()
    logger.info(f"\n{max_length_summary}")

    # Total NBD Length Summary
    # logger.info(f"\n{nbd_label} Total Overlap as % Summary:")
    # total_length_summary = hits.groupby('coding_class')[total_length_col].describe()
    # logger.info(f"\n{total_length_summary}")

    # Coverage Summary (already a fraction)
    logger.info(f"\n{nbd_label} Coverage Summary:")
    coverage_summary = hits.groupby("coding_class")[coverage_col].describe()
    logger.info(f"\n{coverage_summary}")

    logger.info("=" * 60)


def create_visualization(hits: pd.DataFrame, nbd_type: str, output_dir: str) -> None:
    """
    Create multi-panel visualization of NBD metrics using KDE plots.

    Args:
        hits: DataFrame with aggregated transcript metrics
        nbd_type: Type of NBD feature being analyzed
        output_dir: Prefix for output file path
    """
    logger.info(f"Creating KDE visualization plots for {nbd_type.upper()}")

    if hits.empty:
        logger.warning(
            f"No data to plot for {nbd_type.upper()}, skipping visualization"
        )
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"{nbd_type.upper()} — No data", fontsize=14)
        plt.tight_layout()
        plt.savefig(
            f"{output_dir}/{nbd_type}_distributions.png", dpi=150, bbox_inches="tight"
        )
        plt.close()
        return

    # Dynamic column names (now normalized per kb)
    hit_count_col = f"{nbd_type}_hit_count_per_kb"
    max_length_col = f"max_{nbd_type}_length_as_%"
    coverage_col = f"{nbd_type}_coverage"
    nbd_label = nbd_type.upper()

    # Set seaborn style for better aesthetics
    sns.set_style("whitegrid")

    # Define color palette for consistency
    palette = {"coding": "#1f77b4", "lncRNA": "#ff7f0e"}

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Hit Count Distribution (KDE)
    logger.debug(f"Plotting {nbd_label} hit count KDE distribution")
    for coding_class in sorted(hits["coding_class"].unique()):
        data = hits[hits["coding_class"] == coding_class][hit_count_col]
        sns.kdeplot(
            data=data,
            ax=axes[0, 0],
            label=coding_class,
            color=palette.get(coding_class),
            fill=True,
            alpha=0.3,
            linewidth=2,
        )
    # Set x-limits to 5th and 95th percentiles
    all_data = hits[hit_count_col]
    axes[0, 0].set_xlim(-all_data.quantile(0.05), all_data.quantile(0.99))
    axes[0, 0].set_title(
        f"{nbd_label} Hit Count per kb", fontsize=12, fontweight="bold"
    )
    axes[0, 0].set_xlabel(f"{nbd_label} Hit Count / kb", fontsize=10)
    axes[0, 0].set_ylabel("Density", fontsize=10)
    axes[0, 0].legend(title="Coding Class", fontsize=9)
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: Max NBD Length Distribution (KDE)
    logger.debug(f"Plotting max {nbd_label} length KDE distribution")
    for coding_class in sorted(hits["coding_class"].unique()):
        data = hits[hits["coding_class"] == coding_class][max_length_col]
        sns.kdeplot(
            data=data,
            ax=axes[0, 1],
            label=coding_class,
            color=palette.get(coding_class),
            fill=True,
            alpha=0.3,
            linewidth=2,
        )
    # Set x-limits to 5th and 95th percentiles
    all_data = hits[max_length_col]
    axes[0, 1].set_xlim(-all_data.quantile(0.05), all_data.quantile(0.99))
    axes[0, 1].set_title(
        f"Max {nbd_label} Length as coverage of transcript",
        fontsize=12,
        fontweight="bold",
    )
    axes[0, 1].set_xlabel(f"Max {nbd_label} Length coverage (%)", fontsize=10)
    axes[0, 1].set_ylabel("Density", fontsize=10)
    axes[0, 1].legend(title="Coding Class", fontsize=9)
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Coverage Distribution (KDE)
    logger.debug(f"Plotting {nbd_label} coverage KDE distribution")
    for coding_class in sorted(hits["coding_class"].unique()):
        data = hits[hits["coding_class"] == coding_class][coverage_col]
        sns.kdeplot(
            data=data,
            ax=axes[1, 0],
            label=coding_class,
            color=palette.get(coding_class),
            fill=True,
            alpha=0.3,
            linewidth=2,
        )
    # Set x-limits to 5th and 95th percentiles
    all_data = hits[coverage_col]
    axes[1, 0].set_xlim(-all_data.quantile(0.05), all_data.quantile(0.99))
    axes[1, 0].set_title(f"{nbd_label} Coverage", fontsize=12, fontweight="bold")
    axes[1, 0].set_xlabel(f"{nbd_label} Coverage (%)", fontsize=10)
    axes[1, 0].set_ylabel("Density", fontsize=10)
    axes[1, 0].legend(title="Coding Class", fontsize=9)
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Transcript Length Distribution (KDE)
    logger.debug("Plotting transcript length KDE distribution")
    for coding_class in sorted(hits["coding_class"].unique()):
        data = hits[hits["coding_class"] == coding_class]["transcript_length"]
        sns.kdeplot(
            data=data,
            ax=axes[1, 1],
            label=coding_class,
            color=palette.get(coding_class),
            fill=True,
            alpha=0.3,
            linewidth=2,
        )
    # Set x-limits to 5th and 95th percentiles
    all_data = hits["transcript_length"]
    axes[1, 1].set_xlim(-all_data.quantile(0.05), all_data.quantile(0.99))
    axes[1, 1].set_title(
        "Transcript Length Distribution", fontsize=12, fontweight="bold"
    )
    axes[1, 1].set_xlabel("Transcript Length (bp)", fontsize=10)
    axes[1, 1].set_ylabel("Density", fontsize=10)
    axes[1, 1].legend(title="Coding Class", fontsize=9)
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    plot_file = f"{output_dir}/{nbd_type}_distributions.png"
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    logger.info(f"Saved visualization to: {plot_file}")

    plt.close()


def save_transcript_summary(hits: pd.DataFrame, nbd_type: str, output_dir: str) -> None:
    """
    Save per-transcript summary statistics to file.

    Args:
        hits: DataFrame with aggregated transcript metrics
        nbd_type: Type of NBD feature being analyzed
        output_dir: Prefix for output file path
    """
    summary_file = f"{output_dir}/transcript_{nbd_type}_summary.tsv"
    logger.info(f"Writing per-transcript {nbd_type.upper()} summary to: {summary_file}")

    hits.to_csv(summary_file, sep="\t", index=True)
    logger.info(f"Wrote {len(hits)} transcript records to {summary_file}")


def main(
    overlap_file: str,
    output_dir: str,
    pc_id_file: str,
    lnc_id_file: str,
    nbd_type: str = "g4",
    keep_other: bool = False,
) -> None:
    """
    Main analysis pipeline for NBD overlaps with transcripts.

    Args:
        overlap_file: Path to overlap BED file
        output_dir: Prefix for output files
        pc_id_file: Path to protein-coding transcript IDs
        lnc_id_file: Path to lncRNA transcript IDs
        nbd_type: Type of NBD feature (e.g., 'g4', 'z_dna', 'triplex')
    """
    logger.info("=" * 60)
    logger.info(f"STARTING {nbd_type.upper()} TRANSCRIPT ANALYSIS PIPELINE")
    logger.info("=" * 60)
    logger.info(f"NBD Type: {nbd_type}")

    # Load transcript classification IDs
    pc_ids, lnc_ids = load_transcript_ids(pc_id_file, lnc_id_file)

    # Read and process overlap file
    df = read_overlap_file(overlap_file)
    df = filter_features(df)
    df = extract_attributes(df)

    # Classify and calculate metrics
    df = classify_transcripts(df, pc_ids, lnc_ids, version=False, keep_other=keep_other)
    df = calculate_lengths(df)

    # Aggregate per transcript
    hits = aggregate_transcript_metrics(df, nbd_type)

    # Generate outputs
    print_summary_statistics(hits, nbd_type)
    create_visualization(hits, nbd_type, output_dir)
    save_transcript_summary(hits, nbd_type, output_dir)

    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETED SUCCESSFULLY")
    logger.info("=" * 60)
    log_job_completion(logger)


if __name__ == "__main__":
    # Get parameters from Snakemake with error handling
    try:
        overlap_file = snakemake.input.get("overlap")
        if not overlap_file:
            raise ValueError("Required input 'overlap' not provided in snakemake.input")
        if not Path(overlap_file).exists():
            raise FileNotFoundError(f"Overlap file not found: {overlap_file}")

        pc_id_file = snakemake.input.get("pc_id_file")
        if not pc_id_file:
            raise ValueError(
                "Required input 'pc_id_file' not provided in snakemake.input"
            )
        if not Path(pc_id_file).exists():
            raise FileNotFoundError(f"Protein-coding ID file not found: {pc_id_file}")

        lnc_id_file = snakemake.input.get("lnc_id_file")
        if not lnc_id_file:
            raise ValueError(
                "Required input 'lnc_id_file' not provided in snakemake.input"
            )
        if not Path(lnc_id_file).exists():
            raise FileNotFoundError(f"lncRNA ID file not found: {lnc_id_file}")

        output_dir = snakemake.params.get("output_dir", "transcript_nbd_summary")
        logger.info(f"Output directory parameter: {output_dir}")

        # Get NBD type from wildcards or params
        nbd_type = snakemake.wildcards.get("nbd_type") or snakemake.params.get(
            "nbd_type"
        )
        if not nbd_type:
            logger.error(
                "Required parameter 'nbd_type' not provided in snakemake.wildcards or snakemake.params"
            )
            sys.exit(1)
        logger.info(f"NBD type parameter: {nbd_type}")

        # Whether to keep non-pc and non-lncRNA transcripts in the analysis
        keep_other = snakemake.params.get("keep_other", False)

        # Ensure output directory exists
        output_dir = Path(output_dir)
        if output_dir and not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created output directory: {output_dir}")

    except (ValueError, FileNotFoundError) as e:
        logger.error(f"Input validation error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during input validation: {e}")
        sys.exit(1)

    main(overlap_file, output_dir, pc_id_file, lnc_id_file, nbd_type, keep_other)

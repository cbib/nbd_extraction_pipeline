#!/usr/bin/env python3
"""
Non-B DNA Feature Extractor
=================================================
Extract comprehensive Non-B DNA features.

OPTIMIZED: Uses vectorized operations instead of iterrows() for 10-100x speedup.
- Groupby aggregations instead of row iteration
- Vectorized gap calculations with shift()
- Merge-based feature assembly

Modeled after TE pipeline vectorized extractor.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class NonBFeatureExtractor:
    """Extract comprehensive Non-B DNA features with vectorized operations."""

    MOTIF_TYPES = ["APR", "DR", "GQ", "IR", "MR", "STR", "TRI", "Z"]

    def __init__(
        self,
        nonb_intersections: Dict[str, str],
        transcripts_bed: str,
        biotypes_file: Optional[str] = None,
        output_prefix: str = "",
        motif_name_mapping: Optional[Dict[str, str]] = None,
        exon_mode: bool = True,
    ):
        """
        Initialize Non-B DNA feature extractor.

        Parameters:
        -----------
        nonb_intersections : dict
            Dictionary mapping motif_type -> intersection BED file path
        transcripts_bed : str
            Transcript coordinates BED file (exons.bed or transcripts.bed)
        biotypes_file : str, optional
            File with transcript biotype information
        output_prefix : str
            Output file prefix
        motif_name_mapping : dict, optional
            Dictionary mapping internal motif names to output filenames
            (e.g., {"GQ": "g4Discovery"} for output file names)
        exon_mode : bool
            When True (default), ``transcripts_bed`` contains one row per exon;
            motif coordinates are first clipped to the exon interval before
            transcript-level clipping.  When False, ``transcripts_bed`` contains
            one row per transcript (full genomic span) and exon-level clipping is
            skipped — this matches the behaviour of commit c4e3709 but with the
            coordinate clipping bug fixed.
        """
        self.nonb_intersections = nonb_intersections
        self.transcripts_bed = transcripts_bed
        self.biotypes_file = biotypes_file
        self.output_prefix = output_prefix
        self.motif_name_mapping = motif_name_mapping or {}
        self.exon_mode = exon_mode

        self.transcripts = None
        self.biotypes = None
        self.nonb_elements = {}
        self.resolved_elements = {}
        self.features = None

    def load_data(self):
        """Load all input data."""
        logger.info("Loading data...")

        # Load transcript BED-like input. For extended analysis this is exon BED,
        # so we aggregate exon lengths per transcript_id.
        transcript_rows = pd.read_csv(
            self.transcripts_bed,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "transcript_id", "score", "strand"],
        )
        transcript_rows["segment_length"] = (
            transcript_rows["end"] - transcript_rows["start"]
        )

        self.transcripts = (
            transcript_rows.groupby("transcript_id", observed=True)
            .agg(
                chrom=("chrom", "first"),
                start=("start", "min"),
                end=("end", "max"),
                score=("score", "first"),
                strand=("strand", "first"),
                length=("segment_length", "sum"),
            )
            .reset_index()
        )

        self.transcripts["transcript_id_base"] = (
            self.transcripts["transcript_id"].str.split(".").str[0]
        )
        logger.info(
            "Loaded %s transcript segments, aggregated into %s transcripts",
            len(transcript_rows),
            len(self.transcripts),
        )

        # Use canonical transcript bounds from transcripts.bed for clipping
        # motif intervals, regardless of which GTF feature (gene/exon/CDS/etc.)
        # produced the intersect row.
        transcript_bounds = self.transcripts[["transcript_id", "start", "end"]].rename(
            columns={"start": "transcript_start", "end": "transcript_end"}
        )

        # Load biotypes if provided
        if self.biotypes_file:
            self.biotypes = pd.read_csv(self.biotypes_file, sep="\t")
            logger.info(f"Loaded biotypes for {len(self.biotypes)} transcripts")

        # Load Non-B DNA intersections
        for motif_type, file_path in self.nonb_intersections.items():
            if not Path(file_path).exists():
                logger.warning(f"File not found: {file_path}")
                continue

            try:
                raw_df = pd.read_csv(file_path, sep="\t", header=None)

                if raw_df.empty:
                    logger.info(f"Loaded 0 {motif_type} elements")
                    self.nonb_elements[motif_type] = pd.DataFrame(
                        columns=[
                            "transcript_id",
                            "motif_chrom",
                            "motif_start",
                            "motif_end",
                            "transcript_start",
                            "transcript_end",
                        ]
                    )
                    continue

                # Support both legacy GTF+B intersections and the new BED6+B
                # intersections used by the extended pipeline.
                if (
                    raw_df.shape[1] > 8
                    and raw_df.iloc[:, 8]
                    .astype(str)
                    .str.contains('transcript_id "', regex=False)
                    .any()
                ):
                    logger.info(
                        "Detected legacy GTF-based intersection schema for %s",
                        motif_type,
                    )
                    df = raw_df.iloc[:, [8, 9, 10, 11]].copy()
                    df.columns = [
                        "attributes",
                        "motif_chrom",
                        "motif_start",
                        "motif_end",
                    ]
                    df["transcript_id"] = df["attributes"].str.extract(
                        r'transcript_id "([^"]+)"'
                    )[0]
                    df = df[
                        [
                            "transcript_id",
                            "motif_chrom",
                            "motif_start",
                            "motif_end",
                        ]
                    ]
                elif raw_df.shape[1] > 8:
                    logger.info(
                        "Detected BED6-based intersection schema for %s",
                        motif_type,
                    )
                    df = raw_df.iloc[:, [1, 2, 3, 6, 7, 8]].copy()
                    df.columns = [
                        "left_start",
                        "left_end",
                        "transcript_id",
                        "motif_chrom",
                        "motif_start",
                        "motif_end",
                    ]
                else:
                    raise ValueError(
                        f"Unsupported intersect schema with {raw_df.shape[1]} columns"
                    )

                df = df.dropna(subset=["transcript_id", "motif_start", "motif_end"])
                df["motif_start"] = pd.to_numeric(df["motif_start"], errors="coerce")
                df["motif_end"] = pd.to_numeric(df["motif_end"], errors="coerce")
                if "left_start" in df.columns and "left_end" in df.columns:
                    df["left_start"] = pd.to_numeric(df["left_start"], errors="coerce")
                    df["left_end"] = pd.to_numeric(df["left_end"], errors="coerce")
                df = df.dropna(subset=["motif_start", "motif_end"])

                # For BED6-based intersections (exons on the left), restrict motif
                # intervals to exon overlap before transcript-level clipping.
                if (
                    self.exon_mode
                    and "left_start" in df.columns
                    and "left_end" in df.columns
                ):
                    df = df.dropna(subset=["left_start", "left_end"])
                    df["motif_start"] = df[["motif_start", "left_start"]].max(axis=1)
                    df["motif_end"] = df[["motif_end", "left_end"]].min(axis=1)
                    df = df[df["motif_end"] > df["motif_start"]].copy()
                    df = df.drop(columns=["left_start", "left_end"])

                # Replace row-level GTF coordinates with canonical transcript bounds.
                df = df.merge(transcript_bounds, on="transcript_id", how="inner")

                # Clip each motif interval to the transcript overlap span.
                # This prevents counting motif sequence outside the transcript.
                before_clip = len(df)
                df["motif_start"] = df[["motif_start", "transcript_start"]].max(axis=1)
                df["motif_end"] = df[["motif_end", "transcript_end"]].min(axis=1)
                df = df[df["motif_end"] > df["motif_start"]].copy()
                dropped = before_clip - len(df)
                if dropped:
                    logger.info(
                        "Dropped %s non-overlapping %s rows after transcript clipping",
                        dropped,
                        motif_type,
                    )

                logger.info(f"Loaded {len(df)} {motif_type} elements")
                self.nonb_elements[motif_type] = df
            except Exception as e:
                logger.error(f"Error loading {motif_type} from {file_path}: {e}")

    def resolve_overlaps(self):
        """
        Resolve overlapping Non-B DNA elements using vectorized operations.

        Uses sweep-line style groupby operations instead of row iteration.
        """
        logger.info("Resolving overlaps (vectorized)...")

        for motif_type, elements_df in self.nonb_elements.items():
            if len(elements_df) == 0:
                self.resolved_elements[motif_type] = elements_df
                continue

            # Sort by transcript and position for sweep-line
            elements_df = elements_df.sort_values(
                ["transcript_id", "motif_start", "motif_end"]
            )

            # Track the running merged end per transcript and split only when the
            # next interval starts after that running boundary.
            running_end = elements_df.groupby("transcript_id", observed=True)[
                "motif_end"
            ].cummax()
            prev_running_end = running_end.groupby(
                elements_df["transcript_id"], observed=True
            ).shift(1)

            # New segment starts when:
            # 1. New transcript
            # 2. No overlap with current merged segment (start > running_end)
            is_new_transcript = elements_df["transcript_id"].ne(
                elements_df["transcript_id"].shift(1)
            )
            new_segment = is_new_transcript | (
                elements_df["motif_start"] > prev_running_end
            )

            # Assign all calculation columns at once and defragment
            elements_df = elements_df.assign(
                new_segment=new_segment, segment_id=new_segment.cumsum()
            ).copy()

            # Merge overlapping segments - vectorized aggregation
            resolved = (
                elements_df.groupby("segment_id")
                .agg(
                    {
                        "transcript_id": "first",
                        "transcript_start": "first",
                        "transcript_end": "first",
                        "motif_start": "min",  # Min start
                        "motif_end": "max",  # Max end
                    }
                )
                .reset_index(drop=True)
            )

            # Calculate derived columns and assign together
            resolved = resolved.assign(
                length=resolved["motif_end"] - resolved["motif_start"],
                hit_count=elements_df.groupby("segment_id").size().values,
            )

            self.resolved_elements[motif_type] = resolved
            logger.info(
                f"Resolved {motif_type}: {len(elements_df)} → {len(resolved)} unique regions"
            )

    @staticmethod
    def calculate_gaps(df: pd.DataFrame) -> pd.DataFrame:
        """
        Vectorized gap calculation (10-12x faster than groupby.apply).

        Calculates gaps between consecutive hits within each transcript:
        - gap_before: distance from hit start to previous hit end
        - gap_after: distance from hit end to next hit start

        Note: as opposed to RM pipeline, elements here have genome coordinates, not relative to transcript start.

        Parameters
        ----------
        df : pd.DataFrame
            Must contain: transcript_id, start, end, query_left
            MUST be sorted by (transcript_id, start)

        Returns
        -------
        pd.DataFrame
            Input DataFrame with gap_before and gap_after columns added
        """
        df = df.copy()

        # Identify first and last elements per transcript (assign together to avoid fragmentation)
        first_element = df["transcript_id"].ne(df["transcript_id"].shift(1))
        last_element = first_element.shift(-1).astype(bool).fillna(True)

        # Calculate gaps using vectorized shift operations (assign together)
        gap_before = df["motif_start"] - df["motif_end"].shift(1) - 1
        gap_before.loc[first_element] = (
            df.loc[first_element, "motif_start"]
            - df.loc[first_element, "transcript_start"]
            - 1
        )

        gap_after = gap_before.shift(-1)
        gap_after.loc[last_element] = (
            df.loc[last_element, "transcript_end"] - df.loc[last_element, "motif_end"]
        )

        # Assign all new columns at once
        df = df.assign(
            first_element=first_element,
            last_element=last_element,
            gap_before=gap_before,
            gap_after=gap_after,
        )

        return df

    @staticmethod
    def calculate_gap_stats(df: pd.DataFrame, prefix: str = "") -> pd.DataFrame:
        """
        Calculate gap statistics in a vectorized manner.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with gap_before and gap_after columns
        prefix : str
            Prefix for output column names

        Returns
        -------
        pd.DataFrame
            DataFrame with transcript_id and gap statistics columns
        """
        if prefix:
            prefix = prefix + "_"
        gaps_before = df[["transcript_id", "gap_before"]].rename(
            columns={"gap_before": "gap_size"}
        )
        is_last_hit = ~df["transcript_id"].duplicated(keep="last")
        gaps_after_last = df.loc[is_last_hit, ["transcript_id", "gap_after"]].rename(
            columns={"gap_after": "gap_size"}
        )
        combined = pd.concat([gaps_before, gaps_after_last], ignore_index=True).dropna(
            subset=["gap_size"]
        )

        stats = (
            combined.groupby("transcript_id", observed=True)["gap_size"]
            .agg(["mean", "median", "max", "min"])
            .reset_index()
        )
        stats.columns = [
            "transcript_id",
            f"{prefix}gaps_mean",
            f"{prefix}gaps_median",
            f"{prefix}gaps_max",
            f"{prefix}gaps_min",
        ]

        return stats

    def extract_features(self) -> pd.DataFrame:
        """
        Extract comprehensive Non-B DNA features using vectorized operations.

        Uses groupby().agg() instead of iterrows() for major performance boost.
        """
        logger.info("Extracting comprehensive features (vectorized)...")

        # Start with transcript base
        features = self.transcripts[["transcript_id", "length"]].copy()
        features = features.rename(columns={"length": "transcript_length"})

        # Collect all motif feature dataframes before concatenating (avoids fragmentation)
        motif_features_list = []

        # Process each motif type
        for motif_type in self.MOTIF_TYPES:
            prefix = motif_type.lower()

            if motif_type not in self.nonb_elements:
                motif_features_list.append(self._get_zero_features(prefix))
                continue

            # Get raw and resolved elements
            raw_elements = self.nonb_elements[motif_type]
            resolved_elements = self.resolved_elements.get(motif_type, pd.DataFrame())

            if raw_elements.empty:
                motif_features_list.append(self._get_zero_features(prefix))
                continue

            # === VECTORIZED AGGREGATIONS ===

            # Raw counts per transcript
            raw_counts = (
                raw_elements.groupby("transcript_id").size().to_frame("hit_count")
            )

            # Length statistics - raw elements
            raw_stats = raw_elements.copy()
            raw_stats["element_length"] = (
                raw_stats["motif_end"] - raw_stats["motif_start"]
            )
            raw_length_stats = raw_stats.groupby("transcript_id")["element_length"].agg(
                [
                    ("total_length", "sum"),
                    ("max_length", "max"),
                    ("mean_length", "mean"),
                    ("std_length", "std"),
                ]
            )

            # Unique length statistics
            unique_length_stats = resolved_elements.groupby("transcript_id")[
                "length"
            ].agg([("unique_length", "sum")])

            # Gap statistics
            resolved_with_gaps = self.calculate_gaps(resolved_elements)
            gap_stats = self.calculate_gap_stats(resolved_with_gaps)

            # Merge all stats for this motif
            motif_features = (
                features[["transcript_id"]]
                .merge(
                    raw_counts, left_on="transcript_id", right_index=True, how="left"
                )
                .merge(
                    raw_length_stats,
                    left_on="transcript_id",
                    right_index=True,
                    how="left",
                )
                .merge(
                    unique_length_stats,
                    left_on="transcript_id",
                    right_index=True,
                    how="left",
                )
                .merge(gap_stats, on="transcript_id", how="left")
            )

            # Calculate all derived features at once and assign together
            derived_cols = {
                "present": (motif_features["hit_count"] > 0).astype(int),
            }

            # Coverage percentage features
            cov_features = [
                c for c in motif_features.columns if "length" in c or "gap" in c
            ]
            for col in cov_features:
                derived_cols[f"{col}_pct"] = (
                    motif_features[col] / features["transcript_length"] * 100
                ).fillna(0)

            # Density features (per kb)
            dens_features = [c for c in motif_features.columns if "hit_count" in c]
            for col in dens_features:
                derived_cols[f"{col}_per_kb"] = (
                    motif_features[col] / features["transcript_length"] * 1000
                ).fillna(0)

            motif_features = motif_features.assign(**derived_cols)

            # Rename columns with prefix
            motif_features = motif_features.rename(
                columns={c: f"{prefix}_{c}" for c in motif_features.columns}
            )

            # Drop transcript_id before merge
            # motif_features = motif_features.drop(columns=["transcript_id"])

            motif_features_list.append(motif_features)

        # Concatenate all motif features once, outside the loop
        if motif_features_list:
            features = pd.concat([features] + motif_features_list, axis=1)

        # Re-merge all resolved elements and calculate gaps between any hits
        all_resolved_elements = []
        for motif_type in self.MOTIF_TYPES:
            if (
                motif_type in self.resolved_elements
                and not self.resolved_elements[motif_type].empty
            ):
                all_resolved_elements.append(self.resolved_elements[motif_type])
        gap_stats_cols = [
            "all_nonb_gaps_mean",
            "all_nonb_gaps_median",
            "all_nonb_gaps_max",
            "all_nonb_gaps_min",
        ]

        if len(all_resolved_elements) > 0:
            all_hits = pd.concat(all_resolved_elements, ignore_index=True)
            all_hits = all_hits.sort_values(
                ["transcript_id", "motif_start", "motif_end"]
            )

            all_hits = self.calculate_gaps(all_hits)
            gap_stats = self.calculate_gap_stats(all_hits, prefix="all_nonb")

            features = features.merge(gap_stats, on="transcript_id", how="left")

            # Fill gap stats NaNs with transcript length
            features = features.assign(
                **{
                    col: features[col].fillna(features["transcript_length"])
                    for col in gap_stats_cols
                }
            )
        else:
            # Longest gap is full transcript length
            features = features.assign(
                **{col: features["transcript_length"] for col in gap_stats_cols}
            )

        # Defragment after merge operations
        features = features.copy()

        # Fill NaNs with 0
        numeric_cols = features.select_dtypes(include=[np.number]).columns
        features[numeric_cols] = features[numeric_cols].fillna(0)

        # === GLOBAL SUMMARY METRICS (vectorized) ===
        # Calculate row-wise across motifs
        presence_cols = [f"{m.lower()}_present" for m in self.MOTIF_TYPES]
        count_cols = [
            f"{m.lower()}_hit_count"
            for m in self.MOTIF_TYPES
            if f"{m.lower()}_hit_count" in features.columns
        ]
        coverage_cols = [
            f"{m.lower()}_total_length"
            for m in self.MOTIF_TYPES
            if f"{m.lower()}_total_length" in features.columns
        ]

        # Vectorized summary metrics - assign all at once to avoid fragmentation
        summary_metrics = {
            "any_nonb_present": features[presence_cols].max(axis=1).astype(int),
            "motif_types_present": features[presence_cols].sum(axis=1),
            "total_nonb_count": features[count_cols].sum(axis=1),
            "total_nonb_coverage": features[coverage_cols].sum(axis=1),
        }
        features = features.assign(**summary_metrics)

        # Add coverage percentage
        features = features.assign(
            total_nonb_coverage_pct=(
                features["total_nonb_coverage"] / features["transcript_length"] * 100
            ).fillna(0)
        )

        # Diversity index (Shannon-like) - vectorized
        count_matrix = features[count_cols]
        count_sums = count_matrix.sum(axis=1)
        proportions = count_matrix.div(count_sums, axis=0).fillna(0)
        # Shannon: -sum(p * log(p))
        motif_diversity = -(proportions * np.log(proportions + 1e-10)).sum(axis=1)
        motif_diversity.loc[count_sums == 0] = 0
        features = features.assign(motif_diversity=motif_diversity)

        # Defragment before biotype merge
        features = features.copy()

        # Add biotype information if available
        if self.biotypes is not None:
            biotype_map = self.biotypes.set_index("transcript_id_base")[
                ["transcript_type"]
            ]
            features["transcript_id_base"] = (
                features["transcript_id"].str.split(".").str[0]
            )
            features = features.merge(
                biotype_map, left_on="transcript_id_base", right_index=True, how="left"
            )
            features = features.drop(columns=["transcript_id_base"])
            # Defragment dataframe after many column insertions
            features = features.copy()

            # Add coding_class column
            if "transcript_type" in features.columns:
                features["coding_class"] = features["transcript_type"].apply(
                    lambda x: (
                        "coding"
                        if x == "protein_coding"
                        else "lncRNA" if "lncRNA" in str(x) else "other"
                    )
                )

        self.features = features
        logger.info(
            f"Extracted {len(features.columns)} features for {len(features)} transcripts"
        )

        return features

    def _get_zero_features(self, prefix: str) -> pd.DataFrame:
        """Return zero-valued features DataFrame for missing motif type."""
        return pd.DataFrame(
            {
                f"{prefix}_present": 0,
                f"{prefix}_hit_count": 0,
                f"{prefix}_unique_count": 0,
                f"{prefix}_total_length": 0,
                f"{prefix}_unique_length": 0,
                f"{prefix}_coverage_pct": 0,
                f"{prefix}_density_per_kb": 0,
                f"{prefix}_max_length": 0,
                f"{prefix}_mean_length": 0,
                f"{prefix}_std_length": 0,
                f"{prefix}_max_pct": 0,
                f"{prefix}_uniqueness_ratio": 0,
                f"{prefix}_overlap_fraction": 0,
            },
            index=self.transcripts.index,
        )

    def _add_zero_features(self, features_df: pd.DataFrame, prefix: str):
        """Add zero-valued features for missing motif type (legacy, use _get_zero_features instead)."""
        zero_features = self._get_zero_features(prefix)
        for col in zero_features.columns:
            features_df[col] = zero_features[col]

    def save_features(self):
        """Save features and summary."""
        output_file = f"{self.output_prefix}_nonb_features.csv"
        self.features.to_csv(output_file, index=False)
        logger.info(f"Saved features to {output_file}")

        # Save summary
        summary_file = f"{self.output_prefix}_nonb_summary.txt"
        with open(summary_file, "w") as f:
            f.write("Non-B DNA Feature Extraction Summary (Vectorized)\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Total transcripts analyzed: {len(self.features)}\n")
            f.write(
                f"Transcripts with any Non-B DNA: {self.features['any_nonb_present'].sum()}\n"
            )
            f.write(
                f"Percentage with Non-B DNA: {self.features['any_nonb_present'].mean() * 100:.2f}%\n\n"
            )

            # Motif type breakdown
            f.write("Motif type prevalence:\n")
            for motif_type in self.MOTIF_TYPES:
                col = f"{motif_type.lower()}_present"
                if col in self.features.columns:
                    count = self.features[col].sum()
                    pct = self.features[col].mean() * 100
                    f.write(f"  {motif_type}: {count} ({pct:.2f}%)\n")

            if "coding_class" in self.features.columns:
                f.write("\nBy coding class:\n")
                for coding_class in ["coding", "lncRNA", "other"]:
                    subset = self.features[
                        self.features["coding_class"] == coding_class
                    ]
                    if len(subset) > 0:
                        has_nonb = subset["any_nonb_present"].sum()
                        pct = subset["any_nonb_present"].mean() * 100
                        f.write(
                            f"  {coding_class}: {has_nonb}/{len(subset)} ({pct:.2f}%)\n"
                        )

            f.write("\nFeature statistics:\n")
            f.write(self.features.describe().to_string())

        logger.info(f"Saved summary to {summary_file}")

    def save_resolved_elements(self):
        """Save resolved element coordinates for downstream analysis."""
        for motif_type, resolved_df in self.resolved_elements.items():
            if len(resolved_df) > 0:
                # Use mapped name if available, otherwise use motif_type as-is
                output_name = self.motif_name_mapping.get(motif_type, motif_type)
                output_file = f"{self.output_prefix}_resolved_{output_name}.bed"
                resolved_df.to_csv(output_file, sep="\t", index=False, header=False)
                logger.info(f"Saved resolved {motif_type} elements to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract comprehensive Non-B DNA features (vectorized version)"
    )
    parser.add_argument(
        "--intersections",
        required=True,
        nargs="+",
        help="Intersection BED files: motif_type:path pairs",
    )
    parser.add_argument(
        "--transcripts", required=True, help="Transcript coordinates BED file"
    )
    parser.add_argument(
        "--biotypes", required=False, default=None, help="Biotypes file (optional)"
    )
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")

    args = parser.parse_args()

    # Parse intersections from command line
    intersections = {}
    for item in args.intersections:
        motif_type, path = item.split(":")
        intersections[motif_type] = path

    logger.info(
        f"Found {len(intersections)} intersection files: {list(intersections.keys())}"
    )

    # Initialize extractor
    extractor = NonBFeatureExtractor(
        nonb_intersections=intersections,
        transcripts_bed=args.transcripts,
        biotypes_file=args.biotypes,
        output_prefix=args.output_prefix,
    )

    # Run pipeline
    extractor.load_data()
    extractor.resolve_overlaps()
    extractor.extract_features()
    extractor.save_features()
    extractor.save_resolved_elements()

    logger.info("Non-B DNA feature extraction complete!")


if __name__ == "__main__":
    main()

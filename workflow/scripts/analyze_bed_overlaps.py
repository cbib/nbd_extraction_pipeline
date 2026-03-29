#!/usr/bin/env python3
"""
Analyze overlaps within a single BED file and report informative statistics.

This preprocessing tool helps understand the nature and extent of overlapping
regions in Non-B DNA motif annotations, informing downstream overlap resolution
strategies.

Usage:
    python analyze_bed_overlaps.py input.bed [output_report.txt]

Output includes:
    - Summary statistics (total regions, overlap counts, percentages)
    - Overlap distribution (histogram of regions per overlap cluster)
    - Size distributions (for all regions, overlapping vs non-overlapping)
    - Detailed overlap clusters (groups of mutually overlapping regions)
"""

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def read_bed(filepath):
    """Read BED file into pandas DataFrame."""
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
        usecols=[0, 1, 2, 3, 4, 5] if count_columns(filepath) >= 6 else [0, 1, 2],
    )

    # Add columns if not present
    if "name" not in df.columns:
        df["name"] = df.index.astype(str)
    if "score" not in df.columns:
        df["score"] = 0
    if "strand" not in df.columns:
        df["strand"] = "."

    df["length"] = df["end"] - df["start"]
    df["region_id"] = df.index
    return df


def count_columns(filepath):
    """Count columns in BED file."""
    with open(filepath) as f:
        first_line = f.readline().strip()
        return len(first_line.split("\t"))


def find_overlaps_vectorized(df):
    """
    Find overlapping regions using vectorized operations.

    Returns:
        - overlap_matrix: Dict mapping region_id -> set of overlapping region_ids
        - clusters: List of sets, each containing mutually overlapping region_ids
    """
    # Sort by chromosome and start position
    df_sorted = df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

    overlap_graph = defaultdict(set)

    # Process each chromosome separately
    for chrom in df_sorted["chrom"].unique():
        chrom_df = df_sorted[df_sorted["chrom"] == chrom].copy()

        if len(chrom_df) <= 1:
            continue

        # For each region, find overlaps with subsequent regions
        for i, (idx, row) in enumerate(chrom_df.iterrows()):
            region_id = row["region_id"]
            chrom_start = row["start"]
            chrom_end = row["end"]

            # Check subsequent regions until we're past possible overlaps
            for j in range(i + 1, len(chrom_df)):
                other_row = chrom_df.iloc[j]
                other_id = other_row["region_id"]
                other_start = other_row["start"]
                other_end = other_row["end"]

                # If other region starts after current ends, no more overlaps possible
                if other_start >= chrom_end:
                    break

                # Check for overlap
                if other_start < chrom_end and other_end > chrom_start:
                    overlap_graph[region_id].add(other_id)
                    overlap_graph[other_id].add(region_id)

    # Find connected components (overlap clusters)
    clusters = []
    visited = set()

    for region_id in df["region_id"]:
        if region_id in visited:
            continue

        # BFS to find all connected regions
        cluster = set()
        queue = [region_id]

        while queue:
            current = queue.pop(0)
            if current in visited:
                continue

            visited.add(current)
            cluster.add(current)

            # Add neighbors
            for neighbor in overlap_graph.get(current, []):
                if neighbor not in visited:
                    queue.append(neighbor)

        if len(cluster) > 1:
            clusters.append(cluster)
        elif len(cluster) == 1 and region_id not in overlap_graph:
            # Isolated region
            pass

    return overlap_graph, clusters


def calculate_overlap_stats(df, overlap_graph, clusters):
    """Calculate comprehensive overlap statistics."""
    stats = {}

    # Basic counts
    stats["total_regions"] = len(df)
    stats["regions_with_overlaps"] = len([r for r in overlap_graph if overlap_graph[r]])
    stats["isolated_regions"] = stats["total_regions"] - stats["regions_with_overlaps"]
    stats["overlap_percentage"] = (
        (stats["regions_with_overlaps"] / stats["total_regions"] * 100)
        if stats["total_regions"] > 0
        else 0
    )

    # Cluster statistics
    stats["num_clusters"] = len(clusters)
    stats["cluster_sizes"] = [len(c) for c in clusters]
    stats["max_cluster_size"] = (
        max(stats["cluster_sizes"]) if stats["cluster_sizes"] else 0
    )
    stats["mean_cluster_size"] = (
        np.mean(stats["cluster_sizes"]) if stats["cluster_sizes"] else 0
    )
    stats["median_cluster_size"] = (
        np.median(stats["cluster_sizes"]) if stats["cluster_sizes"] else 0
    )

    # Overlaps per region
    overlaps_per_region = [len(overlap_graph.get(r, [])) for r in df["region_id"]]
    stats["mean_overlaps_per_region"] = np.mean(overlaps_per_region)
    stats["median_overlaps_per_region"] = np.median(overlaps_per_region)
    stats["max_overlaps_per_region"] = max(overlaps_per_region)

    # Size distributions
    df["has_overlap"] = df["region_id"].isin(
        [r for r in overlap_graph if overlap_graph[r]]
    )

    stats["mean_length_all"] = df["length"].mean()
    stats["median_length_all"] = df["length"].median()
    stats["mean_length_overlapping"] = (
        df[df["has_overlap"]]["length"].mean()
        if stats["regions_with_overlaps"] > 0
        else 0
    )
    stats["mean_length_isolated"] = (
        df[~df["has_overlap"]]["length"].mean() if stats["isolated_regions"] > 0 else 0
    )

    return stats


def generate_report(filepath, df, overlap_graph, clusters, stats, output_file=None):
    """Generate comprehensive text report."""
    lines = []

    lines.append("=" * 80)
    lines.append(f"BED File Overlap Analysis Report")
    lines.append("=" * 80)
    lines.append(f"File: {filepath}")
    lines.append(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # Summary Statistics
    lines.append("SUMMARY STATISTICS")
    lines.append("-" * 80)
    lines.append(f"Total regions:                    {stats['total_regions']:,}")
    lines.append(
        f"Regions with overlaps:            {stats['regions_with_overlaps']:,} ({stats['overlap_percentage']:.2f}%)"
    )
    lines.append(f"Isolated regions:                 {stats['isolated_regions']:,}")
    lines.append(f"Number of overlap clusters:       {stats['num_clusters']:,}")
    lines.append("")

    # Overlap Characteristics
    lines.append("OVERLAP CHARACTERISTICS")
    lines.append("-" * 80)
    lines.append(
        f"Mean overlaps per region:         {stats['mean_overlaps_per_region']:.2f}"
    )
    lines.append(
        f"Median overlaps per region:       {stats['median_overlaps_per_region']:.1f}"
    )
    lines.append(
        f"Max overlaps per region:          {stats['max_overlaps_per_region']}"
    )
    lines.append("")

    # Cluster Statistics
    lines.append("CLUSTER STATISTICS")
    lines.append("-" * 80)
    lines.append(f"Largest cluster size:             {stats['max_cluster_size']}")
    lines.append(f"Mean cluster size:                {stats['mean_cluster_size']:.2f}")
    lines.append(
        f"Median cluster size:              {stats['median_cluster_size']:.1f}"
    )
    lines.append("")

    # Cluster size distribution
    if stats["cluster_sizes"]:
        lines.append("Cluster Size Distribution:")
        cluster_size_counts = (
            pd.Series(stats["cluster_sizes"]).value_counts().sort_index()
        )
        for size, count in cluster_size_counts.items():
            lines.append(f"  {size:3d} regions: {count:6,} clusters")
    lines.append("")

    # Region Length Statistics
    lines.append("REGION LENGTH STATISTICS")
    lines.append("-" * 80)
    lines.append(f"Mean length (all regions):        {stats['mean_length_all']:.2f} bp")
    lines.append(
        f"Median length (all regions):      {stats['median_length_all']:.2f} bp"
    )
    if stats["regions_with_overlaps"] > 0:
        lines.append(
            f"Mean length (overlapping):        {stats['mean_length_overlapping']:.2f} bp"
        )
    if stats["isolated_regions"] > 0:
        lines.append(
            f"Mean length (isolated):           {stats['mean_length_isolated']:.2f} bp"
        )
    lines.append("")

    # Largest clusters detail
    if clusters:
        lines.append("LARGEST OVERLAP CLUSTERS (Top 10)")
        lines.append("-" * 80)
        sorted_clusters = sorted(clusters, key=len, reverse=True)[:10]

        for i, cluster in enumerate(sorted_clusters, 1):
            cluster_regions = df[df["region_id"].isin(cluster)].sort_values("start")
            total_span = cluster_regions["end"].max() - cluster_regions["start"].min()
            total_bases = cluster_regions["length"].sum()
            chrom = cluster_regions["chrom"].iloc[0]
            start = cluster_regions["start"].min()
            end = cluster_regions["end"].max()

            lines.append(f"\nCluster {i}: {len(cluster)} overlapping regions")
            lines.append(f"  Location: {chrom}:{start}-{end}")
            lines.append(f"  Total span: {total_span:,} bp")
            lines.append(f"  Total bases (summed): {total_bases:,} bp")
            lines.append(f"  Coverage ratio: {total_bases/total_span:.2f}x")

    lines.append("")
    lines.append("=" * 80)

    report = "\n".join(lines)

    # Output
    if output_file:
        with open(output_file, "w") as f:
            f.write(report)
        print(f"Report written to: {output_file}")

    print(report)

    return report


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_bed = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    if not Path(input_bed).exists():
        print(f"Error: File not found: {input_bed}")
        sys.exit(1)

    print(f"Reading BED file: {input_bed}")
    df = read_bed(input_bed)

    print(f"Analyzing overlaps in {len(df):,} regions...")
    overlap_graph, clusters = find_overlaps_vectorized(df)

    print("Calculating statistics...")
    stats = calculate_overlap_stats(df, overlap_graph, clusters)

    print("\nGenerating report...")
    generate_report(input_bed, df, overlap_graph, clusters, stats, output_file)


if __name__ == "__main__":
    main()

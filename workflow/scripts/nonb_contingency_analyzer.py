#!/usr/bin/env python3
"""
Non-B DNA Contingency Analyzer
===============================
Perform chi-square contingency tests for Non-B DNA motif enrichment
differences between protein-coding and lncRNA transcripts.

Modeled after te_contingency_analyzer.py from the TE pipeline.
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class NonBContingencyAnalyzer:
    """Perform contingency analysis on Non-B DNA features."""

    def __init__(self, features_file: str, output_prefix: str):
        """
        Initialize contingency analyzer.

        Parameters:
        -----------
        features_file : str
            CSV file with Non-B DNA features
        output_prefix : str
            Output file prefix
        """
        self.features_file = features_file
        self.output_prefix = output_prefix
        self.features = None

        # Motif types to analyze
        self.motif_types = [
            "apr",
            "dr",
            "gq",
            "ir",
            "mr",
            "str",
            "tri",
            "z",
            "g4discovery",
            "slipped",
        ]

    def load_features(self):
        """Load feature matrix."""
        logger.info(f"Loading features from {self.features_file}")
        self.features = pd.read_csv(self.features_file)
        logger.info(
            f"Loaded {len(self.features)} transcripts with {len(self.features.columns)} features"
        )

        # Verify required columns
        if "coding_class" not in self.features.columns:
            raise ValueError("Feature file must contain 'coding_class' column")

    def perform_overall_test(self) -> Tuple[float, float, bool]:
        """
        Perform overall Non-B DNA presence test.

        Returns:
        --------
        chi2 : float
            Chi-square statistic
        p_value : float
            P-value
        significant : bool
            Whether result is significant (p < 0.05)
        """
        logger.info("Performing overall Non-B DNA presence test...")

        # Filter to coding and lncRNA only
        data = self.features[
            self.features["coding_class"].isin(["coding", "lncRNA"])
        ].copy()

        # Create contingency table
        contingency = pd.crosstab(data["coding_class"], data["any_nonb_present"])

        # Perform chi-square test
        chi2, p_value, dof, expected = chi2_contingency(contingency)

        significant = p_value < 0.05

        logger.info(
            f"Overall test - Chi2: {chi2:.4f}, p-value: {p_value:.4e}, Significant: {significant}"
        )

        return chi2, p_value, significant

    def perform_motif_type_tests(self) -> pd.DataFrame:
        """
        Perform chi-square tests for each motif type.

        Returns:
        --------
        pd.DataFrame
            Results with columns: motif_type, pct_coding, pct_lncrna,
            chi2, p_value, cramers_v, significant
        """
        logger.info("Performing motif type contingency tests...")

        # Filter to coding and lncRNA
        data = self.features[
            self.features["coding_class"].isin(["coding", "lncRNA"])
        ].copy()

        results = []

        for motif_type in self.motif_types:
            presence_col = f"{motif_type}_present"

            if presence_col not in data.columns:
                logger.warning(
                    f"Column {presence_col} not found, skipping {motif_type}"
                )
                continue

            # Calculate percentages
            pct_coding = (
                data[data["coding_class"] == "coding"][presence_col].mean() * 100
            )
            pct_lncrna = (
                data[data["coding_class"] == "lncRNA"][presence_col].mean() * 100
            )

            # Create contingency table
            contingency = pd.crosstab(data["coding_class"], data[presence_col])

            # Perform chi-square test
            try:
                chi2, p_value, dof, expected = chi2_contingency(contingency)

                # Calculate Cramér's V (effect size)
                n = contingency.sum().sum()
                cramers_v = np.sqrt(chi2 / n)

                significant = p_value < 0.05

                results.append(
                    {
                        "motif_type": motif_type.upper(),
                        "pct_coding": pct_coding,
                        "pct_lncrna": pct_lncrna,
                        "chi2": chi2,
                        "p_value": p_value,
                        "cramers_v": cramers_v,
                        "significant": significant,
                    }
                )

            except Exception as e:
                logger.warning(f"Could not perform test for {motif_type}: {e}")

        results_df = pd.DataFrame(results)
        logger.info(f"Completed {len(results_df)} motif type tests")

        return results_df

    def save_results(self, overall_stats: Tuple, motif_results: pd.DataFrame):
        """Save analysis results."""

        # Save motif type results
        motif_output = f"{self.output_prefix}_motif_type_chi_square.csv"
        motif_results.to_csv(motif_output, index=False)
        logger.info(f"Saved motif type results to {motif_output}")

        # Save report
        report_output = f"{self.output_prefix}_contingency_report.txt"
        with open(report_output, "w") as f:
            f.write("Non-B DNA Contingency Analysis Report\n")
            f.write("=" * 80 + "\n\n")

            # Overall statistics
            data = self.features[
                self.features["coding_class"].isin(["coding", "lncRNA"])
            ]
            n_coding = len(data[data["coding_class"] == "coding"])
            n_lncrna = len(data[data["coding_class"] == "lncRNA"])

            f.write(f"Total coding transcripts: {n_coding}\n")
            f.write(f"Total lncRNA transcripts: {n_lncrna}\n\n")

            # Overall Non-B presence
            pct_coding_any = (
                data[data["coding_class"] == "coding"]["any_nonb_present"].mean() * 100
            )
            pct_lncrna_any = (
                data[data["coding_class"] == "lncRNA"]["any_nonb_present"].mean() * 100
            )

            f.write(
                f"Percentage of coding transcripts with Non-B DNA: {pct_coding_any:.2f}%\n"
            )
            f.write(
                f"Percentage of lncRNA transcripts with Non-B DNA: {pct_lncrna_any:.2f}%\n\n"
            )

            # Overall test
            chi2, p_value, significant = overall_stats
            f.write(f"Chi-square statistic: {chi2:.4f}\n")
            f.write(f"P-value: {p_value:.4e}\n")
            f.write(f"Significant difference: {significant}\n\n")

            # Motif type summary
            n_significant = motif_results["significant"].sum()
            n_total = len(motif_results)
            f.write(
                f"Motif Types - Significant differences: {n_significant} / {n_total}\n\n"
            )

            # Top motif differences
            f.write("Top 10 motif types by effect size (Cramér's V):\n")
            top_motifs = motif_results.nlargest(10, "cramers_v")[
                ["motif_type", "pct_coding", "pct_lncrna", "cramers_v", "p_value"]
            ]
            f.write(top_motifs.to_string(index=False))
            f.write("\n\n")

            # Detailed results
            f.write("Detailed Motif Type Results:\n")
            f.write("-" * 80 + "\n")
            for _, row in motif_results.iterrows():
                f.write(f"\n{row['motif_type']}:\n")
                f.write(f"  Coding:  {row['pct_coding']:.2f}%\n")
                f.write(f"  lncRNA:  {row['pct_lncrna']:.2f}%\n")
                f.write(f"  Chi2:    {row['chi2']:.4f}\n")
                f.write(f"  P-value: {row['p_value']:.4e}\n")
                f.write(f"  Effect:  {row['cramers_v']:.4f}\n")
                f.write(f"  Sig:     {row['significant']}\n")

        logger.info(f"Saved report to {report_output}")


def main():
    parser = argparse.ArgumentParser(
        description="Perform Non-B DNA contingency analysis"
    )
    parser.add_argument("--features", required=True, help="Feature CSV file")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = NonBContingencyAnalyzer(
        features_file=args.features, output_prefix=args.output_prefix
    )

    # Run analysis
    analyzer.load_features()
    overall_stats = analyzer.perform_overall_test()
    motif_results = analyzer.perform_motif_type_tests()
    analyzer.save_results(overall_stats, motif_results)

    logger.info("Contingency analysis complete!")


if __name__ == "__main__":
    main()

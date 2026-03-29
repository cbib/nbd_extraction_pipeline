#!/usr/bin/env python3
"""
Non-B DNA Statistical Analyzer
===============================
Perform comprehensive statistical analysis of Non-B DNA features
comparing protein-coding vs lncRNA transcripts.

Includes:
- Univariate tests (t-tests, Mann-Whitney U)
- Effect size calculation (Cohen's d)
- Multiple testing correction (Bonferroni, FDR)
- Feature importance analysis
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class NonBStatisticalAnalyzer:
    """Comprehensive statistical analysis of Non-B DNA features."""

    def __init__(self, features_file: str, output_prefix: str):
        """
        Initialize statistical analyzer.

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
        self.numeric_features = None

    def load_features(self):
        """Load and prepare feature matrix."""
        logger.info(f"Loading features from {self.features_file}")
        self.features = pd.read_csv(self.features_file)
        logger.info(f"Loaded {len(self.features)} transcripts")

        # Verify required columns
        if "coding_class" not in self.features.columns:
            raise ValueError("Feature file must contain 'coding_class' column")

        # Identify numeric feature columns
        exclude_cols = [
            "transcript_id",
            "transcript_id_base",
            "transcript_type",
            "coding_class",
            "transcript_length",
        ]
        self.numeric_features = [
            col
            for col in self.features.columns
            if col not in exclude_cols
            and pd.api.types.is_numeric_dtype(self.features[col])
        ]

        logger.info(
            f"Identified {len(self.numeric_features)} numeric features for analysis"
        )

    def perform_univariate_tests(self) -> pd.DataFrame:
        """
        Perform univariate statistical tests for each feature.

        Returns:
        --------
        pd.DataFrame
            Results with columns: feature, mean_coding, mean_lncrna,
            t_statistic, p_value_ttest, u_statistic, p_value_mannwhitney,
            cohens_d, significant
        """
        logger.info("Performing univariate tests...")

        # Filter to coding and lncRNA
        data = self.features[
            self.features["coding_class"].isin(["coding", "lncRNA"])
        ].copy()
        coding_data = data[data["coding_class"] == "coding"]
        lncrna_data = data[data["coding_class"] == "lncRNA"]

        results = []

        for feature in self.numeric_features:
            coding_values = coding_data[feature].dropna()
            lncrna_values = lncrna_data[feature].dropna()

            if len(coding_values) < 2 or len(lncrna_values) < 2:
                logger.warning(f"Insufficient data for {feature}, skipping")
                continue

            # Calculate means
            mean_coding = coding_values.mean()
            mean_lncrna = lncrna_values.mean()

            # T-test (parametric)
            try:
                t_stat, p_ttest = stats.ttest_ind(coding_values, lncrna_values)
            except Exception as e:
                logger.warning(f"T-test failed for {feature}: {e}")
                t_stat, p_ttest = np.nan, np.nan

            # Mann-Whitney U test (non-parametric)
            try:
                u_stat, p_mann = stats.mannwhitneyu(
                    coding_values, lncrna_values, alternative="two-sided"
                )
            except Exception as e:
                logger.warning(f"Mann-Whitney test failed for {feature}: {e}")
                u_stat, p_mann = np.nan, np.nan

            # Cohen's d (effect size)
            pooled_std = np.sqrt(
                (coding_values.std() ** 2 + lncrna_values.std() ** 2) / 2
            )
            cohens_d = (mean_coding - mean_lncrna) / pooled_std if pooled_std > 0 else 0

            results.append(
                {
                    "feature": feature,
                    "mean_coding": mean_coding,
                    "mean_lncrna": mean_lncrna,
                    "std_coding": coding_values.std(),
                    "std_lncrna": lncrna_values.std(),
                    "t_statistic": t_stat,
                    "p_value_ttest": p_ttest,
                    "u_statistic": u_stat,
                    "p_value_mannwhitney": p_mann,
                    "cohens_d": cohens_d,
                }
            )

        results_df = pd.DataFrame(results)

        # Multiple testing correction
        if len(results_df) > 0:
            # Bonferroni
            results_df["p_bonferroni"] = results_df["p_value_ttest"] * len(results_df)
            results_df["p_bonferroni"] = results_df["p_bonferroni"].clip(upper=1.0)

            # FDR (Benjamini-Hochberg)
            results_df = results_df.sort_values("p_value_ttest")
            n = len(results_df)
            results_df["p_fdr"] = results_df["p_value_ttest"] * n / (np.arange(n) + 1)
            results_df["p_fdr"] = results_df["p_fdr"].clip(upper=1.0)

            # Significance flags
            results_df["significant_ttest"] = results_df["p_value_ttest"] < 0.05
            results_df["significant_mann"] = results_df["p_value_mannwhitney"] < 0.05
            results_df["significant_fdr"] = results_df["p_fdr"] < 0.05

        logger.info(f"Completed univariate tests for {len(results_df)} features")

        return results_df

    def perform_feature_importance(self) -> Tuple[pd.DataFrame, Dict[str, float]]:
        """
        Calculate feature importance using Random Forest.

        Returns:
        --------
        Tuple[pd.DataFrame, Dict[str, float]]
            Feature importance scores and performance metrics
        """
        logger.info("Calculating feature importance with Random Forest...")

        # Prepare data
        data = self.features[
            self.features["coding_class"].isin(["coding", "lncRNA"])
        ].copy()
        X = data[self.numeric_features].fillna(0)
        y = (data["coding_class"] == "coding").astype(int)

        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Train/test split for performance evaluation
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, test_size=0.2, random_state=42, stratify=y
        )

        # Train Random Forest
        rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
        rf.fit(X_train, y_train)

        # Evaluate performance
        y_pred = rf.predict(X_test)
        y_proba = rf.predict_proba(X_test)[:, 1]
        performance = {
            "accuracy": accuracy_score(y_test, y_pred),
            "precision": precision_score(y_test, y_pred, zero_division=0),
            "recall": recall_score(y_test, y_pred, zero_division=0),
            "f1": f1_score(y_test, y_pred, zero_division=0),
            "roc_auc": (
                roc_auc_score(y_test, y_proba) if len(np.unique(y_test)) > 1 else np.nan
            ),
            "pr_auc": (
                average_precision_score(y_test, y_proba)
                if len(np.unique(y_test)) > 1
                else np.nan
            ),
        }
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        performance.update({"tn": tn, "fp": fp, "fn": fn, "tp": tp})

        # Extract importances
        importances = pd.DataFrame(
            {"feature": self.numeric_features, "importance": rf.feature_importances_}
        ).sort_values("importance", ascending=False)

        logger.info("Feature importance calculation complete")

        return importances, performance

    def save_results(
        self,
        univariate_results: pd.DataFrame,
        importance_results: pd.DataFrame,
        performance_results: Dict[str, float],
    ):
        """Save all analysis results."""

        # Save univariate results
        univar_output = f"{self.output_prefix}_univariate_tests.csv"
        univariate_results.to_csv(univar_output, index=False)
        logger.info(f"Saved univariate results to {univar_output}")

        # Save feature importance
        import_output = f"{self.output_prefix}_feature_importance.csv"
        importance_results.to_csv(import_output, index=False)
        logger.info(f"Saved feature importance to {import_output}")

        # Save summary report
        report_output = f"{self.output_prefix}_statistical_report.txt"
        with open(report_output, "w") as f:
            f.write("Non-B DNA Statistical Analysis Report\n")
            f.write("=" * 80 + "\n\n")

            # Dataset summary
            data = self.features[
                self.features["coding_class"].isin(["coding", "lncRNA"])
            ]
            f.write(f"Total transcripts: {len(data)}\n")
            f.write(f"  Coding: {len(data[data['coding_class'] == 'coding'])}\n")
            f.write(f"  lncRNA: {len(data[data['coding_class'] == 'lncRNA'])}\n\n")

            # Univariate test summary
            n_sig_ttest = univariate_results["significant_ttest"].sum()
            n_sig_mann = univariate_results["significant_mann"].sum()
            n_sig_fdr = univariate_results["significant_fdr"].sum()
            n_total = len(univariate_results)

            f.write("Univariate Test Summary:\n")
            f.write(f"  Total features tested: {n_total}\n")
            f.write(f"  Significant (t-test, p<0.05): {n_sig_ttest}\n")
            f.write(f"  Significant (Mann-Whitney, p<0.05): {n_sig_mann}\n")
            f.write(f"  Significant (FDR<0.05): {n_sig_fdr}\n\n")

            # Random Forest performance summary
            f.write("Random Forest Performance (80/20 split):\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Accuracy: {performance_results['accuracy']:.4f}\n")
            f.write(f"  Precision: {performance_results['precision']:.4f}\n")
            f.write(f"  Recall: {performance_results['recall']:.4f}\n")
            f.write(f"  F1: {performance_results['f1']:.4f}\n")
            f.write(
                f"  ROC AUC: {performance_results['roc_auc']:.4f}\n"
                if not np.isnan(performance_results["roc_auc"])
                else "  ROC AUC: nan\n"
            )
            f.write(
                f"  PR AUC: {performance_results['pr_auc']:.4f}\n"
                if not np.isnan(performance_results["pr_auc"])
                else "  PR AUC: nan\n"
            )
            f.write(
                "  Confusion Matrix (tn, fp, fn, tp): "
                f"{performance_results['tn']}, {performance_results['fp']}, "
                f"{performance_results['fn']}, {performance_results['tp']}\n\n"
            )

            # Top features by effect size
            f.write("Top 20 features by absolute effect size (Cohen's d):\n")
            f.write("-" * 80 + "\n")
            top_effects = univariate_results.nlargest(20, "cohens_d", keep="all")[
                ["feature", "mean_coding", "mean_lncrna", "cohens_d", "p_value_ttest"]
            ]
            f.write(top_effects.to_string(index=False))
            f.write("\n\n")

            # Top relative features by effect size (length-normalized features)
            f.write(
                "Top 20 relative features (length-normalized; endswith pct or kb) by absolute effect size (Cohen's d):\n"
            )
            f.write("-" * 80 + "\n")
            relative_mask = univariate_results["feature"].str.endswith(("pct", "kb"))
            relative_features = univariate_results.loc[relative_mask].copy()
            if not relative_features.empty:
                relative_features["abs_cohens_d"] = relative_features["cohens_d"].abs()
                top_relative = relative_features.nlargest(20, "abs_cohens_d")[
                    [
                        "feature",
                        "mean_coding",
                        "mean_lncrna",
                        "cohens_d",
                        "p_value_ttest",
                    ]
                ]
                f.write(top_relative.to_string(index=False))
            else:
                f.write("No length-normalized features found (ending with pct or kb).")
            f.write("\n\n")

            # Top features by importance
            f.write("Top 20 features by Random Forest importance:\n")
            f.write("-" * 80 + "\n")
            top_import = importance_results.head(20)
            f.write(top_import.to_string(index=False))
            f.write("\n\n")

            # Most significant features
            f.write("Most significant features (p-value):\n")
            f.write("-" * 80 + "\n")
            most_sig = univariate_results.nsmallest(20, "p_value_ttest")[
                ["feature", "mean_coding", "mean_lncrna", "p_value_ttest", "cohens_d"]
            ]
            f.write(most_sig.to_string(index=False))
            f.write("\n")

        logger.info(f"Saved statistical report to {report_output}")


def main():
    parser = argparse.ArgumentParser(
        description="Perform comprehensive statistical analysis of Non-B DNA features"
    )
    parser.add_argument("--features", required=True, help="Feature CSV file")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix")

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = NonBStatisticalAnalyzer(
        features_file=args.features, output_prefix=args.output_prefix
    )

    # Run analysis
    analyzer.load_features()
    univariate_results = analyzer.perform_univariate_tests()
    importance_results, performance_results = analyzer.perform_feature_importance()
    analyzer.save_results(univariate_results, importance_results, performance_results)

    logger.info("Statistical analysis complete!")


if __name__ == "__main__":
    main()

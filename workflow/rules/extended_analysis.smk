# Non-B DNA Extended Analysis Rules
# ==================================
# Comprehensive feature extraction and statistical analysis
# Organized by execution order

GFA_MOTIFS = ["APR", "DR", "g4Discovery", "IR", "MR", "STR", "TRI", "Z"]

# ============================================================================
# EXECUTION ORDER 1: Preprocessing - Overlap Analysis
# ============================================================================

rule analyze_motif_overlaps:
    """
    Analyze overlaps within each Non-B DNA motif BED file.

    Generates comprehensive statistics to inform overlap resolution strategy:
    - Overlap frequency and distribution
    - Cluster sizes (groups of mutually overlapping regions)
    - Size distributions (overlapping vs isolated regions)
    - Top overlap clusters with genomic locations

    Use these reports to decide whether overlap resolution is needed and
    which strategy to apply (merge, keep longest, sweep-line, etc.).
    """
    input:
        bed = lambda wildcards: f"resources/GRCh38_NonBDNA/GCA_000001405.15_GRCh38_no_alt_analysis_set.{get_motif_filename(wildcards.motif)}_clean.bed"
    output:
        report = "results/overlap_analysis/{motif}_overlap_report.txt"
    conda:
        "lnc-datasets"
    log:
        "logs/overlap_analysis/{motif}.log"
    shell:
        """
        python workflow/scripts/analyze_bed_overlaps.py {input.bed} {output.report} > {log} 2>&1
        """


rule analyze_all_motif_overlaps:
    """
    Run overlap analysis on all Non-B DNA motif types.

    Output: Individual reports for each motif in results/overlap_analysis/
    """
    input:
        expand(
            "results/overlap_analysis/{motif}_overlap_report.txt",
            motif=GFA_MOTIFS
        )
    output:
        touch("results/overlap_analysis/all_motifs_analyzed.flag")


# ============================================================================
# EXECUTION ORDER 3.0: BED-Based Intersections For Extended Analysis
# ============================================================================

rule extended_bedtools_intersect:
    """
    Intersect transcript/exon BED annotations with Non-B DNA motif BED files.

    The left-hand BED is chosen based on the per-dataset ``feature_source`` config
    key ("exons" or "full_transcripts"; default "exons").
    """
    input:
        left = lambda wildcards: (
            f"results/{wildcards.dataset}/exons.bed"
            if get_feature_source(wildcards.dataset) == "exons"
            else f"results/{wildcards.dataset}/transcripts.bed"
        ),
        right = lambda wildcards: get_dataset_motif_bed(
            wildcards.dataset, wildcards.motif_file
        ),
    output:
        overlap = "results/{dataset}/extended_analysis/isect_{motif_file}.bed",
    params:
        extra = "-wa -wb",
    log:
        "logs/{dataset}/extended_analysis/intersect_{motif_file}.log",
    resources:
        mem_mb = 10000,
        runtime = 60,
    wildcard_constraints:
        motif_file = r"g4Discovery|gfa\.(APR|DR|IR|MR|STR|TRI|Z)",
    wrapper:
        "v7.3.0/bio/bedtools/intersect"


# ============================================================================
# EXECUTION ORDER 3.1: Feature Extraction
# ============================================================================

rule extended_feature_extraction:
    """
    Extract comprehensive Non-B DNA features with overlap resolution.
    Processes intersection BED files to calculate 100+ features per transcript.

    Dependencies: create_exons_bed, create_biotypes_from_fasta,
    and extended_bedtools_intersect
    Downstream: extended_contingency_analysis, extended_statistical_analysis

    NOTE: All files accessed by Python script MUST be declared as inputs for
    Snakemake traceability and dependency tracking.
    """
    input:
        transcripts_bed = lambda wildcards: (
            f"results/{wildcards.dataset}/exons.bed"
            if get_feature_source(wildcards.dataset) == "exons"
            else f"results/{wildcards.dataset}/transcripts.bed"
        ),
        biotypes = "results/{dataset}/biotypes.tsv",
        # Explicit BED-based intersection files for the extended pipeline.
        apr = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "APR"),
        dr = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "DR"),
        gq = lambda wildcards: get_extended_isect_bed(
            wildcards.dataset, "g4Discovery"
        ),
        ir = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "IR"),
        mr = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "MR"),
        _str = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "STR"),
        tri = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "TRI"),
        z = lambda wildcards: get_extended_isect_bed(wildcards.dataset, "Z"),
    output:
        features = "results/{dataset}/extended_analysis/features_nonb_features.csv",
        summary = "results/{dataset}/extended_analysis/features_nonb_summary.txt",
        resolved = expand(
            "results/{{dataset}}/extended_analysis/features_resolved_{motif}.bed",
            motif=GFA_MOTIFS
        ),
    params:
        input_dir = "results/{dataset}",
        output_prefix = "results/{dataset}/extended_analysis/features",
        dataset = "{dataset}",
        exon_mode = lambda wildcards: get_feature_source(wildcards.dataset) == "exons",
    log:
        "logs/{dataset}/extended_feature_extraction.log",
    threads: 1
    resources:
        mem_mb = 100000,
        runtime = 3600,
    script:
        "../scripts/extended_feature_extraction_wrapper.py"


# ============================================================================
# EXECUTION ORDER 3.2: Contingency Analysis
# ============================================================================

rule extended_contingency_analysis:
    """
    Perform chi-square contingency tests for Non-B DNA enrichment.
    Compares coding vs lncRNA transcripts for each motif type.

    Dependencies: extended_feature_extraction
    Downstream: create_summary_report
    """
    input:
        features = "results/{dataset}/extended_analysis/features_nonb_features.csv",
    output:
        chi_square = "results/{dataset}/extended_analysis/contingency_motif_type_chi_square.csv",
        report = "results/{dataset}/extended_analysis/contingency_contingency_report.txt",
    params:
        output_prefix = "results/{dataset}/extended_analysis/contingency",
    log:
        "logs/{dataset}/extended_contingency_analysis.log",
    threads: 1
    resources:
        mem_mb = 20000,
        runtime = 10,
    script:
        "../scripts/extended_contingency_wrapper.py"


# ============================================================================
# EXECUTION ORDER 3.3: Statistical Analysis
# ============================================================================

rule extended_statistical_analysis:
    """
    Perform comprehensive statistical analysis of Non-B DNA features.
    Includes univariate tests, effect sizes, and feature importance.

    Dependencies: extended_feature_extraction
    Downstream: create_summary_report
    """
    input:
        features = "results/{dataset}/extended_analysis/features_nonb_features.csv",
    output:
        univariate = "results/{dataset}/extended_analysis/statistics_univariate_tests.csv",
        importance = "results/{dataset}/extended_analysis/statistics_feature_importance.csv",
        report = "results/{dataset}/extended_analysis/statistics_statistical_report.txt",
    params:
        output_prefix = "results/{dataset}/extended_analysis/statistics",
    log:
        "logs/{dataset}/extended_statistical_analysis.log",
    threads: 20  # Use all available cores for Random Forest
    resources:
        mem_mb = 20000,
        runtime = 3600,
    script:
        "../scripts/extended_statistical_wrapper.py"


# ============================================================================
# EXECUTION ORDER 3.4: Extended Analysis Completion
# ============================================================================

rule extended_analysis_all:
    """
    Run complete extended analysis pipeline: feature extraction, contingency tests, and statistics.
    Creates a completion marker file.

    Dependencies: All extended analysis rules
    Downstream: create_summary_report
    """
    input:
        features = "results/{dataset}/extended_analysis/features_nonb_features.csv",
        chi_square = "results/{dataset}/extended_analysis/contingency_motif_type_chi_square.csv",
        univariate = "results/{dataset}/extended_analysis/statistics_univariate_tests.csv",
    output:
        summary = "results/{dataset}/extended_analysis/analysis_complete.txt",
    log:
        "logs/{dataset}/extended_analysis_summary.log",
    shell:
        """
        {{
            echo "Extended Non-B DNA Analysis Complete"
            echo "Dataset: {wildcards.dataset}"
            echo "Date: $(date)"
            echo ""
            echo "Output files:"
            echo "- Features: {input.features}"
            echo "- Contingency: {input.chi_square}"
            echo "- Statistics: {input.univariate}"
        }} > {output.summary} 2>&1
        cp {output.summary} {log}
        """

# Non-B DNA Analysis Pipeline
# ============================
# Comprehensive Snakemake workflow for Non-B DNA element analysis
#
# Includes:
# - Basic per-motif analysis with visualizations (parse_overlaps.py)
# - Extended comprehensive analysis with statistical tests (new)

# Configuration
configfile: "config/config.yaml"
configfile: "config/samples.yaml"

# Global variables
DATASETS = config.get("datasets", ["toy", "gencode.v47"])

GFA_MOTIFS = ["APR", "DR", "g4Discovery_plus", "g4Discovery_minus", "IR", "MR", "STR", "TRI", "Z"]

# Include rule files
include: "workflow/rules/common.smk"
include: "workflow/rules/extended_analysis.smk"
include: "workflow/rules/upstream.smk"

# Check that all required resources are present before running the pipeline
_resource_preflight(DATASETS)


# ============================================================================
# RULE: all (Default Target)
# ============================================================================

rule all:
    """
    Default target: Run complete pipeline (basic + extended analysis).
    """
    input:
        # Basic analysis outputs
        expand(
            "results/{dataset}/transcript_gfa.{motif}_summary.tsv",
            dataset=DATASETS,
            motif=GFA_MOTIFS
        ),
        # Extended analysis outputs
        expand(
            "results/{dataset}/extended_analysis/analysis_complete.txt",
            dataset=DATASETS
        ),
    default_target: True


# ============================================================================
# EXECUTION ORDER 1: Data Preparation
# ============================================================================

rule create_transcripts_bed:
    """
    Convert GTF to transcript-level BED file.
    Extracts transcript entries and creates 6-column BED format.

    Output format (BED6):
    - chrom, start, end, transcript_id, score, strand
    """
    input:
        gtf = lambda wildcards: (
            f"resources/{wildcards.dataset}/{wildcards.dataset}_chr22.gtf"
            if wildcards.dataset == "toy"
            else config["samples"].get(wildcards.dataset, {}).get(
                "gtf", f"resources/{wildcards.dataset}.annotation.gtf"
            )
        ),
    output:
        bed = "results/{dataset}/transcripts.bed",
    log:
        "logs/{dataset}/create_transcripts_bed.log",
    conda:
        "base"
    shell:
        """
        awk -F'\\t' '
        $3 == "transcript" {{
            # Extract transcript_id from attributes
            match($9, /transcript_id "([^"]+)"/, tid)
            # Print BED6 format: chr, start (0-based), end, name, score, strand
            print $1 "\\t" $4-1 "\\t" $5 "\\t" tid[1] "\\t0\\t" $7
        }}' {input.gtf} > {output.bed} 2>{log}

        echo "Created BED file with $(wc -l < {output.bed}) transcripts" >> {log}
        """


rule create_exons_bed:
    """
    Convert GTF to exon-level BED file.
    Keeps one BED6 row per exon with transcript_id in column 4.

    Output format (BED6):
    - chrom, start, end, transcript_id, score, strand
    """
    input:
        gtf = lambda wildcards: (
            f"resources/{wildcards.dataset}/{wildcards.dataset}_chr22.gtf"
            if wildcards.dataset == "toy"
            else config["samples"].get(wildcards.dataset, {}).get(
                "gtf", f"resources/{wildcards.dataset}.annotation.gtf"
            )
        ),
    output:
        bed = "results/{dataset}/exons.bed",
    log:
        "logs/{dataset}/create_exons_bed.log",
    conda:
        "base"
    shell:
        """
        awk -F'\t' '
        $3 == "exon" {{
            match($9, /transcript_id "([^"]+)"/, tid)
            print $1 "\t" $4-1 "\t" $5 "\t" tid[1] "\t0\t" $7
        }}' {input.gtf} > {output.bed} 2>{log}

        echo "Created BED file with $(wc -l < {output.bed}) exons" >> {log}
        """


rule create_biotypes_from_fasta:
    """
    Create biotypes TSV file from protein-coding and lncRNA FASTA files.
    Extracts transcript IDs and assigns transcript_type based on source file.

    Output format:
    - transcript_id_base: Transcript ID without version (e.g., ENST00000456328)
    - transcript_type: 'protein_coding' or 'lncRNA'
    """
    input:
        coding_fasta = lambda wildcards: (
            "results/pc_transcript_ids.txt"
            if wildcards.dataset == "toy"
            else "resources/gencode.v47.pc_transcripts.fa"
        ),
        lncRNA_fasta = lambda wildcards: (
            "results/lncrna_transcript_ids.txt"
            if wildcards.dataset == "toy"
            else "resources/gencode.v47.lncRNA_transcripts.fa"
        ),
    output:
        biotypes = "results/{dataset}/biotypes.tsv",
    log:
        "logs/{dataset}/create_biotypes.log",
    conda:
        "base"
    shell:
        """
        {{
            echo -e "transcript_id_base\\ttranscript_type"

            # Extract protein-coding transcript IDs and label them
            if [[ -f {input.coding_fasta} ]]; then
                if grep -q ">" {input.coding_fasta} 2>/dev/null; then
                    # It's a FASTA file
                    grep ">" {input.coding_fasta} | cut -d'|' -f 1 | sed 's/>//g' | \
                    awk '{{split($1, a, "."); print a[1] "\\tprotein_coding"}}'
                else
                    # It's a plain ID list
                    awk '{{split($1, a, "."); print a[1] "\\tprotein_coding"}}' {input.coding_fasta}
                fi
            fi

            # Extract lncRNA transcript IDs and label them
            if [[ -f {input.lncRNA_fasta} ]]; then
                if grep -q ">" {input.lncRNA_fasta} 2>/dev/null; then
                    # It's a FASTA file
                    grep ">" {input.lncRNA_fasta} | cut -d'|' -f 1 | sed 's/>//g' | \
                    awk '{{split($1, a, "."); print a[1] "\\tlncRNA"}}'
                else
                    # It's a plain ID list
                    awk '{{split($1, a, "."); print a[1] "\\tlncRNA"}}' {input.lncRNA_fasta}
                fi
            fi
        }} > {output.biotypes} 2>{log}

        echo "Created biotypes file with $(tail -n +2 {output.biotypes} | wc -l) transcripts" >> {log}
        """


rule prepare_transcript_ids:
    """
    Extract transcript IDs for protein-coding and lncRNA from biotypes file.
    Required for downstream classification.
    """
    input:
        biotypes = "results/{dataset}/biotypes.tsv",
    output:
        pc_ids = "results/{dataset}/annotation/pc_transcript_ids.txt",
        lnc_ids = "results/{dataset}/annotation/lncrna_transcript_ids.txt",
    log:
        "logs/{dataset}/prepare_transcript_ids.log",
    threads: 1
    resources:
        mem_mb = 500,
    script:
        "workflow/scripts/prepare_transcript_ids.py"


# ============================================================================
# EXECUTION ORDER 2: Basic Analysis (Per-Motif)
# ============================================================================

rule basic_gtf_intersect:
    """
    Intersect annotation GTF with NBD motif BED for the basic analysis path.
    Produces the 18-column (9 GTF + 9 BED) file that parse_overlaps.py expects.
    """
    input:
        left = lambda wildcards: (
            f"resources/{wildcards.dataset}/{wildcards.dataset}_chr22.gtf"
            if wildcards.dataset == "toy"
            else config["samples"].get(wildcards.dataset, {}).get(
                "gtf", f"resources/{wildcards.dataset}.annotation.gtf"
            )
        ),
        right = lambda wildcards: (
            f"resources/toy/toy_{wildcards.motif_file}_chr22_clean.bed"
            if wildcards.dataset == "toy"
            else f"{_NBD}.{wildcards.motif_file}_clean.bed"
        ),
    output:
        overlap = "results/{dataset}/isect_{motif_file}.bed",
    params:
        extra = "-wa -wb",
    log:
        "logs/{dataset}/basic_isect_{motif_file}.log",
    resources:
        mem_mb = 10000,
        runtime = 60,
    wildcard_constraints:
        motif_file = r"g4Discovery(_plus|_minus)?|gfa\.(APR|DR|IR|MR|STR|TRI|Z)",
    wrapper:
        "v7.3.0/bio/bedtools/intersect"


rule basic_motif_analysis:
    """
    Basic per-motif analysis using parse_overlaps.py.
    Generates summary statistics and KDE distribution plots.
    """
    input:
        overlap = lambda wildcards: get_isect_bed(wildcards.dataset, wildcards.motif),
        pc_id_file = "results/{dataset}/annotation/pc_transcript_ids.txt",
        lnc_id_file = "results/{dataset}/annotation/lncrna_transcript_ids.txt",
    output:
        summary = "results/{dataset}/transcript_gfa.{motif}_summary.tsv",
        plot = "results/{dataset}/gfa.{motif}_distributions.png",
    params:
        output_dir = "results/{dataset}",
        nbd_type = "gfa.{motif}",
        keep_other = False,
    log:
        "logs/{dataset}/basic_analysis_gfa.{motif}.log",
    conda:
        "lnc-datasets"
    threads: 1
    resources:
        mem_mb = 10000,
    script:
        "workflow/scripts/parse_overlaps.py"


# ============================================================================
# EXECUTION ORDER 3: Extended Analysis (see rules/extended_analysis.smk)
# ============================================================================
# Rules in extended_analysis.smk:
#   1. extended_feature_extraction
#   2. extended_contingency_analysis
#   3. extended_statistical_analysis
#   4. extended_analysis_all


# ============================================================================
# EXECUTION ORDER 4: Summary Reports
# ============================================================================

rule create_summary_report:
    """
    Create a summary report combining basic and extended analysis results.
    """
    input:
        basic_summaries = expand(
            "results/{{dataset}}/transcript_gfa.{motif}_summary.tsv",
            motif=GFA_MOTIFS
        ),
        extended_summary = "results/{dataset}/extended_analysis/features_nonb_summary.txt",
        contingency = "results/{dataset}/extended_analysis/contingency_contingency_report.txt",
        statistics = "results/{dataset}/extended_analysis/statistics_statistical_report.txt",
    output:
        report = "results/{dataset}/complete_analysis_report.txt",
    log:
        "logs/{dataset}/create_summary_report.log",
    shell:
        """
        {{
            echo "================================================================================"
            echo "Non-B DNA Complete Analysis Report"
            echo "================================================================================"
            echo "Dataset: {wildcards.dataset}"
            echo "Analysis Date: $(date)"
            echo ""
            echo "================================================================================"
            echo "BASIC ANALYSIS SUMMARIES"
            echo "================================================================================"
            for f in {input.basic_summaries}; do
                echo ""
                echo "--- $(basename $f) ---"
                head -20 $f
            done
            echo ""
            echo "================================================================================"
            echo "EXTENDED ANALYSIS - FEATURE SUMMARY"
            echo "================================================================================"
            cat {input.extended_summary}
            echo ""
            echo "================================================================================"
            echo "EXTENDED ANALYSIS - CONTINGENCY TESTS"
            echo "================================================================================"
            cat {input.contingency}
            echo ""
            echo "================================================================================"
            echo "EXTENDED ANALYSIS - STATISTICAL TESTS"
            echo "================================================================================"
            cat {input.statistics}
            echo ""
            echo "================================================================================"
            echo "Analysis complete. See individual files for detailed results."
            echo "================================================================================"
        }} > {output.report} 2>{log}
        """


# ============================================================================
# EXECUTION ORDER 5: Utility Rules
# ============================================================================

rule clean_results:
    """
    Clean all analysis results for a dataset.
    """
    params:
        dataset_dir = "results/{dataset}",
    shell:
        """
        rm -rf {params.dataset_dir}/extended_analysis
        rm -f {params.dataset_dir}/transcript_gfa.*
        rm -f {params.dataset_dir}/gfa.*
        rm -f {params.dataset_dir}/complete_analysis_report.txt
        echo "Cleaned results for dataset: {wildcards.dataset}"
        """


rule list_outputs:
    """
    List all expected output files.
    """
    shell:
        """
        echo ""
        echo "================================================================================"
        echo "Expected Pipeline Outputs"
        echo "================================================================================"
        echo ""
        for dataset in {DATASETS}; do
            echo "Dataset: $dataset"
            echo "----------------------------------------"
            echo ""
            echo "Basic Analysis:"
            for motif in {GFA_MOTIFS}; do
                echo "  - results/$dataset/transcript_gfa.$motif""_summary.tsv"
                echo "  - results/$dataset/gfa.$motif""_distributions.png"
            done
            echo ""
            echo "Extended Analysis:"
            echo "  - results/$dataset/extended_analysis/features_nonb_features.csv"
            echo "  - results/$dataset/extended_analysis/contingency_motif_type_chi_square.csv"
            echo "  - results/$dataset/extended_analysis/statistics_univariate_tests.csv"
            echo ""
            echo "Combined:"
            echo "  - results/$dataset/complete_analysis_report.txt"
            echo ""
        done
        echo "================================================================================"
        """

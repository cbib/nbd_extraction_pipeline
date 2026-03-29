# Non-B DNA Upstream Generation Rules
# ======================================
# Rules to reproduce the Non-B DNA BED files from the GRCh38 assembly FASTA.
#
# Two tools are used:
#   gfa         — identifies APR, DR, IR, MR, STR, Z-DNA, and TRI motifs
#                 https://github.com/abcsFrederick/non-B_gfa
#   g4Discovery — identifies G-quadruplex (G4) motifs
#                 https://github.com/saswat-km/g4Discovery.PanSN
#
# Execution order:
#   1. decompress_assembly        — decompress GRCh38.p14.genome.fa.gz
#   2. install_gfa / install_g4discovery — clone and build tools
#   3. run_gfa                    — run gfa on whole assembly (outputs TSV per motif)
#   4. gfa_tsv_to_bed             — convert each TSV to sorted BED3
#   5. gfa_extract_triplex        — extract TRI from MR TSV
#   6. clean_chrom (common.smk)   — strip PanSN prefixes → _clean.bed
#   7. run_g4discovery            — run g4Discovery on whole assembly
#   8. g4discovery_clean_chrom    — decompress + strip PanSN prefixes
#
# Convenience target: all_upstream

# GFA motif types that gfa outputs directly (excluding GQ and TRI)
_GFA_BASE_MOTIFS = ["APR", "DR", "IR", "MR", "STR", "Z"]
# All GFA-derived downstream types (TRI extracted separately from MR TSV)
_GFA_ALL_MOTIFS = _GFA_BASE_MOTIFS + ["TRI"]

# Assembly prefix (matches the naming of pre-computed BED files)
_ASSEMBLY_PREFIX = "GCA_000001405.15_GRCh38_no_alt_analysis_set"
_NBD_DIR = f"resources/GRCh38_NonBDNA"


# ── Upstream convenience target ──────────────────────────────────────────────

rule all_upstream:
    """Produce all cleaned Non-B DNA BED files (GFA × 7 + g4Discovery × 1)."""
    input:
        expand(
            f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.gfa.{{motif}}_clean.bed",
            motif=_GFA_ALL_MOTIFS,
        ),
        f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.g4Discovery_clean.bed",


# ── Genome assembly ───────────────────────────────────────────────────────────

rule decompress_assembly:
    """Decompress the GRCh38.p14 genome FASTA for tool input."""
    input:
        gz=lambda wc: config.get(
            "assembly_gz",
            "resources/GRCh38.p14.genome.fa.gz",
        ),
    output:
        fa="resources/GRCh38.p14.genome.fa",
    log:
        "logs/upstream/decompress_assembly.log",
    conda:
        "base"
    shell:
        "gunzip -c {input.gz} > {output.fa} 2> {log}"


# ── GFA installation ──────────────────────────────────────────────────────────

rule install_gfa:
    """Clone and compile the non-B_gfa C tool.

    Produces the 'gfa' binary in software/non-B_gfa/gfa.
    Only needs to run once; subsequent runs are skipped because the output
    already exists (Snakemake cache).
    """
    output:
        binary="software/non-B_gfa/gfa",
    log:
        "logs/upstream/install_gfa.log",
    conda:
        "../envs/gfa.yaml"
    shell:
        """
        mkdir -p software
        if [ ! -d software/non-B_gfa/.git ]; then
            git clone https://github.com/abcsFrederick/non-B_gfa.git software/non-B_gfa \
                2> {log}
        fi
        cd software/non-B_gfa
        make 2>> ../../{log}
        echo "gfa binary built at $(pwd)/gfa" >> ../../{log}
        """


# ── g4Discovery installation ──────────────────────────────────────────────────

rule install_g4discovery:
    """Clone the g4Discovery.PanSN Python/R tool.

    No compilation needed; produces the main Python script.
    """
    output:
        script="software/g4Discovery.PanSN/g4Discovery.py",
    log:
        "logs/upstream/install_g4discovery.log",
    conda:
        "../envs/g4discovery.yaml"
    shell:
        """
        mkdir -p software
        if [ ! -d software/g4Discovery.PanSN/.git ]; then
            git clone https://github.com/saswat-km/g4Discovery.PanSN.git \
                software/g4Discovery.PanSN 2> {log}
        fi
        echo "g4Discovery cloned at software/g4Discovery.PanSN" >> {log}
        """


# ── GFA motif annotation ──────────────────────────────────────────────────────

rule run_gfa:
    """Run gfa on the whole GRCh38 assembly to annotate non-B DNA motifs.

    GQ motifs are skipped (-skipGQ) because G-quadruplexes are annotated
    separately with g4Discovery.

    Outputs one TSV file per motif type: APR, DR, IR, MR, STR, Z.
    (TRI is extracted from MR in the gfa_extract_triplex rule.)
    """
    input:
        fa="resources/GRCh38.p14.genome.fa",
        gfa_bin="software/non-B_gfa/gfa",
    output:
        # gfa writes <prefix>_APR.tsv, <prefix>_DR.tsv, etc.
        tsv=expand(
            f"results/gfa_raw/{_ASSEMBLY_PREFIX}_{{motif}}.tsv",
            motif=_GFA_BASE_MOTIFS,
        ),
    params:
        prefix=f"results/gfa_raw/{_ASSEMBLY_PREFIX}",
    log:
        "logs/upstream/run_gfa.log",
    conda:
        "../envs/gfa.yaml"
    threads: 1
    resources:
        mem_mb=65536,
        runtime=1440,  # 24 h — whole-genome annotation is slow
    shell:
        """
        mkdir -p results/gfa_raw
        ./{input.gfa_bin} \
            -seq {input.fa} \
            -out {params.prefix} \
            -skipGQ \
            > {log} 2>&1
        """


rule gfa_tsv_to_bed:
    """Convert a GFA TSV output file to a sorted BED3 file.

    Column mapping (1-based TSV → 0-based BED):
      $1  = chromosome / sequence name
      $4  = start (1-based)  → BED start = $4 - 1
      $5  = stop  (0-based end included)
    """
    input:
        tsv=f"results/gfa_raw/{_ASSEMBLY_PREFIX}_{{motif}}.tsv",
    output:
        bed=f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.gfa.{{motif}}.bed",
    log:
        "logs/upstream/gfa_tsv_to_bed_{motif}.log",
    conda:
        "base"
    wildcard_constraints:
        # Only base motifs — TRI has its own rule
        motif="APR|DR|IR|MR|STR|Z",
    shell:
        r"""
        awk -v OFS='\t' '(NR>1){{s=$4-1; print $1,s,$5}}' {input.tsv} \
            | sort -k1,1 -k2,2n \
            > {output.bed} 2> {log}
        echo "$(wc -l < {output.bed}) records written to {output.bed}" >> {log}
        """


rule gfa_extract_triplex:
    """Extract triplex-forming (TRI) regions from the MR TSV.

    TRI motifs are mirror repeats where column 12 == 1, indicating
    likelihood of triplex DNA formation.
    """
    input:
        mr_tsv=f"results/gfa_raw/{_ASSEMBLY_PREFIX}_MR.tsv",
    output:
        bed=f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.gfa.TRI.bed",
    log:
        "logs/upstream/gfa_extract_triplex.log",
    conda:
        "base"
    shell:
        r"""
        awk -v OFS='\t' '($12==1){{s=$4-1; print $1,s,$5}}' {input.mr_tsv} \
            | sort -k1,1 -k2,2n \
            > {output.bed} 2> {log}
        echo "$(wc -l < {output.bed}) TRI records written to {output.bed}" >> {log}
        """


# ── g4Discovery G-quadruplex annotation ──────────────────────────────────────

rule run_g4discovery:
    """Run g4Discovery on the GRCh38 assembly with default settings.

    Outputs a gzipped BED file containing G4 positions on both strands,
    with both pqsfinder and G4Hunter scores.
    Default thresholds: pqsfinder >= 40, |G4Hunter| >= 1.5, tetrads >= 3.
    """
    input:
        fa="resources/GRCh38.p14.genome.fa",
        script="software/g4Discovery.PanSN/g4Discovery.py",
    output:
        gz=f"results/g4discovery_raw/{_ASSEMBLY_PREFIX}.g4Discovery.bed.gz",
    log:
        "logs/upstream/run_g4discovery.log",
    conda:
        "../envs/g4discovery.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime=2880,  # 48 h — whole-genome G4 prediction is very slow
    shell:
        """
        mkdir -p results/g4discovery_raw
        python3 {input.script} \
            -fa {input.fa} \
            -o {output.gz} \
            > {log} 2>&1
        """


rule g4discovery_clean_chrom:
    """Decompress g4Discovery output and strip PanSN prefixes (GRCh38#0#chr1 → chr1).

    g4Discovery uses PanSN-format chromosome names; the downstream
    bedtools intersect rules expect standard UCSC names (chr1, chr2, …).
    """
    input:
        gz=f"results/g4discovery_raw/{_ASSEMBLY_PREFIX}.g4Discovery.bed.gz",
    output:
        bed=f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.g4Discovery.bed",
        clean=f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.g4Discovery_clean.bed",
    log:
        "logs/upstream/g4discovery_clean_chrom.log",
    conda:
        "base"
    shell:
        """
        # Decompress
        gunzip -c {input.gz} > {output.bed} 2> {log}
        # Strip PanSN prefix and sort
        sed 's/GRCh38#0#//g' {output.bed} \
            | sort -k1,1 -k2,2n \
            > {output.clean} 2>> {log}
        echo "$(wc -l < {output.clean}) G4 records in clean BED" >> {log}
        """

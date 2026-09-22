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

_ASSEMBLY_PREFIX = config["upstream"]["assembly_prefix"]
_NBD_DIR = config["upstream"]["nbd_dir"]
_TOOLS_DIR = config["upstream"]["tools_dir"]
_G4_CONTIGS = config["upstream"]["g4_contigs"]


# ── Upstream convenience target ──────────────────────────────────────────────


rule all_upstream:
    """Produce all cleaned Non-B DNA BED files (GFA × 7 + g4Discovery × 3)."""
    input:
        expand(
            f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.gfa.{{motif}}_clean.bed",
            motif=_GFA_ALL_MOTIFS,
        ),
        f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.g4Discovery_clean.bed",
        f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.g4Discovery_plus_clean.bed",
        f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}.g4Discovery_minus_clean.bed",


# ── Genome assembly ───────────────────────────────────────────────────────────


rule decompress_assembly:
    """Decompress the GRCh38.p14 genome FASTA for tool input."""
    input:
        gz=config["upstream"]["assembly_gz"],
    output:
        fa=config["upstream"]["assembly_fa"],
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
        binary=f"{_TOOLS_DIR}/non-B_gfa/gfa",
    log:
        "logs/upstream/install_gfa.log",
    conda:
        "../envs/gfa.yaml"
    shell:
        """
        mkdir -p {_TOOLS_DIR}
        if [ ! -d {_TOOLS_DIR}/non-B_gfa/.git ]; then
            git clone {config[upstream][gfa_repository]} {_TOOLS_DIR}/non-B_gfa \
                2> {log}
        fi
        cd {_TOOLS_DIR}/non-B_gfa
        git checkout --detach {config[upstream][gfa_revision]} >> ../../{log} 2>&1
        make 2>> ../../{log}
        echo "gfa binary built at $(pwd)/gfa" >> ../../{log}
        """


# ── g4Discovery installation ──────────────────────────────────────────────────


rule install_g4discovery:
    """Clone the g4Discovery.PanSN Python/R tool.

    No compilation needed; produces the main Python script.
    """
    output:
        script=f"{_TOOLS_DIR}/g4Discovery.PanSN/src/g4Discovery.py",
    log:
        "logs/upstream/install_g4discovery.log",
    conda:
        "../envs/g4discovery.yaml"
    shell:
        """
        mkdir -p {_TOOLS_DIR}
        if [ ! -d {_TOOLS_DIR}/g4Discovery.PanSN/.git ]; then
            git clone {config[upstream][g4discovery_repository]} \
                {_TOOLS_DIR}/g4Discovery.PanSN 2> {log}
        fi
        cd {_TOOLS_DIR}/g4Discovery.PanSN
        git checkout --detach {config[upstream][g4discovery_revision]} >> ../../{log} 2>&1
        test -f src/g4Discovery.py
        echo "g4Discovery cloned at $(pwd)" >> ../../{log}
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
        fa=config["upstream"]["assembly_fa"],
        gfa_bin=f"{_TOOLS_DIR}/non-B_gfa/gfa",
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


rule split_assembly_for_g4discovery:
    """Extract configured GRCh38 contigs as one-record FASTA files."""
    input:
        fa=config["upstream"]["assembly_fa"],
    output:
        expand("results/g4discovery_raw/fasta/{contig}.fa", contig=_G4_CONTIGS),
    params:
        contigs=" ".join(_G4_CONTIGS),
    log:
        "logs/upstream/split_assembly_for_g4discovery.log",
    conda:
        "base"
    shell:
        r"""
        mkdir -p results/g4discovery_raw/fasta
        awk -v outdir="results/g4discovery_raw/fasta" -v contigs="{params.contigs}" '
            BEGIN {{ count = split(contigs, items, " "); for (i = 1; i <= count; i++) wanted[items[i]] = 1 }}
            /^>/ {{ contig = substr($1, 2); write_record = (contig in wanted) }}
            write_record {{ print > (outdir "/" contig ".fa") }}
        ' {input.fa} > {log}
        for contig in {params.contigs}; do test -s "results/g4discovery_raw/fasta/$contig.fa"; done
        """


rule run_g4discovery_contig:
    """Run default G4Discovery on one GRCh38 contig, retaining both scores."""
    input:
        fa="results/g4discovery_raw/fasta/{contig}.fa",
        script=f"{_TOOLS_DIR}/g4Discovery.PanSN/src/g4Discovery.py",
    output:
        gz="results/g4discovery_raw/contigs/{contig}.g4Discovery.bed.gz",
    log:
        "logs/upstream/run_g4discovery_{contig}.log",
    conda:
        "../envs/g4discovery.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime=2880,
    shell:
        """
        mkdir -p results/g4discovery_raw/contigs
        python3 {input.script} -fa {input.fa} -chr {wildcards.contig} -o {output.gz} > {log} 2>&1
        """


rule merge_g4discovery:
    """Merge per-contig default G4Discovery BEDs without dropping score columns."""
    input:
        expand("results/g4discovery_raw/contigs/{contig}.g4Discovery.bed.gz", contig=_G4_CONTIGS),
    output:
        gz=f"results/g4discovery_raw/{_ASSEMBLY_PREFIX}.g4Discovery.bed.gz",
    log:
        "logs/upstream/merge_g4discovery.log",
    conda:
        "base"
    shell:
        """
        zcat {input} | sort -k1,1 -k2,2n | gzip -c > {output.gz}
        echo "$(zcat {output.gz} | wc -l) G4 records merged" > {log}
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


# ponytail: g4discovery_split_strands lives in common.smk so it's available without upstream tools

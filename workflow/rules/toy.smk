# ============================================================================
# Reproducible toy resource rules
# ============================================================================

_TOY = config["toy"]
_TOY_MOTIF_FILES = [
    "resources/toy/toy_gfa.APR_chr22_clean.bed",
    "resources/toy/toy_gfa.DR_chr22_clean.bed",
    "resources/toy/toy_gfa.IR_chr22_clean.bed",
    "resources/toy/toy_gfa.MR_chr22_clean.bed",
    "resources/toy/toy_gfa.STR_chr22_clean.bed",
    "resources/toy/toy_gfa.TRI_chr22_clean.bed",
    "resources/toy/toy_gfa.Z_chr22_clean.bed",
    "resources/toy/toy_g4Discovery_plus_chr22_clean.bed",
    "resources/toy/toy_g4Discovery_minus_chr22_clean.bed",
]


rule download_toy_fasta:
    """Download and verify the hg38 chromosome 22 sequence."""
    output:
        archive=_TOY["fasta_archive"],
    params:
        url=_TOY["fasta_url"],
        md5=_TOY["fasta_md5"],
    log:
        "logs/toy/download_fasta.log",
    conda:
        "base"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.archive})
        tmp="{output.archive}.tmp.$$"
        trap 'rm -f "$tmp"' EXIT
        curl --fail --location --retry 3 --output "$tmp" "{params.url}" > {log} 2>&1
        printf '%s  %s\n' "{params.md5}" "$tmp" | md5sum --check --status - >> {log} 2>&1
        mv "$tmp" {output.archive}
        trap - EXIT
        """


rule create_toy_fasta:
    input:
        archive=_TOY["fasta_archive"],
    output:
        fa=_TOY["fasta"],
    params:
        md5=_TOY["fasta_md5"],
    conda:
        "base"
    shell:
        r"""
        printf '%s  %s\n' "{params.md5}" {input.archive} | md5sum --check --status -
        gzip -dc {input.archive} > {output.fa}
        test "$(head -n 1 {output.fa})" = '>chr22'
        """


rule run_toy_gfa:
    input:
        fa=_TOY["fasta"],
        binary=f"{_TOOLS_DIR}/non-B_gfa/gfa",
    output:
        expand("results/toy/gfa_raw/chr22_{motif}.tsv", motif=_GFA_BASE_MOTIFS),
    log:
        "logs/toy/run_gfa.log",
    conda:
        "../envs/gfa.yaml"
    threads: 1
    resources:
        mem_mb=65536,
    shell:
        """
        mkdir -p results/toy/gfa_raw $(dirname {log})
        ./{input.binary} -seq {input.fa} -out results/toy/gfa_raw/chr22 -skipGQ > {log} 2>&1
        """


rule toy_gfa_bed:
    input:
        "results/toy/gfa_raw/chr22_{motif}.tsv",
    output:
        "resources/toy/toy_gfa.{motif}_chr22_clean.bed",
    conda:
        "base"
    wildcard_constraints:
        motif="APR|DR|IR|MR|STR|Z",
    shell:
        r"""awk -v OFS='\t' 'NR>1 {{print $1,$4-1,$5}}' {input} | sort -k1,1 -k2,2n > {output}"""


rule toy_triplex_bed:
    input:
        "results/toy/gfa_raw/chr22_MR.tsv",
    output:
        "resources/toy/toy_gfa.TRI_chr22_clean.bed",
    conda:
        "base"
    shell:
        r"""awk -v OFS='\t' '$12==1 {{print $1,$4-1,$5}}' {input} | sort -k1,1 -k2,2n > {output}"""


rule run_toy_g4discovery:
    input:
        fa=_TOY["fasta"],
        script=f"{_TOOLS_DIR}/g4Discovery.PanSN/src/g4Discovery.py",
    output:
        "results/toy/chr22.g4Discovery.bed.gz",
    log:
        "logs/toy/run_g4discovery.log",
    conda:
        "../envs/g4discovery.yaml"
    threads: 1
    resources:
        mem_mb=32768,
    shell:
        """
        mkdir -p results/toy $(dirname {log})
        python3 {input.script} -fa {input.fa} -chr chr22 -o {output} > {log} 2>&1
        """


rule toy_g4discovery_bed:
    input:
        "results/toy/chr22.g4Discovery.bed.gz",
    output:
        "resources/toy/toy_g4Discovery_chr22_clean.bed",
    conda:
        "base"
    shell:
        """gzip -dc {input} | sed 's/GRCh38#0#//g' | sort -k1,1 -k2,2n > {output}"""


rule download_toy_annotation:
    """Download the pinned Gencode annotation and verify its MD5 before publish."""
    output:
        archive=_TOY["archive"],
    params:
        url=_TOY["annotation_url"],
        md5=_TOY["annotation_md5"],
    log:
        "logs/toy/download_annotation.log",
    conda:
        "base"
    shell:
        r"""
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.archive})
        tmp="{output.archive}.tmp.$$"
        trap 'rm -f "$tmp"' EXIT
        curl --fail --location --retry 3 --output "$tmp" "{params.url}" \
            > {log} 2>&1
        printf '%s  %s\n' "{params.md5}" "$tmp" | md5sum --check --status - \
            >> {log} 2>&1
        mv "$tmp" {output.archive}
        trap - EXIT
        """


rule verify_toy_annotation:
    """Recheck the downloaded archive before it enters the extraction chain."""
    input:
        archive=_TOY["archive"],
    output:
        verified=_TOY["archive"] + ".md5.ok",
    params:
        md5=_TOY["annotation_md5"],
    log:
        "logs/toy/verify_annotation.log",
    conda:
        "base"
    shell:
        r"""
        mkdir -p $(dirname {log})
        printf '%s  %s\n' "{params.md5}" {input.archive} | md5sum --check - \
            > {log} 2>&1
        touch {output.verified}
        """


rule create_toy_annotation:
    """Extract all chr22 GTF records while preserving source headers."""
    input:
        archive=_TOY["archive"],
        verified=_TOY["archive"] + ".md5.ok",
    output:
        gtf=_TOY["extracted_gtf"],
    params:
        chromosome=_TOY["chromosome"],
    log:
        "logs/toy/create_annotation.log",
    conda:
        "base"
    shell:
        r"""
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.gtf})
        gzip -dc {input.archive} \
            | awk -v chromosome="{params.chromosome}" \
                '/^#/ {{ print; next }} $1 == chromosome {{ print }}' \
            > {output.gtf} 2> {log}
        test -s {output.gtf}
        """


rule validate_toy_resources:
    """Record checksums and record counts for generated and fixture inputs."""
    input:
        archive=_TOY["archive"],
        gtf=_TOY["extracted_gtf"],
        fasta_archive=_TOY["fasta_archive"],
        fasta=_TOY["fasta"],
        motifs=_TOY_MOTIF_FILES,
    output:
        provenance=_TOY["provenance"],
    params:
        release=_TOY["release"],
        assembly=_TOY["assembly"],
        chromosome=_TOY["chromosome"],
        annotation_url=_TOY["annotation_url"],
        checksum_url=_TOY["checksum_url"],
        annotation_md5=_TOY["annotation_md5"],
        fasta_url=_TOY["fasta_url"],
        fasta_md5=_TOY["fasta_md5"],
        gfa_revision=config["upstream"]["gfa_revision"],
        g4discovery_revision=config["upstream"]["g4discovery_revision"],
    log:
        "logs/toy/validate_resources.log",
    conda:
        "base"
    shell:
        r"""
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.provenance})
        {{
            printf 'release: "%s"\n' "{params.release}"
            printf 'assembly: "%s"\n' "{params.assembly}"
            printf 'chromosome: "%s"\n' "{params.chromosome}"
            printf 'annotation_url: "%s"\n' "{params.annotation_url}"
            printf 'checksum_url: "%s"\n' "{params.checksum_url}"
            printf 'annotation_md5: "%s"\n' "{params.annotation_md5}"
            printf 'fasta_url: "%s"\n' "{params.fasta_url}"
            printf 'fasta_md5: "%s"\n' "{params.fasta_md5}"
            printf 'gfa_revision: "%s"\n' "{params.gfa_revision}"
            printf 'g4discovery_revision: "%s"\n' "{params.g4discovery_revision}"
            printf 'fasta_sha256: "%s"\n' "$(sha256sum {input.fasta} | awk '{{print $1}}')"
            printf 'annotation_archive_sha256: "%s"\n' \
                "$(sha256sum {input.archive} | awk '{{print $1}}')"
            printf 'extracted_gtf_sha256: "%s"\n' \
                "$(sha256sum {input.gtf} | awk '{{print $1}}')"
            printf 'extracted_gtf_records: %s\n' \
                "$(grep -vc '^#' {input.gtf})"
            printf 'motif_files:\n'
            for motif in {input.motifs}; do
                printf '  - path: "%s"\n' "$motif"
                printf '    sha256: "%s"\n' "$(sha256sum "$motif" | awk '{{print $1}}')"
                printf '    records: %s\n' "$(wc -l < "$motif")"
            done
        }} > {output.provenance} 2> {log}
        """

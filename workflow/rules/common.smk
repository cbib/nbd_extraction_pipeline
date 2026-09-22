# ============================================================================
# Helper Functions
# ============================================================================

from pathlib import Path
from snakemake.exceptions import WorkflowError


def get_feature_source(sample: str) -> str:
    """Return the configured feature source for a dataset."""
    value = config.get("samples", {}).get(sample, {}).get("feature_source", "exons")
    return "full_transcripts" if value == "transcripts" else value


def _resource_preflight(datasets) -> None:
    """Fail before DAG expansion when configured external inputs are absent."""
    missing = []
    buildable = []
    samples = config.get("samples", {})

    def require_file(path: str, purpose: str) -> None:
        if not Path(path).is_file():
            missing.append(f"  - {path} ({purpose})")

    for dataset in datasets:
        sample = samples.get(dataset)
        if sample is None:
            missing.append(f"  - config/samples.yaml entry for dataset '{dataset}'")
            continue

        gtf = sample.get("gtf")
        if not gtf:
            missing.append(f"  - samples.{dataset}.gtf in config/samples.yaml")
        elif dataset == "toy" and not Path(gtf).is_file():
            buildable.append(f"  - {gtf} (built by rule create_toy_annotation)")
        else:
            require_file(gtf, f"annotation GTF for {dataset}")

        if dataset == "toy":
            # Toy transcript biotypes are derived directly from the checked-in GTF.
            pass
        else:
            require_file(
                "resources/gencode.v47.pc_transcripts.fa",
                f"protein-coding transcript FASTA for {dataset}",
            )
            require_file(
                "resources/gencode.v47.lncRNA_transcripts.fa",
                f"lncRNA transcript FASTA for {dataset}",
            )

    if any(dataset != "toy" for dataset in datasets):
        upstream = config.get("upstream", {})
        require_file(upstream.get("assembly_gz", ""), "configured assembly archive")
        for key in ("gfa_repository", "gfa_revision", "g4discovery_repository", "g4discovery_revision"):
            if not upstream.get(key):
                missing.append(f"  - config/config.yaml upstream.{key}")

        gfa_binary = Path(upstream.get("tools_dir", "software")) / "non-B_gfa/gfa"
        g4_script = Path(upstream.get("tools_dir", "software")) / "g4Discovery.PanSN/src/g4Discovery.py"
        if not gfa_binary.is_file():
            buildable.append(f"  - {gfa_binary} (built by rule install_gfa)")
        if not g4_script.is_file():
            buildable.append(f"  - {g4_script} (fetched by rule install_g4discovery)")

    if missing:
        message = [
            "Resource preflight failed before DAG expansion.",
            "Missing configured external prerequisites:",
            *missing,
        ]
        if buildable:
            message.extend(["", "Missing upstream outputs that Snakemake can build:", *buildable])
        message.extend([
            "",
            "Provide the listed resources or update the declared resource manifest;",
            "do not silently change the configured dataset.",
        ])
        raise WorkflowError("\n".join(message))

def get_motif_filename(motif):
    """
    Map motif names to their corresponding filename patterns.

    Args:
        motif: Motif name (e.g., 'APR', 'g4Discovery')

    Returns:
        Filename pattern for the motif:
        - 'g4Discovery' for g4Discovery (GQ)
        - 'gfa.{motif}' for all other motifs
    """
    motif = str(motif)
    if motif in ("g4Discovery", "g4Discovery_plus", "g4Discovery_minus"):
        return motif
    if motif.startswith("gfa."):
        return motif
    else:
        return f"gfa.{motif}"


def get_dataset_motif_bed(dataset, motif):
    """
    Get the motif BED file path for a dataset.

    Toy datasets use chromosome-restricted BED files, while full datasets use
    the genome-wide precomputed resources.
    """
    filename = get_motif_filename(motif)
    if dataset == "toy":
        return f"resources/toy/toy_{filename}_chr22_clean.bed"
    return f"{_NBD}.{filename}_clean.bed"


def get_isect_bed(dataset, motif):
    """
    Get the intersection BED file path for a given dataset and motif.

    Args:
        dataset: Dataset name (e.g., 'toy', 'gencode.v47')
        motif: Motif name (e.g., 'APR', 'g4Discovery')

    Returns:
        Path to intersection BED file
    """
    filename = get_motif_filename(motif)
    return f"results/{dataset}/isect_{filename}.bed"


def get_extended_isect_bed(dataset, motif):
    """
    Get the extended-analysis BED-based intersection file path.
    """
    filename = get_motif_filename(motif)
    return f"results/{dataset}/extended_analysis/isect_{filename}.bed"


# ============================================================================
# Common Rules
# ============================================================================

_ASSEMBLY_PREFIX = config["upstream"]["assembly_prefix"]
_NBD_DIR = config["upstream"]["nbd_dir"]
_NBD = f"{_NBD_DIR}/{_ASSEMBLY_PREFIX}"


rule g4discovery_split_strands:
    """Split the cleaned G4 BED into plus- and minus-strand files (column 6)."""
    input:
        f"{_NBD}.g4Discovery_clean.bed",
    output:
        plus=f"{_NBD}.g4Discovery_plus_clean.bed",
        minus=f"{_NBD}.g4Discovery_minus_clean.bed",
    log:
        "logs/upstream/g4discovery_split_strands.log",
    shell:
        """
        mkdir -p $(dirname {log})
        awk '$6=="+"' {input} > {output.plus}
        awk '$6=="-"' {input} > {output.minus}
        echo "plus: $(wc -l < {output.plus}), minus: $(wc -l < {output.minus})" > {log}
        """

rule split_toy_g4discovery_strands:
    """Split the checked-in toy G4 BED into strand-specific fixtures."""
    input:
        "resources/toy/toy_g4Discovery_chr22_clean.bed",
    output:
        plus="resources/toy/toy_g4Discovery_plus_chr22_clean.bed",
        minus="resources/toy/toy_g4Discovery_minus_chr22_clean.bed",
    log:
        "logs/toy/split_g4discovery_strands.log",
    shell:
        """
        mkdir -p $(dirname {log})
        awk '$6=="+"' {input} > {output.plus}
        awk '$6=="-"' {input} > {output.minus}
        echo "plus: $(wc -l < {output.plus}), minus: $(wc -l < {output.minus})" > {log}
        """


rule gunzip:
    input:
        "{file}.gz",
    output:
        "{file}",
    log:
        "logs/gunzip/{file}.log",
    conda:
        "base"
    wildcard_constraints:
        file=r".+(gtf|bed)",
    shell:
        "gunzip -c {input} > {output}"

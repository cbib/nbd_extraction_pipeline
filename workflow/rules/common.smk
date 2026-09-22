# ============================================================================
# Helper Functions
# ============================================================================

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

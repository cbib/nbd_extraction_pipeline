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
    if motif == "g4Discovery":
        return "g4Discovery"
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
    return (
        "resources/GRCh38_NonBDNA/"
        f"GCA_000001405.15_GRCh38_no_alt_analysis_set.{filename}_clean.bed"
    )


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

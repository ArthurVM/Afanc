from os import path

from .compute_mash_taxonomic_matrices import make_mash_taxonomic_matrices


def makeVariantIndex(args, base_nodes=None):
    """ Generate the Mash Variant Index (MVI).

    This replaces the older parent-child-only index construction. The screen
    read deconvolution now expects the variant-index directory to contain:

        clean FASTAs + ncbi taxonomy
              |
              V
        all-vs-all Mash distances
              |
              +--> variant_index.json
              +--> mash_all_vs_all.tsv
              +--> mash_genus_matrix.tsv
              +--> mash_species_matrix.tsv

    """

    ## compatibility argument
    _ = base_nodes

    names_dmp = path.join(args.autoDB_WDir, "ncbi_taxonomy", "names.dmp")
    nodes_dmp = path.join(args.autoDB_WDir, "ncbi_taxonomy", "nodes.dmp")

    return make_mash_taxonomic_matrices(
        fasta_dir=args.cleanFasta_WDir,
        names_dmp=names_dmp,
        nodes_dmp=nodes_dmp,
        mash_bin="mash",
        out_dir=args.variant_index_WDir,
        threads=args.threads,
        kmer=16,
        sketch_size=10000,
        prefix="all_genomes",
    )

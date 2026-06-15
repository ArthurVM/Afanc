import math


ELASTIC_DISTANCE_SCALE = 0.0015


def score_lca_nodes(family, scoring_ranks, species_rank, global_threshold, min_reads, total_reads, variant_index):
    """Score nodes inside a family subtree for LCA redistribution."""
    scoring = []
    scoring_species = set()

    ## first pass: global calls at genus/species level
    for node in family.traverse():
        if node.taxon_level not in scoring_ranks:
            continue
        if node.raw_clade_reads <= 0 or node.raw_clade_reads < min_reads:
            continue
        if total_reads and (node.raw_clade_reads / total_reads) >= global_threshold:
            node.scoring_rule = "global_threshold"
            scoring.append(node)
            if node.taxon_level == species_rank:
                scoring_species.add(node)

    ## second pass: elastic calls below a globally-scoring species
    for node in family.traverse():
        if not rank_is_below(node.taxon_level, species_rank):
            continue
        if node.raw_clade_reads <= 0 or node.raw_clade_reads < min_reads:
            continue

        species = node.ancestor_at_taxon_level(species_rank)
        if species not in scoring_species:
            continue

        apply_elastic_threshold(node, variant_index)
        if node.raw_clade_reads >= node.weighted_threshold:
            node.scoring_rule = "elastic_threshold"
            scoring.append(node)

    return scoring


def apply_elastic_threshold(node, variant_index):
    """Calculate and store the parent-child elastic read threshold for a node."""
    if node.parent is None:
        w = 0.0
        raw_w = w
        weighted_threshold = 0.0
        node.threshold_func = "root"
        node.variables = None

    elif str(node.ncbi_taxID) in variant_index:
        parent_similarity = variant_index[str(node.ncbi_taxID)]["parent_index"]["mean"]
        child_similarity = variant_index[str(node.ncbi_taxID)]["sibling_index"]["mean"]
        parent_child_distance = max(0.0, 100.0 - parent_similarity)
        within_child_distance = max(0.0, 100.0 - child_similarity)

        raw_w = 1.0 - math.exp(-parent_child_distance / ELASTIC_DISTANCE_SCALE)
        w = min(1.0, max(0.0, raw_w))

        parent_reads = getattr(node.parent, "raw_taxon_reads", node.parent.taxon_reads)
        node.threshold_func = "w = 1 - exp(-parent_child_distance / elastic_distance_scale)"
        node.variables = (
            f"parent_similarity={parent_similarity} child_similarity={child_similarity} "
            f"parent_child_distance={parent_child_distance} "
            f"within_child_distance={within_child_distance} raw_w={raw_w} "
            f"elastic_distance_scale={ELASTIC_DISTANCE_SCALE} "
            f"parent_taxon_reads={parent_reads}"
        )
        weighted_threshold = parent_reads * w

    else:
        w = 0.1
        raw_w = w
        parent_reads = getattr(node.parent, "raw_taxon_reads", node.parent.taxon_reads)
        node.threshold_func = "w = 0.1"
        node.variables = f"parent_taxon_reads={parent_reads}"
        weighted_threshold = parent_reads * w

    node.raw_lb_weight = raw_w
    node.lb_weight = w
    node.raw_weighted_threshold = weighted_threshold
    node.weighted_threshold = max(1.0, weighted_threshold)


def canonical_rank_index(taxon_level):
    """Convert Kraken rank labels to sortable rank indices."""
    ranks = ["R", "K", "D", "P", "C", "O", "F", "G", "S"]
    taxon_level = str(taxon_level)
    if taxon_level in ranks:
        return ranks.index(taxon_level)
    if len(taxon_level) > 1 and taxon_level[0] in ranks and taxon_level[1:].isdigit():
        return ranks.index(taxon_level[0]) + int(taxon_level[1:]) + 0.1
    return None


def rank_is_below(taxon_level, reference_rank):
    """Return True where ``taxon_level`` is below ``reference_rank``."""
    taxon_index = canonical_rank_index(taxon_level)
    reference_index = canonical_rank_index(reference_rank)
    return taxon_index is not None and reference_index is not None and taxon_index > reference_index


def normalise_fraction(value):
    """Accept fractions or percentages; 5 becomes 0.05."""
    value = float(value)
    if value > 1:
        return value / 100.0
    return value

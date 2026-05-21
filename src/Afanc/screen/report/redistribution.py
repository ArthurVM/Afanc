from pathlib import Path

from .scoring import normalise_fraction, score_lca_nodes
from .similarity import read_mvi_mash_distances, read_mvi_variant_index, similarity_between_subtrees


def redistribute_lca_hierarchical(
    root_node,
    mvi,
    family_rank="F",
    species_rank="S",
    scoring_ranks=("G", "S"),
    global_threshold=0.05,
    min_reads=0,
    mash_power=4.0,
    evidence_power=1.0,
    pseudocount=1.0,
    no_recipient_policy="unassigned",
    allow_evidence_fallback=True,
):
    """Redistribute non-scoring direct reads by LCA routing."""
    mvi_path = Path(mvi) if not isinstance(mvi, dict) else None
    similarity_lookup = read_mvi_mash_distances(mvi)
    variant_index = read_mvi_variant_index(mvi)

    root = root_node.root()
    root.snapshot_counts()
    root.redistribution_unassigned_reads = 0
    total_reads = root.raw_clade_reads
    threshold_input = global_threshold
    global_threshold = normalise_fraction(global_threshold)
    scoring_ranks = tuple(scoring_ranks)

    family_roots = root_node.find_by_taxon_level(family_rank)
    if root_node.taxon_level == family_rank and root_node not in family_roots:
        family_roots.insert(0, root_node)
    if not family_roots:
        family_roots = [root_node]

    audit = {
        "algorithm": "lca_hierarchical",
        "parameters": {
            "family_rank": family_rank,
            "species_rank": species_rank,
            "scoring_ranks": scoring_ranks,
            "global_threshold_input": threshold_input,
            "global_threshold": global_threshold,
            "min_reads": min_reads,
            "mvi": None if mvi_path is None else str(mvi_path),
            "metric": "mash_distance",
            "mash_similarity": "1 - mash_distance",
            "mash_power": mash_power,
            "evidence_power": evidence_power,
            "pseudocount": pseudocount,
            "no_recipient_policy": no_recipient_policy,
            "allow_evidence_fallback": allow_evidence_fallback,
            "subspecies_threshold": "elastic_parent_child_threshold",
        },
        "families": [],
        "total_deducted": 0,
        "total_added": 0,
        "total_unassigned": 0,
    }

    for family in family_roots:
        scoring_nodes = score_lca_nodes(
            family=family,
            scoring_ranks=scoring_ranks,
            species_rank=species_rank,
            global_threshold=global_threshold,
            min_reads=min_reads,
            total_reads=total_reads,
            variant_index=variant_index,
        )
        scoring_set = set(scoring_nodes)
        descendant_cache = {}

        family_record = {
            "family": family.name,
            "family_taxid": family.ncbi_taxID,
            "scoring_nodes": [
                {
                    "name": node.name,
                    "taxid": node.ncbi_taxID,
                    "rank": node.taxon_level,
                    "raw_clade_reads": node.raw_clade_reads,
                    "raw_fraction": node.raw_clade_reads / total_reads if total_reads else 0.0,
                    "scoring_rule": getattr(node, "scoring_rule", None),
                    "weighted_threshold": getattr(node, "weighted_threshold", None),
                }
                for node in scoring_nodes
            ],
            "donors": [],
        }

        donors = collect_direct_donors(
            family=family,
            scoring_set=scoring_set,
            min_reads=min_reads,
        )

        for donor in donors:
            donor_reads = int(round(donor.taxon_reads))
            if donor_reads <= 0:
                continue

            ancestor = redistribution_ancestor(
                donor=donor,
                family=family,
                scoring_set=scoring_set,
                descendant_cache=descendant_cache,
            )

            donor_record = {
                "donor": donor.name,
                "donor_taxid": donor.ncbi_taxID,
                "donor_rank": donor.taxon_level,
                "deducted_reads": donor_reads,
                "ancestor": None if ancestor is None else ancestor.name,
                "ancestor_taxid": None if ancestor is None else ancestor.ncbi_taxID,
                "routes": [],
            }

            donor.taxon_reads = 0
            audit["total_deducted"] += donor_reads

            if ancestor is None:
                if no_recipient_policy == "keep_at_family":
                    family.taxon_reads += donor_reads
                    audit["total_added"] += donor_reads
                    donor_record["unresolved"] = {
                        "target": family.name,
                        "target_taxid": family.ncbi_taxID,
                        "added_reads": donor_reads,
                    }
                else:
                    audit["total_unassigned"] += donor_reads
                    root.redistribution_unassigned_reads += donor_reads
                    donor_record["unassigned"] = {
                        "reads": donor_reads,
                        "reason": "no_scoring_node_within_family",
                    }
                family_record["donors"].append(donor_record)
                continue

            route_result = trickle_down(
                donor=donor,
                current=ancestor,
                reads=donor_reads,
                scoring_set=scoring_set,
                descendant_cache=descendant_cache,
                similarity_lookup=similarity_lookup,
                mash_power=mash_power,
                evidence_power=evidence_power,
                pseudocount=pseudocount,
                allow_evidence_fallback=allow_evidence_fallback,
                routes=donor_record["routes"],
            )

            audit["total_added"] += route_result["added"]
            audit["total_unassigned"] += route_result["unassigned"]
            root.redistribution_unassigned_reads += route_result["unassigned"]
            family_record["donors"].append(donor_record)

        root.recalculate_clade_counts(total_reads=total_reads)
        family_record["total_scoring_nodes"] = len(family_record["scoring_nodes"])
        family_record["total_donors"] = len(family_record["donors"])
        audit["families"].append(family_record)

    return audit


def collect_direct_donors(family, scoring_set, min_reads):
    """Collect direct-read bins which do not represent accepted signal."""
    donors = []
    for node in family.traverse():
        if node in scoring_set:
            continue
        if node.raw_taxon_reads <= 0 or node.raw_taxon_reads < min_reads:
            continue
        donors.append(node)
    donors.sort(key=lambda node: node.level_int, reverse=True)
    return donors


def redistribution_ancestor(donor, family, scoring_set, descendant_cache):
    """Find the nearest ancestor which has scoring signal in its subtree."""
    node = donor
    while node is not None:
        if node is not donor and node.parent is None:
            return None
        if scoring_descendants(node, scoring_set, descendant_cache):
            return node
        if node is family:
            return None
        node = node.parent
    return None


def scoring_descendants(node, scoring_set, descendant_cache):
    """Return scoring nodes contained inside ``node``."""
    key = id(node)
    if key not in descendant_cache:
        descendant_cache[key] = [child for child in node.traverse() if child in scoring_set]
    return descendant_cache[key]


def scoring_child_branches(node, scoring_set, descendant_cache):
    """Return child branches which contain at least one scoring node."""
    branches = []
    for child in node.children:
        if scoring_descendants(child, scoring_set, descendant_cache):
            branches.append(child)
    return branches


def trickle_down(
    donor,
    current,
    reads,
    scoring_set,
    descendant_cache,
    similarity_lookup,
    mash_power,
    evidence_power,
    pseudocount,
    allow_evidence_fallback,
    routes,
):
    """Recursively assign reads to scoring nodes through supported branches."""
    if reads <= 0:
        return {"added": 0, "unassigned": 0}

    child_branches = scoring_child_branches(current, scoring_set, descendant_cache)
    if not child_branches:
        if current in scoring_set:
            current.taxon_reads += reads
            routes.append(
                {
                    "target": current.name,
                    "target_taxid": current.ncbi_taxID,
                    "target_rank": current.taxon_level,
                    "added_reads": reads,
                    "reason": "terminal_scoring_node",
                }
            )
            return {"added": reads, "unassigned": 0}

        routes.append(
            {
                "target": current.name,
                "target_taxid": current.ncbi_taxID,
                "target_rank": current.taxon_level,
                "unassigned_reads": reads,
                "reason": "no_scoring_child_branch",
            }
        )
        return {"added": 0, "unassigned": reads}

    weighted_branches = weight_branches(
        donor=donor,
        branches=child_branches,
        similarity_lookup=similarity_lookup,
        mash_power=mash_power,
        evidence_power=evidence_power,
        pseudocount=pseudocount,
        allow_evidence_fallback=allow_evidence_fallback,
    )
    allocations = integer_allocations(reads, weighted_branches)
    if not allocations:
        routes.append(
            {
                "target": current.name,
                "target_taxid": current.ncbi_taxID,
                "target_rank": current.taxon_level,
                "unassigned_reads": reads,
                "reason": "no_weighted_scoring_branch",
            }
        )
        return {"added": 0, "unassigned": reads}

    added = 0
    unassigned = 0
    for branch, metric, weight, allocated_reads in allocations:
        route_record = {
            "from": current.name,
            "from_taxid": current.ncbi_taxID,
            "to": branch.name,
            "to_taxid": branch.ncbi_taxID,
            "to_rank": branch.taxon_level,
            "allocated_reads": allocated_reads,
            "similarity": None if metric is None else metric["similarity"],
            "distance": None if metric is None else metric["distance"],
            "metric": metric,
            "weight": weight,
        }
        routes.append(route_record)
        result = trickle_down(
            donor=donor,
            current=branch,
            reads=allocated_reads,
            scoring_set=scoring_set,
            descendant_cache=descendant_cache,
            similarity_lookup=similarity_lookup,
            mash_power=mash_power,
            evidence_power=evidence_power,
            pseudocount=pseudocount,
            allow_evidence_fallback=allow_evidence_fallback,
            routes=routes,
        )
        route_record["descendant_added"] = result["added"]
        route_record["descendant_unassigned"] = result["unassigned"]
        added += result["added"]
        unassigned += result["unassigned"]

    return {"added": added, "unassigned": unassigned}


def weight_branches(
    donor,
    branches,
    similarity_lookup,
    mash_power,
    evidence_power,
    pseudocount,
    allow_evidence_fallback,
):
    """Weight scoring child branches with taxon distance and read evidence."""
    weighted = []
    for branch in branches:
        metric = similarity_between_subtrees(donor, branch, similarity_lookup)
        evidence = max(float(getattr(branch, "raw_clade_reads", branch.clade_reads)), pseudocount)

        if metric is None:
            if not allow_evidence_fallback:
                continue
            weight = evidence ** evidence_power
        else:
            weight = (metric["similarity"] ** mash_power) * (evidence ** evidence_power)

        if weight > 0:
            weighted.append((branch, metric, weight))

    return weighted


def integer_allocations(reads, weighted_recipients):
    if reads <= 0 or not weighted_recipients:
        return []

    total_weight = sum(weight for _, _, weight in weighted_recipients)
    if total_weight <= 0:
        return []

    raw_allocations = []
    allocated = 0
    for recipient, metric, weight in weighted_recipients:
        exact = reads * (weight / total_weight)
        whole = int(exact)
        allocated += whole
        raw_allocations.append([recipient, metric, weight, whole, exact - whole])

    for row in sorted(raw_allocations, key=lambda item: item[4], reverse=True)[: reads - allocated]:
        row[3] += 1

    return [(recipient, metric, weight, whole) for recipient, metric, weight, whole, _ in raw_allocations if whole > 0]

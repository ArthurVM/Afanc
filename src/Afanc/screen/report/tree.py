import csv
import json
import math
from pathlib import Path


class Tree(object):
    'Tree node.'
    ELASTIC_DISTANCE_SCALE = 0.0015

    def __init__(self, line, name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID, children=None, parent=None):

        self.line = line

        ## initialise empty attributes
        self.name = name
        self.level_int = level_int
        self.clade_perc = clade_perc
        self.clade_reads = clade_reads
        self.taxon_reads = taxon_reads
        self.taxon_level = taxon_level
        self.ncbi_taxID = ncbi_taxID
        self.probability = None

        ## parent-child structure
        self.children = []
        self.parent = parent

        if children is not None:
            for child in children:
                self.add_child(child)


    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)


    def makeJsonLine(self, dbdict):
        """ generate a JSON output field for a node
        """
        if str(self.ncbi_taxID) in dbdict:
            acc = dbdict[str(self.ncbi_taxID)][1]
        elif self.ncbi_taxID in dbdict:
            acc = dbdict[self.ncbi_taxID][1]
        else: acc = "NA"

        json_data = {
            "name" : self.name,
            "reads" : self.taxon_reads,
            "percentage" : self.clade_perc,
            "taxon_id" : self.ncbi_taxID,
            "accession" : acc,
            }

        if hasattr(self, "weighted_threshold"):
            json_data["local_threshold_calc"] = {
                "threshold_func" : self.threshold_func,
                "variables" : self.variables,
                "raw_lower_bound_weighting" : getattr(self, "raw_lb_weight", self.lb_weight),
                "lower_bound_weighting" : self.lb_weight,
                "raw_weighted_threshold" : getattr(self, "raw_weighted_threshold", self.weighted_threshold),
                "weighted_threshold" : self.weighted_threshold
                 }

        return json_data


    def redistribute_lca_hierarchical(
        self,
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
        """ Redistribute non-scoring direct reads by LCA routing.

        This is a tree-native alternative to the older species-pass/child-pass
        redistribution model. Genus and species nodes are called with the same
        global threshold. Below-species nodes are called with the elastic
        parent-child threshold. Direct reads assigned to non-scoring nodes are
        moved to the nearest ancestor whose subtree contains a scoring node,
        then trickled back down through scoring branches using taxon distance
        weights.

        The key idea is that reads should only be commuted as far up the tree as
        necessary to find a local scoring context. Once at this context, reads
        are trickled back down the tree through scoring branches. This avoids
        forcing genus-level ambiguity directly into species-level bins, and
        avoids moving reads across distant taxonomic branches where no local
        evidence exists.

        Flow:

            non-scoring donor
                    |
                    v
            nearest ancestor with scoring signal
                    |
                    v
            weighted split between scoring child branches
                    |
                    v
            recurse until no scoring descendant remains
                    |
                    v
            assign to scoring genus/species/sub-species node

        Rules:

        - The donor unit is direct ``taxon_reads``.
        - Non-scoring internal nodes donate only their direct reads.
        - Scoring descendants are not zeroed by a non-scoring ancestor.
        - Reads are only finally assigned to scoring taxa.
        - Internal scoring taxa may retain reads when they have no scoring
          descendants.

        e.g.

                    F
                    |
              -------------
              |           |
             G1          G2*
              |
          ---------
          |       |
         S1      S2*
          |
         S1.1

        where * indicates a scoring node. Direct reads assigned to non-scoring
        S1.1 are commuted to the nearest ancestor with scoring signal (G1 or F,
        depending on the surrounding calls), then redistributed down into the
        scoring branch S2 rather than globally across the full database.
        """
        ## read MVI inputs
        mvi_path = Path(mvi) if not isinstance(mvi, dict) else None
        similarity_lookup = self._read_mvi_mash_distances(mvi)
        variant_index = self._read_mvi_variant_index(mvi)

        ## snapshot original counts before mutating direct assignments
        root = self.root()
        root.snapshot_counts()
        root.redistribution_unassigned_reads = 0
        total_reads = root.raw_clade_reads
        threshold_input = global_threshold
        global_threshold = self._normalise_fraction(global_threshold)
        scoring_ranks = tuple(scoring_ranks)

        ## constrain redistribution to local family subtrees
        family_roots = self.find_by_taxon_level(family_rank)
        if self.taxon_level == family_rank and self not in family_roots:
            family_roots.insert(0, self)
        if not family_roots:
            family_roots = [self]

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
            ## identify the nodes which are allowed to retain final read assignments
            scoring_nodes = self._lca_score_nodes(
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

            ## initialise a per-family audit block
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

            donors = self._lca_collect_direct_donors(
                family=family,
                scoring_set=scoring_set,
                min_reads=min_reads,
            )

            ## commute each non-scoring direct-read bin independently
            for donor in donors:
                donor_reads = int(round(donor.taxon_reads))
                if donor_reads <= 0:
                    continue

                ## find the closest upstream node that can route reads to scoring signal
                ancestor = self._lca_redistribution_ancestor(
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

                ## remove direct reads from the non-scoring donor
                donor.taxon_reads = 0
                audit["total_deducted"] += donor_reads

                ## if there is no local scoring context, reads cannot be defensibly assigned
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

                ## trickle the donor reads back down through scoring child branches
                route_result = self._lca_trickle_down(
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

            ## rebuild clade counts after direct assignments have been altered
            root.recalculate_clade_counts(total_reads=total_reads)
            family_record["total_scoring_nodes"] = len(family_record["scoring_nodes"])
            family_record["total_donors"] = len(family_record["donors"])
            audit["families"].append(family_record)

        return audit


    def _lca_score_nodes(self, family, scoring_ranks, species_rank, global_threshold, min_reads, total_reads, variant_index):
        """ Score nodes inside a family subtree.

        Genus and species level nodes are scored against the global threshold.
        Nodes below species are not assessed globally, since these are expected
        to be sparsely populated by Kraken2. Instead, below-species nodes are
        scored only where their species-level ancestor scores globally, and are
        then tested against an elastic parent-child threshold.

        This gives:

            G / S     : global threshold
            S1 / S2   : elastic threshold, conditional on parent species call
        """
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
            if not self._rank_is_below(node.taxon_level, species_rank):
                continue
            if node.raw_clade_reads <= 0 or node.raw_clade_reads < min_reads:
                continue

            species = node.ancestor_at_taxon_level(species_rank)
            if species not in scoring_species:
                continue

            node._elastic_threshold(variant_index)
            if node.raw_clade_reads >= node.weighted_threshold:
                node.scoring_rule = "elastic_threshold"
                scoring.append(node)

        return scoring


    def _lca_collect_direct_donors(self, family, scoring_set, min_reads):
        """ Collect direct-read bins which do not represent accepted signal.

        Importantly, this function only looks at ``taxon_reads``. An internal
        non-scoring node can donate the reads assigned directly to it without
        disturbing any scoring descendants in its subtree.
        """
        donors = []
        for node in family.traverse():
            if node in scoring_set:
                continue
            if node.raw_taxon_reads <= 0 or node.raw_taxon_reads < min_reads:
                continue
            donors.append(node)
        donors.sort(key=lambda node: node.level_int, reverse=True)
        return donors


    def _lca_redistribution_ancestor(self, donor, family, scoring_set, descendant_cache):
        """ Find the nearest ancestor which has scoring signal in its subtree.

        The donor reads are not sent to the root of the whole tree. They are
        commuted only to the first local taxonomic context from which there is a
        scoring route back down. If no such context exists inside the family,
        the reads remain unresolved.
        """
        node = donor
        while node is not None:
            if node is not donor and node.parent is None:
                return None
            if self._lca_scoring_descendants(node, scoring_set, descendant_cache):
                return node
            if node is family:
                return None
            node = node.parent
        return None


    def _lca_scoring_descendants(self, node, scoring_set, descendant_cache):
        """ Return scoring nodes contained inside ``node``.

        This is cached because the same descendant query is used repeatedly
        while routing each donor through the tree.
        """
        key = id(node)
        if key not in descendant_cache:
            descendant_cache[key] = [child for child in node.traverse() if child in scoring_set]
        return descendant_cache[key]


    def _lca_scoring_child_branches(self, node, scoring_set, descendant_cache):
        """ Return child branches which contain at least one scoring node."""
        branches = []
        for child in node.children:
            if self._lca_scoring_descendants(child, scoring_set, descendant_cache):
                branches.append(child)
        return branches


    def _lca_trickle_down(
        self,
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
        """ Recursively assign reads to scoring nodes through supported branches.

        At each node, reads are split between child branches which contain
        scoring signal. If there are no scoring child branches, reads stop at
        the current node only if the current node itself scores.

        This allows internal Kraken2 assignments to remain at scoring internal
        taxa, rather than forcing all reads to terminal tips.
        """
        if reads <= 0:
            return {"added": 0, "unassigned": 0}

        ## identify routes from this node to lower scoring calls
        child_branches = self._lca_scoring_child_branches(current, scoring_set, descendant_cache)
        if not child_branches:
            ## terminal scoring context: assign reads here
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

            ## no scoring route remains
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

        ## weight each scoring child branch and split reads by integer allocation
        weighted_branches = self._lca_weight_branches(
            donor=donor,
            branches=child_branches,
            similarity_lookup=similarity_lookup,
            mash_power=mash_power,
            evidence_power=evidence_power,
            pseudocount=pseudocount,
            allow_evidence_fallback=allow_evidence_fallback,
        )
        allocations = self._integer_allocations(reads, weighted_branches)
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
            ## record this branch-level decision before recursing down the tree
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
            result = self._lca_trickle_down(
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


    def _lca_weight_branches(
        self,
        donor,
        branches,
        similarity_lookup,
        mash_power,
        evidence_power,
        pseudocount,
        allow_evidence_fallback,
    ):
        """ Weight scoring child branches with taxon distance and read evidence.

        The matrix lookup is taxID based. Where no distance is available, the
        function can optionally fall back to abundance-only weighting. This is
        useful while testing incomplete MVI matrices, but should be treated as a
        lower confidence route in the audit output.
        """
        weighted = []
        for branch in branches:
            metric = self._similarity_between_subtrees(donor, branch, similarity_lookup)
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


    def _similarity_between_subtrees(self, donor, recipient, similarity_lookup):
        """ Return direct or best descendant-pair similarity metric for two nodes.

        This supports comparing internal taxa even where the matrix only
        contains lower-rank leaves. A direct taxID match is preferred. If this
        is absent, the highest similarity between descendant leaves is used.
        """
        donor_taxid = int(donor.ncbi_taxID)
        recipient_taxid = int(recipient.ncbi_taxID)

        if donor_taxid == recipient_taxid:
            return {"similarity": 1.0, "distance": 0.0, "source": "self"}

        direct = similarity_lookup.get((donor_taxid, recipient_taxid))
        if direct is not None:
            return direct

        metrics = []
        for donor_leaf in donor.leaf_descendants():
            for recipient_leaf in recipient.leaf_descendants():
                metric = similarity_lookup.get((int(donor_leaf.ncbi_taxID), int(recipient_leaf.ncbi_taxID)))
                if metric is not None:
                    metrics.append(metric)

        if metrics:
            return max(metrics, key=lambda metric: metric["similarity"])

        return None


    def snapshot_counts(self):
        """ Store pre-redistribution counts for every node below this node."""
        for node in self.traverse():
            node.raw_taxon_reads = node.taxon_reads
            node.raw_clade_reads = node.clade_reads
            node.raw_clade_perc = node.clade_perc


    def find_by_taxon_level(self, taxon_level):
        """ Return descendants whose Kraken taxon_level exactly matches ``taxon_level``."""
        return [node for node in self.traverse() if node.taxon_level == taxon_level]


    def ancestor_at_taxon_level(self, taxon_level):
        """ Return this node or its nearest ancestor at ``taxon_level``."""
        node = self
        while node is not None:
            if node.taxon_level == taxon_level:
                return node
            node = node.parent
        return None


    def _elastic_threshold(self, variant_index):
        """ Calculate the parent-child elastic read threshold for this node.

        Kraken2 assigns reads to the lowest taxonomic node supported by the
        database. Where a child is very similar to its parent, only a small
        fraction of reads are expected to resolve specifically to the child.
        The elastic threshold therefore scales with parent-child distance.

        The denominator here is the number of reads assigned directly to the
        parent, not the parent clade total:

            parent_clade_reads = parent_taxon_reads + child_clade_reads

        Using parent clade reads would make abundant child calls
        self-penalising. The direct parent bin is a better approximation of the
        ambiguous read pool.
        """
        if self.parent is None:
            w = 0.0
            raw_w = w
            weighted_threshold = 0.0
            self.threshold_func = "root"
            self.variables = None

        elif str(self.ncbi_taxID) in variant_index:
            ## calculate absolute parent-child Mash distance on the 0-100 scale
            parent_similarity = variant_index[str(self.ncbi_taxID)]["parent_index"]["mean"]
            child_similarity = variant_index[str(self.ncbi_taxID)]["sibling_index"]["mean"]
            parent_child_distance = max(0.0, 100.0 - parent_similarity)
            within_child_distance = max(0.0, 100.0 - child_similarity)

            ## saturating distance transform: close children require fewer reads
            raw_w = 1.0 - math.exp(-parent_child_distance / self.ELASTIC_DISTANCE_SCALE)
            w = min(1.0, max(0.0, raw_w))

            ## Use direct parent assignments as the ambiguous-read denominator.
            ## The parent clade total includes this child, which makes abundant
            ## sub-species calls self-penalising and can erase true sparse tips.
            parent_reads = getattr(self.parent, "raw_taxon_reads", self.parent.taxon_reads)
            self.threshold_func = "w = 1 - exp(-parent_child_distance / elastic_distance_scale)"
            self.variables = (
                f"parent_similarity={parent_similarity} child_similarity={child_similarity} "
                f"parent_child_distance={parent_child_distance} "
                f"within_child_distance={within_child_distance} raw_w={raw_w} "
                f"elastic_distance_scale={self.ELASTIC_DISTANCE_SCALE} "
                f"parent_taxon_reads={parent_reads}"
            )
            weighted_threshold = parent_reads * w

        else:
            ## if the child is absent from the MVI, fall back to a conservative
            ## direct-parent threshold
            w = 0.1
            raw_w = w
            parent_reads = getattr(self.parent, "raw_taxon_reads", self.parent.taxon_reads)
            self.threshold_func = "w = 0.1"
            self.variables = f"parent_taxon_reads={parent_reads}"
            weighted_threshold = parent_reads * w

        self.raw_lb_weight = raw_w
        self.lb_weight = w
        self.raw_weighted_threshold = weighted_threshold
        self.weighted_threshold = max(1.0, weighted_threshold)


    @staticmethod
    def _canonical_rank_index(taxon_level):
        """ Convert Kraken rank labels to sortable rank indices."""
        ranks = ["R", "K", "D", "P", "C", "O", "F", "G", "S"]
        taxon_level = str(taxon_level)
        if taxon_level in ranks:
            return ranks.index(taxon_level)
        if len(taxon_level) > 1 and taxon_level[0] in ranks and taxon_level[1:].isdigit():
            return ranks.index(taxon_level[0]) + int(taxon_level[1:]) + 0.1
        return None


    def _rank_is_below(self, taxon_level, reference_rank):
        """ Return True where ``taxon_level`` is below ``reference_rank``."""
        taxon_index = self._canonical_rank_index(taxon_level)
        reference_index = self._canonical_rank_index(reference_rank)
        return taxon_index is not None and reference_index is not None and taxon_index > reference_index


    def leaf_descendants(self):
        """ Return leaf descendants below this node, including self if it is a leaf."""
        return [node for node in self.traverse() if not node.children]


    def recalculate_clade_counts(self, total_reads=None):
        """ Recompute clade_reads and clade_perc from taxon_reads bottom-up."""
        root = self.root()
        recalculated_total = root._recalculate_clade_reads()
        if total_reads is None:
            total_reads = recalculated_total

        root._set_clade_perc(total_reads)
        return total_reads


    def root(self):
        node = self
        while node.parent is not None:
            node = node.parent
        return node


    def to_kraken_report(self, include_unassigned=True, prune_zero=False):
        """ Generate Kraken-like report text from the current tree state."""
        root = self.root()
        lines = []

        unassigned = int(round(getattr(root, "redistribution_unassigned_reads", 0)))
        if include_unassigned and unassigned > 0:
            denominator = getattr(root, "raw_clade_reads", root.clade_reads + unassigned)
            perc = (unassigned / denominator) * 100 if denominator else 0.0
            lines.append("\t".join([f"{perc:.2f}", str(unassigned), str(unassigned), "U", "0", "unassigned_after_redistribution"]))

        lines.extend(root.to_kraken_report_lines(prune_zero=prune_zero))
        return "\n".join(lines) + "\n"


    def to_kraken_report_lines(self, prune_zero=False):
        """ Generate Kraken-like report lines from this node downward."""
        if prune_zero and self.clade_reads <= 0 and self.taxon_reads <= 0:
            return []

        indent = "  " * self.level_int
        lines = [
            "\t".join(
                [
                    f"{self.clade_perc:.2f}",
                    str(int(round(self.clade_reads))),
                    str(int(round(self.taxon_reads))),
                    str(self.taxon_level),
                    str(self.ncbi_taxID),
                    f"{indent}{self.name}",
                ]
            )
        ]

        for child in self.children:
            lines.extend(child.to_kraken_report_lines(prune_zero=prune_zero))

        return lines


    def write_kraken_report(self, outfile, include_unassigned=True, prune_zero=False):
        """ Write Kraken-like report text from the current tree state."""
        with open(outfile, "w") as fout:
            fout.write(self.to_kraken_report(include_unassigned=include_unassigned, prune_zero=prune_zero))


    def _recalculate_clade_reads(self):
        self.clade_reads = self.taxon_reads + sum(child._recalculate_clade_reads() for child in self.children)
        return self.clade_reads


    def _set_clade_perc(self, total_reads):
        if total_reads:
            self.clade_perc = (self.clade_reads / total_reads) * 100
        else:
            self.clade_perc = 0.0

        for child in self.children:
            child._set_clade_perc(total_reads)


    def traverse(self):
        """ Generator function for yielding all subnodes of a given node
        """
        for child in self.children:
            for grandchild in child.traverse():
                yield grandchild
        else:
            yield self


    @staticmethod
    def _normalise_fraction(value):
        """ Accept fractions or percentages; 5 becomes 0.05."""
        value = float(value)
        if value > 1:
            return value / 100.0
        return value


    @staticmethod
    def _read_mvi_variant_index(mvi):
        """ Read variant_index.json from an MVI directory or mapping."""
        if isinstance(mvi, dict):
            if "variant_index" in mvi:
                return mvi["variant_index"]
            return mvi

        mvi_path = Path(mvi)
        variant_index_path = mvi_path / "variant_index.json" if mvi_path.is_dir() else None
        if variant_index_path is None or not variant_index_path.exists():
            raise ValueError("variant_index was not supplied and variant_index.json was not found in the MVI.")

        with open(variant_index_path, "r") as fin:
            data = json.load(fin)

        return data["variant_index"]


    @staticmethod
    def _read_mvi_mash_distances(mvi):
        """ Read Mash distance matrices from an MVI path into a taxID-pair lookup.

        MVI input may be a TSV matrix, a directory containing Mash matrix TSVs,
        or a prebuilt metric lookup. Matrix values are native Mash distances and
        are converted to weighting similarities as ``1 - distance``.
        """
        if isinstance(mvi, dict):
            return mvi

        matrix_path = Path(mvi)
        if matrix_path.is_dir():
            matrix_files = [
                path for path in sorted(matrix_path.glob("*.tsv"))
                if "matrix" in path.name or "all_vs_all" in path.name
            ]
        else:
            matrix_files = [matrix_path]

        lookup = {}
        for matrix_file in matrix_files:
            Tree._read_mash_matrix_file(matrix_file, lookup)

        return lookup


    @staticmethod
    def _read_mash_matrix_file(matrix_file, lookup):
        with open(matrix_file, "r", newline="") as fin:
            reader = csv.reader(fin, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                return

            col_labels = header[1:]
            col_taxids = [Tree._taxid_from_matrix_label(label) for label in col_labels]

            for row in reader:
                if not row:
                    continue

                row_taxid = Tree._taxid_from_matrix_label(row[0])
                if row_taxid is None:
                    continue

                for col_taxid, value in zip(col_taxids, row[1:]):
                    if col_taxid is None or value in {None, ""}:
                        continue

                    raw_value = float(value)
                    metric = Tree._metric_from_mash_distance(raw_value, source=str(matrix_file))
                    Tree._store_best_metric(lookup, row_taxid, col_taxid, metric)


    @staticmethod
    def _taxid_from_matrix_label(label):
        label = str(label).strip()
        if not label:
            return None

        if "|" in label:
            label = label.split("|", 1)[0]
        else:
            label = Path(label).name
            label = label.split("_", 1)[0]

        try:
            return int(label)
        except ValueError:
            return None


    @staticmethod
    def _metric_from_mash_distance(value, source):
        distance = max(float(value), 0.0)
        similarity = max(1.0 - distance, 0.0)

        return {
            "similarity": similarity,
            "distance": distance,
            "source": source,
        }


    @staticmethod
    def _store_best_metric(lookup, taxid_a, taxid_b, metric):
        key = (int(taxid_a), int(taxid_b))
        reverse_key = (int(taxid_b), int(taxid_a))

        existing = lookup.get(key)
        if existing is None or metric["similarity"] > existing["similarity"]:
            lookup[key] = metric
            lookup[reverse_key] = metric


    @staticmethod
    def _integer_allocations(reads, weighted_recipients):
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

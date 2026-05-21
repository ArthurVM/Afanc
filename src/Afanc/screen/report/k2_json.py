import json


def get_scoring_nodes(root_node):
    """Return all nodes called as true signal by the tree deconvolution model."""
    return [node for node in root_node.traverse() if hasattr(node, "scoring_rule")]


def get_terminal_scoring_nodes(root_node):
    """Return scoring nodes which do not contain a lower scoring descendant."""
    scoring_nodes = set(get_scoring_nodes(root_node))
    terminal_nodes = []

    for node in scoring_nodes:
        scoring_descendants = [child for child in node.traverse() if child is not node and child in scoring_nodes]
        if not scoring_descendants:
            terminal_nodes.append(node)

    return sorted(terminal_nodes, key=lambda node: node.level_int)


def node_to_json_line(node, dbdict):
    """Generate a JSON output field for a node."""
    if str(node.ncbi_taxID) in dbdict:
        acc = dbdict[str(node.ncbi_taxID)][1]
    elif node.ncbi_taxID in dbdict:
        acc = dbdict[node.ncbi_taxID][1]
    else:
        acc = "NA"

    json_data = {
        "name" : node.name,
        "reads" : node.taxon_reads,
        "percentage" : node.clade_perc,
        "taxon_id" : node.ncbi_taxID,
        "accession" : acc,
        }

    if hasattr(node, "weighted_threshold"):
        json_data["local_threshold_calc"] = {
            "threshold_func" : node.threshold_func,
            "variables" : node.variables,
            "raw_lower_bound_weighting" : getattr(node, "raw_lb_weight", node.lb_weight),
            "lower_bound_weighting" : node.lb_weight,
            "raw_weighted_threshold" : getattr(node, "raw_weighted_threshold", node.weighted_threshold),
            "weighted_threshold" : node.weighted_threshold
             }

    return json_data


def makeJson(root_node, output_prefix, reportsDir, pct_threshold, num_threshold, dbdict, audit):
    """Generate a JSON report from terminal deconvolved tree calls."""
    out_json = f"{reportsDir}/{output_prefix}.k2.json"

    json_dict = {
        "Thresholds" : { "reads" : num_threshold, "percentage" : pct_threshold },
        "Deconvolution" : audit,
        "Detection_events" : []
        }

    for node in get_terminal_scoring_nodes(root_node):
        species_node = node.ancestor_at_taxon_level("S")

        if species_node is not None and node != species_node:
            json_line = node_to_json_line(species_node, dbdict)
            json_line["closest_variant"] = node_to_json_line(node, dbdict)
            json_dict["Detection_events"].append(json_line)

        else:
            json_dict["Detection_events"].append(node_to_json_line(node, dbdict))

    with open(out_json, "w") as fout:
        json.dump(json_dict, fout, indent = 4, default=str)

    return out_json

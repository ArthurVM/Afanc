""" Parse and deconvolve a kraken2 report using the Afanc tree model.
"""
import json
from collections import defaultdict
from os import path

from .tree import Tree


def read_variant_index(similarity_index):
    """ Reads the similarity index json
    """
    with open(similarity_index, 'r') as fin:
        variant_index = json.load(fin)

    return variant_index["variant_index"]


def parseK2line(line):
    """ parses a kraken2 report line
    """
    sline = line.strip("\n").split('\t')

    if len(sline) < 5:
        return []
    try:
        int(sline[1])
    except ValueError:
        return []

    #Extract relevant information
    clade_perc = float(sline[0])
    clade_reads =  int(sline[1])
    taxon_reads = int(sline[2])
    taxon_level = sline[-3]
    ncbi_taxID = int(sline[-2])

    #Get name and spaces
    spaces = 0
    name = sline[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break

    #Determine which level based on number of spaces
    level_int = int(spaces/2)

    return name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID


def readK2report(report):
    """ Read the kraken2 report and filter according to user defined values
    """

    resultsdict = defaultdict(list)

    main_lvls = ['R','K','D','P','C','O','F','G','S']

    ## initialise an empty root_node
    ## instances where the wrong database is used will throw an error since a root node may not be present in the kraken2 report
    root_node = None
    base_nodes = {}
    prev_node = -1

    with open(report, "r") as fin:

        for line in fin.readlines():
            parsed = parseK2line(line)
            if not parsed:
                continue

            name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID = parsed

            if name == "unclassified":
                continue

            ## handle tree root
            if ncbi_taxID == 1:
                root_node = Tree(line, name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID)
                prev_node = root_node

                base_nodes[ncbi_taxID] = root_node
                continue

            #move to correct parent
            while level_int != (prev_node.level_int + 1):
                prev_node = prev_node.parent

            #determine correct level ID
            if taxon_level == '-' or len(taxon_level) > 1:
                if prev_node.taxon_level in main_lvls:
                    taxon_level = prev_node.taxon_level + '1'
                else:
                    num = int(prev_node.taxon_level[-1]) + 1
                    taxon_level = prev_node.taxon_level[:-1] + str(num)

            #make node
            curr_node = Tree(line, name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID, None, prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node

            base_nodes[ncbi_taxID] = curr_node

    return base_nodes, root_node


def get_scoring_nodes(root_node):
    """ Return all nodes called as true signal by the tree deconvolution model.
    """
    return [node for node in root_node.traverse() if hasattr(node, "scoring_rule")]


def get_terminal_scoring_nodes(root_node):
    """ Return scoring nodes which do not contain a lower scoring descendant.

    These nodes represent the final deconvolved calls. Internal scoring nodes
    are retained only when no more specific scoring taxon is present below
    them.
    """
    scoring_nodes = set(get_scoring_nodes(root_node))
    terminal_nodes = []

    for node in scoring_nodes:
        scoring_descendants = [child for child in node.traverse() if child is not node and child in scoring_nodes]
        if not scoring_descendants:
            terminal_nodes.append(node)

    return sorted(terminal_nodes, key=lambda node: node.level_int)


def makeJson(root_node, output_prefix, reportsDir, pct_threshold, num_threshold, dbdict, audit):
    """ Generate a JSON report from terminal deconvolved tree calls.

    Below-species calls are reported as a species-level event with the called
    lower taxon captured in ``closest_variant``. This preserves the downstream
    mapping behaviour used by Afanc screen.
    """
    out_json = f"{reportsDir}/{output_prefix}.k2.json"

    json_dict = {
        "Thresholds" : { "reads" : num_threshold, "percentage" : pct_threshold },
        "Deconvolution" : audit,
        "Detection_events" : []
        }

    ## create json report dict
    for node in get_terminal_scoring_nodes(root_node):

        species_node = node.ancestor_at_taxon_level("S")

        ## if this is a below-species hit, keep the species as the mapping
        ## context and record the more specific taxon as the closest variant
        if species_node is not None and node != species_node:
            json_line = species_node.makeJsonLine(dbdict)
            json_line["closest_variant"] = node.makeJsonLine(dbdict)
            json_dict["Detection_events"].append(json_line)

        else:
            json_dict["Detection_events"].append(node.makeJsonLine(dbdict))

    with open(out_json, "w") as fout:
        json.dump(json_dict, fout, indent = 4, default=str)

    return out_json


def parseK2reportMain(args, dbdict):
    """ main function """
    report_path = f"{args.k2WDir}/{args.output_prefix}.k2.report.txt"

    if not path.exists(report_path):
        return None

    base_nodes, root_node = readK2report(report_path)

    ## escape if there i no root node present in the kraken2 report file
    if root_node == None:
        return None

    audit = root_node.redistribute_lca_hierarchical(
        mvi=getattr(args, "mvi_path", args.database),
        global_threshold=args.pct_threshold,
        min_reads=0,
    )

    commute_stats_path = path.join(args.k2WDir, f"{args.output_prefix}.commute_stats.json")
    filtered_report_path = path.join(args.k2WDir, f"{args.output_prefix}.filtered.k2.report.txt")

    with open(commute_stats_path, "w") as fout:
        json.dump(audit, fout, indent = 4, default=str)

    root_node.write_kraken_report(filtered_report_path, prune_zero=True)

    if len(get_terminal_scoring_nodes(root_node)) == 0:
        return None

    out_json = makeJson(root_node, args.output_prefix, args.reportsDir, args.pct_threshold, args.num_threshold, dbdict, audit)

    return out_json

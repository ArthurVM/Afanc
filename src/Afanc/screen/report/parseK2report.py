""" Uses Bayes rule to commute reads from taxa which fall below 40% of the calculated elastic threshold to taxa which dont.
"""
import json
from collections import defaultdict
from os import path, listdir

from Afanc.utilities.generalUtils import gendbdict
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
            name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID = parseK2line(line)

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


def get_local_max_list(branch_box):
    """ Generate a set of local max weighting nodes for each node in branch_box
    """

    best_hits = []
    for node in branch_box:
        top_hit = node.find_local_max()
        best_hits.append(top_hit)

    return set(best_hits)


def find_best_hit(root_node, variant_index, pct_threshold, num_threshold):
    """ Find the lowest level nodes on each branch which are weighted greater than the weight threshold,
    then generate a set of local max weighting nodes for each node in branch_box.
    """

    best_hits = []

    ## find lowest level scoring nodes
    for node in root_node.traverse():
        ## check if node scores above the pct_threshold and is a species complex level or lower node
        if node.clade_perc >= pct_threshold and node.clade_reads >= num_threshold and node.level_int >= 9:
            ## construct scorebox then strip off the head node
            scorebox = [ n.clade_perc for n in node.traverse() ][:-1]

            ## check if any node in the box exceeds the pct_threshold
            ## if so, the head node is not the lowest level node which exceeds pct_threshold on this branch
            if any(score >= pct_threshold for score in scorebox):
                continue

            ## else it is assumed to be the lowest level scoring node on this branch
            else:
                ## find the local max for this scoring node
                top_hit = node.find_local_max(variant_index)
                best_hits.append(top_hit)

        else:
            # print(colored(line, 'red'))
            continue

    return best_hits


def commute_reads(args, root_node, variant_index, pct_threshold, num_threshold):

    commute_stats = defaultdict(dict)

    ## find lowest level scoring nodes
    for node in root_node.traverse():
        ## check if node scores above the pct_threshold
        # print(node.name, node.level_int, node.clade_perc, node.clade_reads)
        if node.clade_perc >= pct_threshold and node.clade_reads >= num_threshold and node.level_int >= 9:
            ## construct scorebox then strip off the head node
            scorebox = [ n.clade_perc for n in node.traverse() ][:-1]

            ## check if any node in the box exceeds the pct_threshold
            ## if so, the head node is not the lowest level node which exceeds pct_threshold on this branch
            if any(score >= pct_threshold for score in scorebox):
                continue

            ## else it is assumed to be the lowest level scoring node on this branch
            else:
                ## find the local max for this scoring node
                commute_dict = node.bayes_commute(variant_index, args.lower_bound, args.upper_bound)
                commute_stats.update(commute_dict)

        else:
            continue

    with open(f"./{args.output_prefix}.commute_stats.json", "w") as fout:
        json.dump(commute_stats, fout, indent = 4)


def makeJson(branch_box, output_prefix, reportsDir, pct_threshold, num_threshold, dbdict):
    """ takes the results dict and generates a json report.
    {
    F : [tax1, tax2, ..., taxn],
    G : [tax1, tax2, ..., taxn],
    S : [tax1, tax2, ..., taxn],
    S1 : [tax1, tax2, ..., taxn],
    ...,
    Sn : [tax1, tax2, ..., taxn]
    }
    """

    taxon_level_key = { \
    "P" : "Phylum",
    "C" : "Class",
    "O" : "Order",
    "F" : "Family",
    "G" : "Genus",
    "G1": "Species Complex",
    "S" : "Species"
     }

    out_json = f"{reportsDir}/{output_prefix}.k2.json"

    json_dict = { "Thresholds" : { "reads" : num_threshold, "percentage" : pct_threshold }, "Detection_events" : [] }

    ## create json report dict
    for node in branch_box:

        ## check if the node is its own mother_clade, and therefore has no scoring subclades
        if node != node.mother_clade:
            json_line = node.mother_clade.makeJsonLine(dbdict)
            json_line["closest_variant"] = node.makeJsonLine(dbdict)
            json_dict["Detection_events"].append(json_line)

        else:
            json_dict["Detection_events"].append(node.makeJsonLine(dbdict))

    with open(out_json, "w") as fout:
        json.dump(json_dict, fout, indent = 4)

    return out_json


def parseK2reportMain(args, dbdict):
    """ main function """
    report_path = f"{args.k2WDir}/{args.output_prefix}.k2.report.txt"

    if not path.exists(report_path):
        return None

    variant_index = read_variant_index(args.variant_index_path)

    base_nodes, root_node = readK2report(report_path)

    ## escape if there i no root node present in the kraken2 report file
    if root_node == None:
        return None

    commute_reads(args, root_node, variant_index, args.pct_threshold, args.num_threshold)

    best_hits = find_best_hit(root_node, variant_index, args.pct_threshold, args.num_threshold)

    if len(best_hits) == 0:
        return None

    out_json = makeJson(best_hits, args.output_prefix, args.reportsDir, args.pct_threshold, args.num_threshold, dbdict)

    return out_json

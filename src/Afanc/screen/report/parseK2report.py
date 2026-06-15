""" Parse and deconvolve a kraken2 report using the Afanc tree model.
"""
import json
from os import path

from .k2_json import get_scoring_nodes, get_terminal_scoring_nodes, makeJson
from .kraken_tree_io import parseK2line, readK2report, write_kraken_report
from .redistribution import redistribute_lca_hierarchical


def read_variant_index(similarity_index):
    """ Reads the similarity index json
    """
    with open(similarity_index, 'r') as fin:
        variant_index = json.load(fin)

    return variant_index["variant_index"]


def parseK2reportMain(args, dbdict):
    """ main function """
    report_path = f"{args.k2WDir}/{args.output_prefix}.k2.report.txt"

    if not path.exists(report_path):
        return None

    base_nodes, root_node = readK2report(report_path)

    ## escape if there i no root node present in the kraken2 report file
    if root_node == None:
        return None

    audit = redistribute_lca_hierarchical(
        root_node,
        mvi=getattr(args, "mvi_path", args.database),
        global_threshold=args.pct_threshold,
        min_reads=0,
    )

    commute_stats_path = path.join(args.k2WDir, f"{args.output_prefix}.commute_stats.json")
    filtered_report_path = path.join(args.k2WDir, f"{args.output_prefix}.filtered.k2.report.txt")

    with open(commute_stats_path, "w") as fout:
        json.dump(audit, fout, indent = 4, default=str)

    write_kraken_report(root_node, filtered_report_path, prune_zero=True)

    if len(get_terminal_scoring_nodes(root_node)) == 0:
        return None

    out_json = makeJson(root_node, args.output_prefix, args.reportsDir, args.pct_threshold, args.num_threshold, dbdict, audit)

    return out_json

from collections import defaultdict

from .tree import Tree


def parseK2line(line):
    """Parse a Kraken2 report line."""
    sline = line.strip("\n").split('\t')

    if len(sline) < 5:
        return []
    try:
        int(sline[1])
    except ValueError:
        return []

    clade_perc = float(sline[0])
    clade_reads = int(sline[1])
    taxon_reads = int(sline[2])
    taxon_level = sline[-3]
    ncbi_taxID = int(sline[-2])

    spaces = 0
    name = sline[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break

    level_int = int(spaces / 2)

    return name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID


def readK2report(report):
    """Read a Kraken2 report into a tree."""
    resultsdict = defaultdict(list)
    main_lvls = ['R', 'K', 'D', 'P', 'C', 'O', 'F', 'G', 'S']

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

            if ncbi_taxID == 1:
                root_node = Tree(line, name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID)
                prev_node = root_node
                base_nodes[ncbi_taxID] = root_node
                continue

            while level_int != (prev_node.level_int + 1):
                prev_node = prev_node.parent

            if taxon_level == '-' or len(taxon_level) > 1:
                if prev_node.taxon_level in main_lvls:
                    taxon_level = prev_node.taxon_level + '1'
                else:
                    num = int(prev_node.taxon_level[-1]) + 1
                    taxon_level = prev_node.taxon_level[:-1] + str(num)

            curr_node = Tree(line, name, level_int, clade_perc, clade_reads, taxon_reads, taxon_level, ncbi_taxID, None, prev_node)
            prev_node.add_child(curr_node)
            prev_node = curr_node

            base_nodes[ncbi_taxID] = curr_node

    return base_nodes, root_node


def to_kraken_report(root_node, include_unassigned=True, prune_zero=False):
    """Generate Kraken-like report text from the current tree state."""
    root = root_node.root()
    lines = []

    unassigned = int(round(getattr(root, "redistribution_unassigned_reads", 0)))
    if include_unassigned and unassigned > 0:
        denominator = getattr(root, "raw_clade_reads", root.clade_reads + unassigned)
        perc = (unassigned / denominator) * 100 if denominator else 0.0
        lines.append("\t".join([f"{perc:.2f}", str(unassigned), str(unassigned), "U", "0", "unassigned_after_redistribution"]))

    lines.extend(to_kraken_report_lines(root, prune_zero=prune_zero))
    return "\n".join(lines) + "\n"


def to_kraken_report_lines(node, prune_zero=False):
    """Generate Kraken-like report lines from this node downward."""
    if prune_zero and node.clade_reads <= 0 and node.taxon_reads <= 0:
        return []

    indent = "  " * node.level_int
    lines = [
        "\t".join(
            [
                f"{node.clade_perc:.2f}",
                str(int(round(node.clade_reads))),
                str(int(round(node.taxon_reads))),
                str(node.taxon_level),
                str(node.ncbi_taxID),
                f"{indent}{node.name}",
            ]
        )
    ]

    for child in node.children:
        lines.extend(to_kraken_report_lines(child, prune_zero=prune_zero))

    return lines


def write_kraken_report(root_node, outfile, include_unassigned=True, prune_zero=False):
    """Write Kraken-like report text from the current tree state."""
    with open(outfile, "w") as fout:
        fout.write(to_kraken_report(root_node, include_unassigned=include_unassigned, prune_zero=prune_zero))

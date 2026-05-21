from Afanc.screen.report.k2_json import get_terminal_scoring_nodes
from Afanc.screen.report.redistribution import redistribute_lca_hierarchical
from Afanc.screen.report.tree import Tree


def _add_child(parent, name, level_int, clade_reads, taxon_reads, rank, taxid):
    child = Tree("", name, level_int, 0.0, clade_reads, taxon_reads, rank, taxid, parent=parent)
    parent.add_child(child)
    return child


def test_lca_redistribution_remains_available_after_tree_split():
    root = Tree("", "root", 0, 100.0, 100, 0, "R", 1)
    family = _add_child(root, "family", 1, 100, 5, "F", 10)
    genus = _add_child(family, "genus", 2, 95, 30, "G", 20)
    species = _add_child(genus, "species", 3, 65, 65, "S", 30)

    audit = redistribute_lca_hierarchical(
        root,
        mvi={"variant_index": {}},
        global_threshold=0.05,
    )

    assert audit["algorithm"] == "lca_hierarchical"
    assert audit["total_deducted"] == 5
    assert audit["total_added"] == 5
    assert species.taxon_reads == 70
    assert get_terminal_scoring_nodes(root) == [species]


def test_tree_compatibility_wrapper_delegates_to_redistribution():
    root = Tree("", "root", 0, 100.0, 100, 0, "R", 1)
    family = _add_child(root, "family", 1, 100, 5, "F", 10)
    genus = _add_child(family, "genus", 2, 95, 30, "G", 20)
    species = _add_child(genus, "species", 3, 65, 65, "S", 30)

    root.redistribute_lca_hierarchical(
        mvi={"variant_index": {}},
        global_threshold=0.05,
    )

    assert species.taxon_reads == 70

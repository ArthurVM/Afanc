class Tree(object):
    """Kraken taxonomy tree node.

    This class intentionally contains only tree state and traversal/counting
    helpers. Parsing, scoring, redistribution, similarity loading, and report
    rendering live in sibling modules.
    """

    def __init__(
        self,
        line,
        name,
        level_int,
        clade_perc,
        clade_reads,
        taxon_reads,
        taxon_level,
        ncbi_taxID,
        children=None,
        parent=None,
    ):
        self.line = line
        self.name = name
        self.level_int = level_int
        self.clade_perc = clade_perc
        self.clade_reads = clade_reads
        self.taxon_reads = taxon_reads
        self.taxon_level = taxon_level
        self.ncbi_taxID = ncbi_taxID
        self.probability = None

        self.children = []
        self.parent = parent

        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)

    def snapshot_counts(self):
        """Store pre-redistribution counts for every node below this node."""
        for node in self.traverse():
            node.raw_taxon_reads = node.taxon_reads
            node.raw_clade_reads = node.clade_reads
            node.raw_clade_perc = node.clade_perc

    def find_by_taxon_level(self, taxon_level):
        """Return descendants whose Kraken taxon_level exactly matches ``taxon_level``."""
        return [node for node in self.traverse() if node.taxon_level == taxon_level]

    def ancestor_at_taxon_level(self, taxon_level):
        """Return this node or its nearest ancestor at ``taxon_level``."""
        node = self
        while node is not None:
            if node.taxon_level == taxon_level:
                return node
            node = node.parent
        return None

    def leaf_descendants(self):
        """Return leaf descendants below this node, including self if it is a leaf."""
        return [node for node in self.traverse() if not node.children]

    def recalculate_clade_counts(self, total_reads=None):
        """Recompute clade_reads and clade_perc from taxon_reads bottom-up."""
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

    def traverse(self):
        """Yield all subnodes below this node, including self last."""
        for child in self.children:
            for grandchild in child.traverse():
                yield grandchild
        else:
            yield self

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

    def redistribute_lca_hierarchical(self, *args, **kwargs):
        """Compatibility wrapper for the LCA hierarchical redistribution strategy."""
        from .redistribution import redistribute_lca_hierarchical

        return redistribute_lca_hierarchical(self, *args, **kwargs)

    def to_kraken_report(self, *args, **kwargs):
        """Compatibility wrapper for Kraken-like report rendering."""
        from .kraken_tree_io import to_kraken_report

        return to_kraken_report(self, *args, **kwargs)

    def to_kraken_report_lines(self, *args, **kwargs):
        """Compatibility wrapper for Kraken-like report line rendering."""
        from .kraken_tree_io import to_kraken_report_lines

        return to_kraken_report_lines(self, *args, **kwargs)

    def write_kraken_report(self, *args, **kwargs):
        """Compatibility wrapper for Kraken-like report writing."""
        from .kraken_tree_io import write_kraken_report

        return write_kraken_report(self, *args, **kwargs)

    def makeJsonLine(self, *args, **kwargs):
        """Compatibility wrapper for detection-event JSON rendering."""
        from .k2_json import node_to_json_line

        return node_to_json_line(self, *args, **kwargs)

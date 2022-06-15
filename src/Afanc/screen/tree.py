class Tree(object):
    'Tree node.'
    def __init__(self, line, name, level_int, clade_perc, clade_reads, taxa_reads, taxa_level, ncbi_taxID, children=None, parent=None):

        self.line = line

        ## initialise empty attributes
        self.name = name
        self.level_int = level_int
        self.clade_perc = clade_perc
        self.clade_reads = clade_reads
        self.taxa_reads = taxa_reads
        self.taxa_level = taxa_level
        self.ncbi_taxID = ncbi_taxID

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
        if self.ncbi_taxID in dbdict:
            acc = dbdict[self.ncbi_taxID][1]
        else: acc = "NA"

        json_data = { "name" : self.name, "reads" : self.taxa_reads, "percentage" : self.clade_perc, "taxon_id" : self.ncbi_taxID, "accession" : acc }

        return json_data

    def find_local_max(self, local_threshold, return_all_nodes=False):
        """ Finds the local maximum scoring tip node. Takes a branch node which is assumed to be the lowest level node
        which weights higher than some threshold (clade percentage threshold), and returns the tip node (not necessarily
        the lowest level node) with the highest score.

        e.g.
                N  -  -  -  -  -  0
              /   \
             n0   n1  -  -  -  -  1
            /  \    \
           n2  n3   n4   -  -  -  2
          /  \        \
         n5  n6       n7 -  -  -  3

        where N is the lowest level node for which N.clade_perc >= pct_threshold. The most likely lowest level hit is
        the highest weight in the set containing all tips: {n5, n6, n3, n7}.

        Returns the node if it exceeds the provided local_threshold.

        """

        ## construct a list of all tip nodes in this clade
        tips = sorted([c for c in self.traverse() if len(c.children) == 0], key=lambda x: x.clade_perc, reverse=True)

        for tip in tips:
            tip.mother_clade = self

        if return_all_nodes:
            ## return the full ordered list of nodes
            return tips
        else:
            ## else return only the top hit as default

            ## checks if the node weight exceeds the provided threshold
            if tips[0].clade_perc > local_threshold:
                return tips[0]

            ## else return the base node
            else:
                self.mother_clade = self
                return self

    def traverse(self):
        """ Generator function for yielding all subnode of a given node
        """
        for child in self.children:
            for grandchild in child.traverse():
                yield grandchild
        else:
            yield self

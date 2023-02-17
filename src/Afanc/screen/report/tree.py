class Tree(object):
    'Tree node.'
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
        if self.ncbi_taxID in dbdict:
            acc = dbdict[self.ncbi_taxID][1]
        else: acc = "NA"

        if self.probability == None:
            json_data = {
                "name" : self.name,
                "reads" : self.taxon_reads,
                "percentage" : self.clade_perc,
                "taxon_id" : self.ncbi_taxID,
                "accession" : acc,
                }
        else:
            json_data = {
                "name" : self.name,
                "reads" : self.taxon_reads,
                "percentage" : self.clade_perc,
                "score" : self.probability,
                "taxon_id" : self.ncbi_taxID,
                "accession" : acc }

        if hasattr(self, "weighted_threshold"):
            json_data["local_threshold_calc"] = {
                "threshold_func" : self.threshold_func,
                "variables" : self.variables,
                "lower_bound_weighting" : self.lb_weight,
                "weighted_threshold" : self.weighted_threshold
                 }

        return json_data


    def find_local_max(self, variant_index, return_all_nodes=False):
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

        Returns the node if it exceeds the weighted lower bound read count.
        """

        ## construct a list of all tips which exceed the local threshold in this clade
        tips = []
        branch_mapped_reads = sum([c.taxon_reads for c in self.traverse()][:-1])

        for c in self.traverse():

            ## skip the root node
            if c == self:
                continue

            ## naive probability of this variant being correct
            ## simply a proportion of reads which map to this particular variant
            ## TODO: this block makes me nervous and needs some work
            if branch_mapped_reads != 0 and c.clade_reads <= branch_mapped_reads:
                c.probability = c.clade_reads/branch_mapped_reads

            # print(c.name)
            # print(variant_index[str(c.ncbi_taxID)])

            ## check to see if the parent and child ncbi taxonomy IDs exists within the variant index
            if str(c.parent.ncbi_taxID) in variant_index:
                """
                Where the parent taxon exist within the variant index, calculate the ratio between the
                normalised parent & child ANI distance as

                    R(c) = (sI(p) - sP(p))  /  (sI(c) - sP(c))
                           |______________|    |______________|
                                  |                   |
                              norm s(p)           norm s(c)
                where

                    s is some similarity function
                    sI(n) = the mean intrataxon similarity of taxon at node n
                    sP(n) = the mean parent-child similarity of taxon at node n
                    norm s(n) = the child-parent similarity of node n normalised against the mean intrataxon similarity of the parent taxon of n
                    c = child node
                    p = parent node of c

                The lower bound threshold weight is calculated as

                    w = 1 - R(c)

                The lower bound threshold for read count at a varian level taxon to call a hit is therefore

                    w * parent_taxon_reads
                """
                norm_s_p = variant_index[str(c.parent.ncbi_taxID)]["sibling_index"]["mean"] - variant_index[str(c.parent.ncbi_taxID)]["parent_index"]["mean"]  ## normalised parent ANI
                norm_s_c = variant_index[str(c.ncbi_taxID)]["sibling_index"]["mean"] - variant_index[str(c.ncbi_taxID)]["parent_index"]["mean"] ## normalised child ANI
                w = 1 - ((100*norm_s_p) * norm_s_c)

                c.threshold_func = f"w = 1 - ((100*norm_s_p) * norm_s_c)"
                c.variables = f"norm_s_p={norm_s_p} norm_s_c={norm_s_c}"

                # print("norm_s_c", norm_s_c)
                # print("norm_s_p", norm_s_p)
                # print(variant_index[str(c.ncbi_taxID)]["sibling_index"]["mean"], variant_index[str(c.ncbi_taxID)]["parent_index"]["mean"])

            elif str(c.parent.ncbi_taxID) not in variant_index and str(c.ncbi_taxID) in variant_index:
                """
                Where the parent taxon does not exist within the variant index but this node does (orphan node), it is assumed to be
                at a species level or higher (and therefore there exists no genome to calculate distance from).
                The lower bound threshold weight for orphan nodes is simply

                    (sI(c) - sP(c)) / 10

                where sP(n) is the mean parent-child similarity of taxon at node n.
                """
                norm_s_c = variant_index[str(c.ncbi_taxID)]["sibling_index"]["mean"] - variant_index[str(c.ncbi_taxID)]["parent_index"]["mean"] ## normalised child ANI
                w = norm_s_c/10

                c.threshold_func = f"w = norm_s_c/10"
                c.variables = f"norm_s_c={norm_s_c}"

                # print("norm_s_c", norm_s_c)
                # print(variant_index[str(c.ncbi_taxID)]["sibling_index"]["mean"], variant_index[str(c.ncbi_taxID)]["parent_index"]["mean"])

            elif str(c.parent.ncbi_taxID) not in variant_index and str(c.ncbi_taxID) not in variant_index:
                """
                Where neither the parent nor the child taxa exist within the variant index, the node is assumed to be genus level or higher.
                The threshold weighting is therefore 0.1.
                """
                w = 0.1

                c.threshold_func = f"""w = 0.1"""
                c.variable = None

            weighted_threshold = c.parent.clade_reads*w

            ## capture weighting values
            c.lb_weight = w
            c.weighted_threshold = weighted_threshold

            # print(f"w={w}\nparent reads={c.parent.clade_reads}\nTHRESHOLD : {weighted_threshold}\n")

            # print(c.name, c.level_int, weighted_threshold)

            if c.clade_reads >= weighted_threshold:
                tips.append(c)

        ## sort tips first by their taxonomic rank (by level_int) and then by the percentage of reads which map to that clade
        tips = sorted(tips, key=lambda x: f"{x.level_int}.{x.clade_perc}", reverse=True)

        for tip in tips:
            # print(tip.name, tip.clade_perc, tip.level_int)
            ## assign the mother clade for this node tip
            tip.mother_clade = self

        if return_all_nodes:
            ## return the full ordered list of nodes
            return tips
        else:
            ## else return only the top hit as default

            ## checks if there are any nodes which exceed the local threshold
            if len(tips) > 0:
                return tips[0]

            ## else return the base node
            else:
                self.mother_clade = self

                return self


    def traverse(self):
        """ Generator function for yielding all subnodes of a given node
        """
        for child in self.children:
            for grandchild in child.traverse():
                yield grandchild
        else:
            yield self

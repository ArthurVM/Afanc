from statistics import stdev
from collections import defaultdict
import pprint

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
                "lower_bound_weighting" : self.lb_weight,
                "weighted_threshold" : self.weighted_threshold
                 }

        return json_data


    def bayes_commute(self, variant_index, l, u):
        """ Redistribute reads from taxa which fall below the l*e to taxa which fall above it, where e is the elastic threshold.
        A key point of this function is that it handles each taxon independantly. Reads cannot be commuted up the tree, only across.
        For example, take the tree:

        0.            M.tb
                       |
              -------------------
              |        |        |
        1.  M.bovis  M.tb L1  M.tb L2
              |                 |
        2.  M.b BCG           M.tb L2.1

        Reads can be commuted between all taxa on level 1 (e.g from M.bovis to M.tb L1/L2) since they share the same parent, but not between taxa on level 2 (e.g. M.b BCG to M. tb L2.1)
        since they do not share the same parent, and are cousin taxa. However, reads which are redistributed from M.bovis to M.tb L2 are trickled into L2.1 in a number which conserves
        the proportion of reads from M.tb L2 which were assigned to subtaxa L2.1 prior to commuting. This means that monotypic tip level taxa which do not exceed the elastic threshold
        prior to commuting cannot exceed it after commuting, even if the parent does.
        """
        commute_dict = defaultdict(dict)    ## dictionary to collect commuting values

        ## construct a dictionary of form { parent.name :  [ child nodes ] } to ensure reads are not commuted across cousin taxa
        taxon_dict = defaultdict(list)

        for c in self.traverse():
            taxon_dict[c.parent.name].append(c)

        # for key, value in taxon_dict.items():
        #     print(key, [v.name for v in value])

        ## handle each taxon independantly
        for parent, children in taxon_dict.items():

            miss_tips = []  ## box for tips which fall below the lower bound of the commuting threshold
            hit_tips = []   ## box for tips which fall above the upper bound of the commuting threshold
            total_miss = 0    ## total number of reads from nodes within the miss_tips box
            total_hit = 0    ## total number of reads from the nodes within the hit_tips box

            for c in children:

                ## skip the root node
                if c == self:
                    continue

                ## calculate elastic threshold of this node
                c._elastic_threshold(variant_index)

                ## check to see if the node belongs to the symmetric set defined by the lower and upper bound commuting thresholds
                if c.clade_reads <= (c.weighted_threshold*l):
                    miss_tips.append(c)
                    total_miss += c.taxon_reads/c.weighted_threshold

                elif c.clade_reads >= (c.weighted_threshold*u):
                    hit_tips.append(c)
                    total_hit += c.taxon_reads/c.weighted_threshold

                ## capture the pre-commuting taxon read stats to calculate probabilities
                c.pc_taxon_reads = c.taxon_reads
                c.pc_clade_reads = c.clade_reads
                c.pc_clade_perc = c.clade_perc

            ## remove orphan nodes
            hit_tips = [i for i in hit_tips if i.parent in hit_tips or i.parent == self]

            t = 0
            rounded = 0

            for i in miss_tips:
                for j in hit_tips:
                    pos_prob = 1
                    j_prob = (j.pc_taxon_reads/j.weighted_threshold)/total_hit

                    bayes_prob = pos_prob * j_prob

                    commuted_reads = round(i.taxon_reads*bayes_prob)

                    commute_dict[j.name][i.name] = { "i_ID" : i.ncbi_taxID, "j_ID" : j.ncbi_taxID, "j_prob" : j_prob, "B_prob" : bayes_prob, "commuted_reads" : commuted_reads }

                    j.commuted_reads = commuted_reads

                    ## update read stats at this node
                    j.taxon_reads += commuted_reads
                    j.clade_reads += commuted_reads
                    j.clade_perc += (j.pc_clade_perc/j.pc_clade_reads)*commuted_reads

                    # for child in j.traverse():
                    #     print(child.name, child.clade_reads, child.clade_perc)
                    #     child.clade_reads = (child.clade_perc/100)*j.clade_reads
                    #     child.taxon_reads = (child.clade_perc/100)*j.clade_reads

                    # print(f"{j.commuted_reads} reads commuted from {i.name} to {j.name}. New total = {j.clade_reads} (B={bayes_prob}\t wt={j.weighted_threshold})")

            # print(f"{parent} HITS")
            # for c in hit_tips:
            #     print(c.name, c.clade_reads, c.weighted_threshold, c.variables, c.parent.clade_reads)
            #
            # print(f"{parent} MISSES")
            # for c in miss_tips:
            #     print(c.name, c.clade_reads, c.weighted_threshold, c.variables, c.parent.clade_reads)

        return commute_dict


    def find_local_max(self, variant_index, bayes_threshold=False, return_all_nodes=False):
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

        for c in self.traverse():

            ## skip the root node
            if c == self:
                continue

            c._elastic_threshold(variant_index)

            if c.clade_reads >= c.weighted_threshold:
                # print("HIT", c.name, "\t\t", c.clade_reads, c.weighted_threshold)
                tips.append(c)
            else:
                None
                # print(c.name, "\t\t", c.clade_reads, c.weighted_threshold)

        ## sort tips first by their taxonomic rank (by level_int) and then by the percentage of reads which map to that clade
        unfiltered_tips = sorted(tips, key=lambda x: f"{x.level_int}.{x.taxon_reads/x.weighted_threshold}", reverse=True)

        ## strip out any orphan tips
        ## i.e. tips where the parent did not pass the weighting threshold
        tips = [tip for tip in unfiltered_tips if tip.parent in tips or tip.parent == self]

        for tip in tips:
            ## assign the mother clade for this node tip
            tip.mother_clade = self
            # print(tip.name, "e=", tip.taxon_reads/tip.weighted_threshold)

        if return_all_nodes:
            return tips ## return the full ordered list of nodes

        else:
            ## else return only the top hit as default

            ## checks if there are any nodes which exceed the local threshold
            if len(tips) > 0:
                return tips[0]

            ## else return the base node
            else:
                self.mother_clade = self

                return self


    def _elastic_threshold(self, variant_index):
        ## check to see if the parent and child ncbi taxonomy IDs exists within the variant index
        if str(self.parent.ncbi_taxID) in variant_index:
            """
            Where the parent taxon exist within the variant index, calculate the ratio between the
            normalised parent & child ANI distance as

                R(c) = (fI(p) - fP(p))  /  (fI(c) - fP(c))
                       |______________|    |______________|
                              |                   |
                          norm f(p)           norm f(c)
            where

                f = some similarity function
                fI(n) = the mean intrataxon similarity of taxon at node n
                fP(n) = the mean parent-child similarity of taxon at node n
                norm f(n) = the child-parent similarity of node n normalised against the mean intrataxon similarity of the parent taxon of n
                c = child node
                p = parent node of c

            The lower bound threshold weight is calculated as

                w = 1 - R(c)

            The lower bound threshold for read count at a varian level taxon to call a hit is therefore

                w * parent_taxon_reads
            """
            # print("CONDITION 1", self.name, self.parent.name)

            norm_f_p =  variant_index[str(self.parent.ncbi_taxID)]["sibling_index"]["mean"] \
                        - variant_index[str(self.parent.ncbi_taxID)]["parent_index"]["mean"]  ## normalised parent ANI
            norm_f_c =  variant_index[str(self.ncbi_taxID)]["sibling_index"]["mean"]\
                        - variant_index[str(self.ncbi_taxID)]["parent_index"]["mean"] ## normalised child ANI
            w = 1 - ((100*norm_f_p) * norm_f_c)

            self.threshold_func = f"w = 1 - ((100*norm_f_p) * norm_f_c)"
            self.variables = f"norm_f_p={norm_f_p} norm_f_c={norm_f_c}"

            self.parent._elastic_threshold(variant_index)
            weighted_threshold = self.parent.weighted_threshold*w

            # print(variant_index[str(self.ncbi_taxID)]["sibling_index"]["mean"], variant_index[str(self.ncbi_taxID)]["parent_index"]["mean"])

        elif str(self.parent.ncbi_taxID) not in variant_index and str(self.ncbi_taxID) in variant_index:
            """
            Where the parent taxon does not exist within the variant index but this node does (orphan node), it is assumed to be
            at a species level or higher (and therefore there exists no genome to calculate distance from).
            The lower bound threshold weight for orphan nodes is simply

                (sI(c) - sP(c)) / 10

            where sP(n) is the mean parent-child similarity of taxon at node n.
            """
            # print("ORPHAN", self.name, self.parent.name)

            norm_f_c =  variant_index[str(self.ncbi_taxID)]["sibling_index"]["mean"] \
                        - variant_index[str(self.ncbi_taxID)]["parent_index"]["mean"] ## normalised child ANI
            w = norm_f_c/10

            self.threshold_func = f"w = norm_f_c/10"
            self.variables = f"norm_f_c={norm_f_c}"

            weighted_threshold = self.parent.clade_reads*w

            # print(variant_index[str(self.ncbi_taxID)]["sibling_index"]["mean"], variant_index[str(self.ncbi_taxID)]["parent_index"]["mean"])

        elif str(self.parent.ncbi_taxID) not in variant_index and str(self.ncbi_taxID) not in variant_index:
            """
            Where neither the parent nor the child taxa exist within the variant index, the node is assumed to be genus level or higher.
            The threshold weighting is therefore 0.1.
            """
            # print("CONDITION 3", self.name)

            w = 0.1
            self.threshold_func = f"""w = 0.1"""
            self.variables = None

            weighted_threshold = self.parent.clade_reads*w

        ## capture weighting values
        self.lb_weight = w
        self.weighted_threshold = weighted_threshold


    def traverse(self):
        """ Generator function for yielding all subnodes of a given node
        """
        for child in self.children:
            for grandchild in child.traverse():
                yield grandchild
        else:
            yield self

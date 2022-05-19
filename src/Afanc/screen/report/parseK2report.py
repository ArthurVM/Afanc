import json
from os import path
from collections import defaultdict

from Afanc.accdict import dbdict


class Taxa(object):
    """docstring for Taxa."""

    def __init__(self, line):
        super(Taxa, self).__init__()

        self.line = line

    def parseLine(self):
        """ parses a kraken2 report line
        """
        sline = self.line.strip("\n").split('\t')

        if len(sline) < 5:
            return []
        try:
            int(sline[1])
        except ValueError:
            return []

        #Extract relevant information
        perc = float(sline[0])
        clade_reads =  int(sline[1])
        taxa_reads = int(sline[2])
        taxa_level = sline[-3]
        ncbitaxID = sline[-2]

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

        self.name = name
        self.level_int = spaces
        self.clade_perc = perc
        self.clade_reads = clade_reads
        self.taxa_reads = taxa_reads
        self.taxa_level = taxa_level
        self.ncbi_taxID = ncbitaxID

    def makeJsonLine(self):
        """ generate a JSON output field
        """
        if self.ncbi_taxID in dbdict:
            acc = dbdict[self.ncbi_taxID][1]
        else: acc = "NA"

        return { "reads" : self.taxa_reads, "percentage" : self.clade_perc, "name" : self.name, "taxon_id" : self.ncbi_taxID, "accession" : acc }


def readK2report(report, pct_threshold, num_threshold):
    """ Read the kraken2 report and filter according to user defined values
    """

    resultsdict = defaultdict(list)

    with open(report, "r") as fin:
        for line in fin.readlines():

            taxa = Taxa(line)
            taxa.parseLine()

            if (taxa.clade_perc >= pct_threshold) & (taxa.clade_reads >= num_threshold):
                if taxa.level_int >= 7:
                    resultsdict[taxa.taxa_level].append(taxa)

    return resultsdict


def makeJson(resultsdict, pct_threshold, num_threshold, output_prefix, reportsDir):
    """ takes the results dict and generates a json report.
    {
    F : [taxa1, taxa2, ..., taxan],
    G : [taxa1, taxa2, ..., taxan],
    S : [taxa1, taxa2, ..., taxan],
    S1 : [taxa1, taxa2, ..., taxan],
    ...,
    Sn : [taxa1, taxa2, ..., taxan]
    }
    """

    taxa_level_key = { \
    "P" : "Phylum",
    "C" : "Class",
    "O" : "Order",
    "F" : "Family",
    "G" : "Genus",
    "G1": "Species Complex",
    "S" : "Species"
     }

    out_file = f"{reportsDir}/{output_prefix}.k2.json"

    json_dict = { "Thresholds" : { "reads" : num_threshold, "percentage" : pct_threshold } }

    for taxa_level, taxa_box in resultsdict.items():

        ## deal with taxa naming
        if ( taxa_level.startswith("S") ) & ( len(taxa_level) > 1 ):
            tl = f"Subspecies {taxa_level[1:]}"
        else:
            tl = taxa_level_key[taxa_level]

        json_dict[tl] = []

        for taxa in taxa_box:
            json_dict[tl].append(taxa.makeJsonLine())

    with open(out_file, "w") as fout:
        json.dump(json_dict, fout, indent = 4)


def parseK2reportMain(args):
    """ main function """
    report_path = f"{args.k2WDir}/{args.output_prefix}.k2.report.txt"
    resultsdict = readK2report(report_path, args.pct_threshold, args.num_threshold)
    makeJson(resultsdict, args.pct_threshold, args.num_threshold, args.output_prefix, args.reportsDir)

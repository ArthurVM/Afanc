import gzip
import json
from collections import defaultdict
from os import path
from Bio import SeqIO

from Afanc.accdict import dbdict


def parseK2output(k2out):
    """ parses a kraken2 output containing taxID information for
    each read.
    """
    readdict = defaultdict(str)
    with open(k2out, 'r') as fin:
        for line in fin.readlines():
            sline = line.split("\t")
            readID = sline[1]
            taxID = sline[2]

            readdict[readID] = (taxID)

    return readdict


def parseAccKeys(accKeys):
    """ parses an accesion keys tsv, with lines of form:
    taxID    GCA_key    Accession
    """
    taxIDdict = defaultdict(list)
    with open(accKeys, 'r') as fin:
        for line in fin.readlines():
            sline = line.strip("\n").split("\t")
            taxIDdict[sline[0]] = sline[1:]

    return taxIDdict


def parseJSONreport(json_report):
    """ parses kraken2 JSON report to define which taxIDs to
    extract from the readset.
    """
    taxIDs = []
    with open(json_report, 'r') as fin:
        for line in fin.readlines():
            if "taxon_id" in line:
                sline = line.split('"')
                taxIDs.append(sline[-2])

    return taxIDs


def splitFASTQ(fastq, dbdict, readdict, taxIDs):
    """ parses a fastq file and uses the readdict to split
    it into bins according to their kraken2 taxID designation.
    Outputs files for each:

    ./{GCA_key}_{Accession}_{1,2}.fq

    according to the taxIDdict.
    """
    readbins = defaultdict(list)

    with gzip.open(fastq, "rt") as fin:
        for read in SeqIO.parse(fin, 'fastq'):
            readID = read.id[:-2]
            if readID in readdict:
                taxID = readdict[readID]

                if taxID in taxIDs and taxID in dbdict:
                    outkey = "_".join(dbdict[taxID])
                    readbins[outkey].append(read)
                else:
                    readbins["unclassified"].append(read)
            else:
                readbins["unclassified"].append(read)

    for key, value in readbins.items():
        print(key, len(value))

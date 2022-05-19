import argparse
import sys
import time
import os
from numpy import mean, sort, cumsum


def genomeSize(fasta_file):
    """ takes a fasta file and sums the length of the sequences within
    """
    import gzip
    from Bio import SeqIO
    from numpy import sum

    if fasta_file.endswith(".gz"):
        fin = gzip.open(fasta_file, "rt")
    else:
        fin = open(fasta_file, "r")

    genlen = sum([len(rec.seq) for rec in SeqIO.parse(fin, "fasta")])
    fin.close()

    return genlen


def gini(depth_stdout, w=1):
    """Calculates the Gini coefficient of a coverage array isolated from a coverage .bed file
    Maps pretty well to A/(A+B) in its functionality, but effectively calculates Gini as the mean absolute difference.
    """
    lines = [int(line[2]) for line in depth_stdout if len(line) == 3]

    if w == 1:
        cov_array = lines
    else:
        cov_array = [mean(lines[i:i+w]) for i in range(0, len(lines), w)]

    s_cov_array=sort(cov_array)

    cov_sum = cumsum(s_cov_array)
    height_sum = cumsum(cov_sum)
    l_area_array = height_sum - cov_sum / 2

    eq_area = cov_sum[-1] * len(cov_array) / 2  # where eq_area is the area under the line of equality

    return (eq_area - l_area_array[-1]) / eq_area


def meanDOC(covarray):
    """ calculate the mean depth of coverage using a coverage array produced by samtools depth
    """
    return sum([int(i[2]) for i in covarray if len(i) == 3])/len(covarray)


def breadthofCoverage(covarray, genomesize):
    """ calculate the breadth of coverage across a genome
    """
    return len(covarray)/genomesize

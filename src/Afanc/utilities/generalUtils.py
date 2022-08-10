""" General utility funtions for Afanc
"""

import shutil
import re
import time
import sys
from os import path, listdir
from time import localtime, strftime


def isFile(filename):
    """ Checks if a path is an existing file """

    if not path.isfile(filename):
        vprint("MAIN", f"No file found at {filename}", "prRed", )
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))


def isDir(dirname):
    """ Checks if a path is an existing directory """

    if not path.isdir(dirname):
        vprint("MAIN", f"No file found at {dirname}", "prRed")
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(dirname)))


def checkDate(date):
    """ checks the ncbi date is of the correct form """

    sdate = date.split("-")

    if len(sdate) != 3 or len(sdate[0]) != 4 or sdate[1] != "05" or len(sdate[2]) != 2:
        vprint("MAIN", f"Date {date} is invalid. Please ensure the date is of the form YYYY-05-MM.")
        exit(3)
    else:
        return date


def parseBT2out(bt2out):
    """ parses the stderr from bowtie2 mapping
    """
    from collections import defaultdict

    bt2data = defaultdict(str)

    sbox = bt2out.split('\n')
    numreads = sbox[0].split(' ')[0]
    aln_con_0 = ' '.join(sbox[2].split(' ')[4:5])
    aln_con_1 = ' '.join(sbox[3].split(' ')[4:5])
    aln_con_gr1 = ' '.join(sbox[4].split(' ')[4:5])
    overallaln = sbox[-3].split('%')[0]     ## assumes there is a [bam_sort_core] line at the end of the output

    return overallaln


def gendbdict(cleanFasta_WDir):
    """ Generates a dict of form:
    {
        taxID : [ accession, assemblyID ],
        ...
    }
    """

    ## TODO: this is completely dependant on assembly filenames being of the form taxID_assemblyID_accession_genomic.fna
    ## this needs to be generalised
    from collections import defaultdict

    dbdict = defaultdict(list)
    for f in listdir(cleanFasta_WDir):
        fs = f.split("_")
        taxID = fs[0]
        assemblyID = "_".join(fs[1:3])
        accession = "_".join(fs[3:])
        dbdict[taxID] = [assemblyID, accession.split("_genomic")[0]]

    return dbdict


def reformat_mapping_arg(argument):
    """ takes the mapping argument and reformats it
    """
    return argument.replace("_", "-")


###################
# PRINT UTILITIES #
###################

colfuncs = {}
colfunc = lambda f: colfuncs.setdefault(f.__name__, f)

@colfunc
def prRed(sp):
    return f"\033[91m{sp}\033[00m"

@colfunc
def prGreen(sp):
    return f"\033[92m{sp}\033[00m"

@colfunc
def prYellow(sp):
    return f"\033[93m{sp}\033[00m"


def vprint(subprocess, info_text, colour, f=sys.stdout, end="\n"):
    """ controls process output
    """

    time = strftime("%H:%M:%S", localtime())

    print(
        f"\n{time} {colfuncs[colour](subprocess)} :: {info_text}",
        end=end,
        file=f,
        flush=True,
    )

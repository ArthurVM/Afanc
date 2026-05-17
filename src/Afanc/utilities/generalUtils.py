""" General utility funtions for Afanc
"""

import sys
from os import path, listdir
from pathlib import Path
from time import localtime, strftime

from .exceptions import FileNotFoundErrorAfanc, DirectoryNotFoundErrorAfanc, InvalidFileFormatError


def isFile(filename):
    """ Checks if a path is an existing file """

    if not path.isfile(filename):
        # vprint("MAIN", f"No file found at {filename}", "prRed")
        raise FileNotFoundErrorAfanc(filename)
    else:
        return Path(path.abspath(path.realpath(path.expanduser(filename))))


def isDir(dirname):
    """ Checks if a path is an existing directory """

    if not path.isdir(dirname):
        # vprint("MAIN", f"No directory found at {dirname}", "prRed")
        raise DirectoryNotFoundErrorAfanc(dirname)
    else:
        return Path(path.abspath(path.realpath(path.expanduser(dirname))))
    

def checkNcbiTaxDB(ncbi_taxdb):
    """ Checks that the ncbi taxonomy database is correctly formed """
    ncbi_taxdb_abs = isDir(ncbi_taxdb)
    required_files = ["nodes.dmp", "names.dmp", "merged.dmp"]

    for f in required_files:
        isFile(path.join(ncbi_taxdb, f))
    
    return ncbi_taxdb_abs


def checkDate(date):
    """ checks the ncbi date is of the correct form """

    sdate = date.split("-")

    if len(sdate) != 3 or len(sdate[0]) != 4 or sdate[1] != "05" or len(sdate[2]) != 2:
        # vprint("MAIN", f"Date {date} is invalid. Please ensure the date is of the form YYYY-05-MM.")
        # exit(3)
        raise InvalidFileFormatError(file_path=date, details="Date is invalid. Please ensure the date is of the form YYYY-05-MM (with 05 as the month).")

    else:
        return date


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


def iupac(n):
    iupac_designation = {
        "A":["A"],
        "C":["C"],
        "G":["G"],
        "T":["T"],
        "R":["A","G"],
        "Y":["C","T"],
        "S":["G","C"],
        "W":["A","T"],
        "K":["G","T"],
        "M":["A","C"],
        "B":["C","G","T"],
        "D":["A","G","T"],
        "H":["A","C","T"],
        "V":["A","C","G"],
        "N":["A","C","G","T"]
    }
    return iupac_designation[n]


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

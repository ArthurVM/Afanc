""" makes the Afanc working directory and directory structure
"""

import sys
import random, string
from os import mkdir, path, listdir, chdir

from Afanc.utilities.generalUtils import isDir, vprint

subprocessID = "MAIN"


def genKey():
    return "".join(random.choices(string.ascii_letters + string.digits, k=6))


def genScreenDirStructure(args):
    """ Generates the file structure to deposit Afanc screen run files.
    This file structure is as follows:

         working_directory
                 |
            Afanc.{runkey}
                 |
        -------------------
        |        |        |
     reports    k2WDir    bt2WD

    Generates files with a randomly generated run key, which is (probably) unique to this run.

    """
    ## generates a unique run code to name directories, a bit hacky...

    if args.run_key:
        ## generate the run directory using a randomly generated run key
        key = genKey()
        vprint(subprocessID, f"Run key for this run: {key}\n", "prYellow")

        args.runWDir = path.join(isDir("./"), f"{args.output_prefix}.{key}")

    else:
        args.runWDir = path.join(isDir("./"), f"{args.output_prefix}")

    vprint(
        subprocessID,
        f"Generating directory structure for this Afanc-screen run in working directory {args.runWDir}\n",
        "prYellow",
    )

    if not path.exists(args.runWDir):
        mkdir(args.runWDir)
    else:
        vprint(
            subprocessID,
            f"{args.runWDir} exists! Aborting to prevent accidental data loss...\n",
            "prRed",
            sys.stderr,
        )
        exit(3)

    args.reportsDir = path.join(args.runWDir, "reports")
    mkdir(args.reportsDir)  ## make a directory to dump reports

    args.k2WDir = path.join(args.runWDir, "k2WDir")
    mkdir(args.k2WDir)  ## make directory to dump misc output files

    args.bt2WDir = path.join(args.runWDir, "bt2WD")
    mkdir(args.bt2WDir)  ## make directory to contain data for running bowtie2 mapping


def initAutoDBDirStructure(args):
    """ Initialise (but does not create) the file structure to deposit Autodatabase run files.
    """
    args.autoDB_WDir = path.abspath(f"./{args.output_prefix}")
    args.fasta_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_Fasta")
    args.cleanFasta_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_cleanFasta")
    args.mash_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_mash")
    args.qc_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_qc")
    args.kraken2_WDir = path.abspath(f"{args.autoDB_WDir}/krakenBuild_autoDatabase_kraken2Build")
    args.krona_WDir = path.abspath(f"{args.autoDB_WDir}/krakenBuild_autoDatabase_krona")


def checkautodbWD(args):
    """ if fetch_assemblies is False, check the autodatabase directory for required files
    """
    wd_box = listdir(args.database)

    if "krakenBuild_autoDatabase_kraken2Build" in wd_box:
        args.k2_database = path.join(args.database, "krakenBuild_autoDatabase_kraken2Build")
    else:
        vprint(
            subprocessID,
            f"Kraken2 database not found in the autodatabase results directory! Exiting...",
            "prRed",
            sys.stderr,
        )
        exit(4)
    if "selectFasta_autoDatabase_cleanFasta" in wd_box:
        args.cleanFasta = path.join(args.database, "selectFasta_autoDatabase_cleanFasta/")
    else:
        vprint(
            subprocessID,
            f"Fasta file database not found in the autodatabase results directory! Exiting...",
            "prRed",
            sys.stderr,
        )
        exit(4)


def mkchdir(dir, ch=True):
    """ make and change to a dir
    """
    if not path.isdir(dir): mkdir(dir)
    if ch: chdir(dir)

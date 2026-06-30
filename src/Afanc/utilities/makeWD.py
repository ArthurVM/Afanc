""" makes the Afanc working directory and directory structure
"""

import sys
import random, string
from os import mkdir, path, listdir, chdir

from .generalUtils import isDir, vprint

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
     reports    k2WDir    mappingWD

    Generates files with a randomly generated run key, which is (probably) unique to this run.

    """

    args.baseRunDir = path.abspath(isDir("./"))

    ## generates a unique run code to name directories
    if args.run_key:
        ## generate the run directory using a randomly generated run key
        key = genKey()
        vprint(subprocessID, f"Run key for this run: {key}\n", "prYellow")

        args.runWDir = path.join(args.baseRunDir, f"{args.output_prefix}.{key}")

    else:
        args.runWDir = path.join(args.baseRunDir, f"{args.output_prefix}")

    vprint(
        subprocessID,
        f"Generating directory structure for this Afanc-screen run in working directory {args.runWDir}",
        "prYellow",
    )

    if not path.exists(args.runWDir):
        mkdir(args.runWDir)
    else:
        vprint(
            subprocessID,
            f"{args.runWDir} exists! Aborting to prevent accidental data loss...",
            "prRed",
            sys.stderr,
        )
        exit(3)

    args.reportsDir = path.join(args.runWDir, "reports")
    mkdir(args.reportsDir)  ## make a directory to dump reports

    args.k2WDir = path.join(args.runWDir, "k2WDir")
    mkdir(args.k2WDir)  ## make directory to dump misc output files

    args.mappingWDir = path.join(args.runWDir, "mappingWD")
    mkdir(args.mappingWDir)  ## make directory to contain data for read mapping

    args.profilerWDir = path.join(args.runWDir, "profilerWD")
    mkdir(args.profilerWDir)  ## make directory to contain SNP and lineage profiling data

    args.kronaWDir = path.join(args.runWDir, "kronaWD")
    mkdir(args.kronaWDir)  ## make directory to contain Krona reports


def initAutoDBDirStructure(args):
    """ Initialise (but does not create) the file structure to deposit Autodatabase run files.
    """
    args.autoDB_WDir = path.abspath(f"./{args.output_prefix}")
    args.fasta_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_Fasta")
    args.mash_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_mash")
    args.qc_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_qc")
    args.cleanFasta_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_cleanFasta")
    args.variant_index_WDir = path.abspath(f"{args.autoDB_WDir}/selectFasta_autoDatabase_variantIndex")
    args.kraken2_WDir = path.abspath(f"{args.autoDB_WDir}/krakenBuild_autoDatabase_kraken2Build")
    args.krona_WDir = path.abspath(f"{args.autoDB_WDir}/krakenBuild_autoDatabase_krona")


def checkautodbWD(args):
    """ if fetch_assemblies is False, check the autodatabase directory for required files
    """
    from Afanc.screen.profiles import load_profiles_manifest

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

    if "assemblies.json" in wd_box:
        args.assemblies_json = path.join(args.database, "assemblies.json")
    else:
        vprint(
            subprocessID,
            f"assemblies.json not found in the autodatabase results directory! Exiting...",
            "prRed",
            sys.stderr,
        )
        exit(4)

    if "fastas_in_DB.json" in wd_box:
        args.db_fastas_json = path.join(args.database, "fastas_in_DB.json")
    else:
        vprint(
            subprocessID,
            f"fastas_in_DB.json not found in the autodatabase results directory! Exiting...",
            "prRed",
            sys.stderr,
        )
        exit(4)

    mvi_dir = path.join(args.database, "selectFasta_autoDatabase_variantIndex")

    if not path.isdir(mvi_dir):
        vprint(
            subprocessID,
            f"Mash Variant Index directory not found in the autodatabase results directory: {mvi_dir}. Exiting...",
            "prRed",
            sys.stderr,
        )
        exit(4)

    mvi_box = listdir(mvi_dir)
    mvi_files = [
        "variant_index.json",
        "mash_all_vs_all.tsv",
        "mash_genus_matrix.tsv",
        "mash_species_matrix.tsv",
    ]
    missing_mvi_files = [mvi_file for mvi_file in mvi_files if mvi_file not in mvi_box]

    if len(missing_mvi_files) == 0:
        args.variant_index_path = path.join(mvi_dir, "variant_index.json")
        args.mvi_path = mvi_dir
    else:
        vprint(
            subprocessID,
            f"Mash Variant Index files not found in {mvi_dir}: {', '.join(missing_mvi_files)}. Exiting...",
            "prRed",
            sys.stderr,
        )
        exit(4)

    args.profiles_dir, args.profiles_manifest = load_profiles_manifest(
        args.database,
        profiles_dir=getattr(args, "profiles_dir", None),
    )
    args.lineage_profiles_by_accession = {}

def mkchdir(dir, ch=True):
    """ make and change to a dir
    """
    if not path.isdir(dir): mkdir(dir)
    if ch: chdir(dir)

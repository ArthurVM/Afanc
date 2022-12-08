"""
cat mykrobeDB_species.txt | while read -r acc ; do esearch -db assembly -query 'Mycobacterium_abscessus_subsp_bolletii[organism] AND latest[filter]' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | head -n 1 | while read -r url ; do
fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
wget "$url/$fname" ;
done ; done ;
"""

import os
import sys
import subprocess
import shutil
import json
from collections import defaultdict
from os import path

from Afanc.utilities.runCommands import command


def getAccessions(report, dbdict):

    ## TODO: fix this

    accessions = []

    with open(report, 'r') as fin:
        json_data = json.load(fin)

        for subdict in json_data["Detection_events"]:
            accessions.append(dbdict[str(subdict["taxon_id"])][0])

    return accessions


def get_hitIDs(out_json):
    """ get hit IDS from the kraken2 JSON report file
    """

    with open(out_json, 'r') as fin:
        json_data = json.load(fin)

        hit_ids = [ subdict["closest_variant"]["name"] if "closest_variant" in subdict else subdict["name"] for subdict in json_data["Detection_events"] ]

    return hit_ids


def getGenomesbyAcc(genomes, args):
    """ Get genomes using accessions
    """

    stdoutstr = ""
    stderrstr = ""

    for g in genomes:
        download_genome(g, args.stdout, args.stderr)


def getGenomesbyName(assembly_ids, args):
    """ Get genome assemblies for each hit from genbank using the Entrez suite.
    Takes assembly names from the kraken2 report json.
    """

    stdoutstr = ""
    stderrstr = ""

    for g in assembly_ids:
        download_genome(g, args.stdout, args.stderr)


def download_genome(assembly, stdout, stderr, taxID=False):
    """ Download a genome from genbank using the Entrez suite.
    If taxID is provided, this is appended to the front of the assembly filename
    followed by a trailing underscore.
    """

    runline = f"esearch -db assembly -query \'{assembly}[organism] AND latest[filter]\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
    tmp_ftp_dir = command(runline, "DOWNLOAD_HITS").run_comm_quiet(1, stdout, stderr)

    ## get all ftp directories belonging to the given assembly
    ftp_dirs = [ f for f in tmp_ftp_dir[0].decode().split("\n") if f != '']

    ## check anything was found
    if len(ftp_dirs) == 0:
        print(f"{assembly} NOT FOUND ON GENBANK! SKIPPING...")
        return 2

    ## get first ftp directory
    ftp_dir = ftp_dirs[0]

    base = os.path.basename(ftp_dir.split("/")[-1])

    ## add taxID to the start of the output assembly if provided
    if not taxID:
        outfile = f"{base}_genomic.fna.gz"
    else:
        outfile = f"{taxID}_{base}_genomic.fna.gz"

    if os.path.exists(outfile):
        print(f"FILE ALREADY EXISTS! SKIPPING...")
        return 3
    else:
        grepline = f"curl {ftp_dir}/{base}_genomic.fna.gz --output {outfile}"
        stdout, stderr = command(grepline, "DOWNLOAD_HITS").run_comm(1, stdout, stderr)


def getLocalGenomes(out_json, args):
    """ if args.fetch_assemblies is False, get genomes from the autodatabase results directory:
            autodb_results/selectFasta_autoDatabase_cleanFasta/
    """
    from Afanc.utilities.runCommands import command

    with open(args.assemblies_json, "r") as fin:
        assembly_dict = json.load(fin)

    ## construct a dictionary from taxID : local_assemblyID pairs
    db_pathdict = { c.split("_")[0] : c for c in os.listdir(args.cleanFasta) }
    missing_assemblies = []

    ## read in kraken2 report json to get hit genomes for mapping to
    with open(out_json, 'r') as fin:
        json_data = json.load(fin)

        for subdict in json_data["Detection_events"]:

            ## check to see if there is variant information
            ## if so, make this the target
            if "closest_variant" in subdict:
                taxID = str(subdict["closest_variant"]["taxon_id"])
                assemblyID = subdict["closest_variant"]["name"]

            else:
                taxID = str(subdict["taxon_id"])
                assemblyID = subdict["name"]

            ## if the local assembly cannot be found using the ncbi taxID, download genome from GenBank
            ##
            ## TODO : this is the result of a bug, where the initial directory structure is incorrect,
            ## leading to incorrect designation of an ncbi taxID as a file name prefix, but correct
            ## designation by kraken2. This leads to inconsistent taxID handling. Downloading from
            ## GenBank as shown here is a non-optimal solution, but works.
            if taxID not in db_pathdict:
                print(f"Cannot find local assembly for {assemblyID} with ncbi taxID {taxID}! Now attempting to download from GenBank...")
                ## collect missing assemblies for unknown purpose
                missing_assemblies.append(assemblyID)

                download_genome(assemblyID, args.stdout, args.stderr, taxID)

            ## if the assembly can be found, copy to the bt2 working directory for bt2 db construction
            else:
                assembly_path = path.join(args.cleanFasta, db_pathdict[taxID])
                print(f"Copying {assembly_path}...")
                shutil.copy(assembly_path, "./")

"""
cat mykrobeDB_species.txt | while read -r acc ; do esearch -db assembly -query 'Mycobacterium_abscessus_subsp_bolletii[organism] AND latest[filter]' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | head -n 1 | while read -r url ; do
fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
wget "$url/$fname" ;
done ; done ;
"""

import os
import sys
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

        return None

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

        return outfile

    else:
        grepline = f"curl {ftp_dir}/{base}_genomic.fna.gz --output {outfile}"
        stdout, stderr = command(grepline, "DOWNLOAD_HITS").run_comm(1, stdout, stderr)

        return outfile


def check_hits(args, json_data):
    """ Checks to see if the species hits called should be profiled for variants.
    """

    species_box =  [ f["name"].replace(" ", "_") for f in json_data["Detection_events"] ]

    if args.variant_profile != False:
        variant_species = {}
        for species in species_box:
            if species in args.variant_profile:
                variant_species[species] = args.variant_profile[species][0]

        if len(variant_species) != 0:
            return variant_species

    return {False : False}


def getLocalGenomes(out_json, args):
    """ if args.fetch_assemblies is False, get genomes from the autodatabase results directory:
            autodb_results/selectFasta_autoDatabase_cleanFasta/
    """
    from Afanc.utilities.runCommands import command

    with open(args.db_fastas_json, "r") as fin:
        assembly_dict = json.load(fin)

    missing_assemblies = []

    ## read in species report json to get hit genomes for mapping to
    with open(out_json, 'r') as fin:
        json_data = json.load(fin)

        variant_species = check_hits(args, json_data)

        for subdict in json_data["Detection_events"]:

            species_ID = subdict["name"].replace(" ", "_")

            ## check if the hit species needs to be passed to the variant profiling module
            ## if so, return the reference fasta captured within the variants tsv as the fasta to map to
            if species_ID in variant_species:
                assembly_path = variant_species[species_ID]
                print(f"Copying {assembly_path}...")
                shutil.copy(assembly_path, "./")
                subdict["assembly"] = assembly_path

                continue

            ## check to see if there is clustering variant information
            ## if so, make this the target
            elif "closest_variant" in subdict and species_ID not in variant_species:
                subdict = subdict["closest_variant"]
                taxID = str(subdict["taxon_id"])
                assemblyID = subdict["name"]

            else:
                taxID = str(subdict["taxon_id"])
                assemblyID = subdict["name"]

            species_ID = assemblyID.replace(" ", "_")

            ## copy assembly to the bt2 working directory for bt2 db construction
            if taxID in assembly_dict:
                assembly_path = assembly_dict[taxID][0]
            elif species_ID in assembly_dict:
                assembly_path = assembly_dict[species_ID][0]
            else:
                print(taxID, species_ID)

            print(f"Copying {assembly_path}...")
            shutil.copy(assembly_path, "./")

            if not assembly_path == None:
                subdict["assembly"] = assembly_path.split("/")[-1]
            else:
                subdict["assembly"] = assembly_path

    with open(out_json, "w") as fout:
        json.dump(json_data, fout, indent=4)

    return variant_species

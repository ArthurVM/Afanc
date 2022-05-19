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
from collections import defaultdict
from os import path

from Afanc.utilities.runCommands import command


def getNames(report):
    """ parses the kraken2 report and extracts species/variant names to download
    """
    None


def getAccessions(report, dbdict):

    acessions = []

    with open(report, 'r') as fin:
        for line in fin.readlines():
            if "taxon_id" in line:
                ncbitaxID = line.split('\"')[-2]
                if ncbitaxID in dbdict:
                    acessions.append(dbdict[ncbitaxID][1])

    return acessions


def getGenomesbyAcc(genomes):

    stdoutstr = ""
    stderrstr = ""

    for g in genomes:
        runline = f"esearch -db assembly -query \'{g}\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | head -1"
        tmp_ftp_dir = command(runline, "DOWNLOAD_HITS").run_comm_quiet(1, args.stdout, args.stderr)

        ftp_dir = tmp_ftp_dir.decode().strip("\n")

        base = os.path.basename(ftp_dir.split("/")[-1])

        if os.path.exists(f"./{base}_genomic.fna"):
            print(f"FILE ALREADY EXISTS! SKIPPING...")
            continue

        grepline = f"curl {ftp_dir}/{base}_genomic.fna.gz --output {base}_genomic.fna.gz"
        stdout, stderr = command(grepline, "DOWNLOAD_HITS").run_comm(1, args.stdout, args.stderr)

        stdoutstr += stdout.decode() + "\n"
        stderrstr += stderr.decode() + "\n"

    return stdoutstr, stderrstr


def getGenomesbyName(accessions):

    stdoutstr = ""
    stderrstr = ""

    for g in accessions:
        runline = f"esearch -db assembly -query \'{g}[organism] AND latest[filter]\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
        ftp_dir = command(runline, "DOWNLOAD_HITS").run_comm_quiet(1, args.stdout, args.stderr)

        ftp_dir = ftp_dir.decode().strip("\n")

        base = os.path.basename(ftp_dir)

        if os.path.exists(f"./{base}_genomic.fna"):
            print(f"FILE ALREADY EXISTS! SKIPPING...")
            continue

        grepline = f"curl {ftp_dir}/{base}_genomic.fna.gz --output {base}_genomic.fna.gz"
        stdout, stderr = command(grepline, "DOWNLOAD_HITS").run_comm(1, args.stdout, args.stderr)

        stdoutstr += stdout.decode() + "\n"
        stderrstr += stderr.decode() + "\n"

    return stdoutstr, stderrstr


def getLocalGenomes(accessions, autoDB_fasta_dir):
    """ if args.fetch_assemblies is False, get genomes from the autodatabase results directory:
            autodb_results/selectFasta_autoDatabase_cleanFasta/
    """
    from Afanc.utilities.runCommands import command

    db_genomes = os.listdir(autoDB_fasta_dir)
    for g in accessions:
        for dbg in db_genomes:
            if g in dbg:
                assembly_path = path.join(autoDB_fasta_dir, dbg)
                print(f"Copying {assembly_path}...")
                shutil.copy(assembly_path, "./")

    return "", ""

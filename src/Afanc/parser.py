"""
Parse arguments for Afanc
"""

import sys
import argparse

from Afanc.utilities.generalUtils import isDir, isFile, checkDate


def run_subtool(args):
    ## level 0 run function

    ## initialise log files to deposit stdout and stderr when necessary
    initLogFiles(args)

    if args.command == "autodatabase":
        ## run autodatabase from a fasta directory structure
        from Afanc.autodatabase.runFuncs import runAutoDB
        runAutoDB(args)

    elif args.command == "screen":
        ## screen fastq using a database
        from Afanc.screen.runFuncs import runScreen
        runScreen(args)


def initLogFiles(args):
    """ log files within the args object to dump stderr and stdout
    """
    args.stdout = open(f"{args.output_prefix}.stdout.txt", "a")
    args.stderr = open(f"{args.output_prefix}.stderr.txt", "a")


base_parser = argparse.ArgumentParser(add_help=True)

subparsers = base_parser.add_subparsers(
    title="[sub-commands]", dest="command"
)

## autodatabase args parser
parser_autodb = subparsers.add_parser(
    "autodatabase",
    help="Generate a database from a FASTA directory structure with autodatabase.",
)

parser_autodb.add_argument('fastaDir',
    type=isDir,
    action='store',
    help='Path to the species delimited directory structure used to construct the database.')

parser_autodb.add_argument('-o', '--output_prefix',
    type=str,
    default="Afanc_autodb",
    action='store',
    help='Output prefix for this run. Default=Afanc_autodb.')

parser_autodb.add_argument('-d', '--ncbi_date',
    type=checkDate,
    default="2020-05-01",
    action='store',
    help='The date of NCBI taxonomy to download. Must be of the form YYYY-05-MM. Default=2020-05-01.')

parser_autodb.add_argument('-m', '--mode_range',
    type=float,
    default=0.1,
    action='store',
    help='Range to take around the mode of the average mash distance. Default=0.1.')

parser_autodb.set_defaults(func=run_subtool)

## Afanc-screen parser
parser_screen = subparsers.add_parser(
    "screen",
    help="High-resolution metagenomic screening of short read data using a database constructed by autodatabase.",
)

parser_screen.add_argument('database',
    type=isDir,
    action='store',
    help='Path to the results directory created by running autodatabase. If fetch_assemblies is provided, this should just be the path to a Kraken2 database.')

parser_screen.add_argument('fastq',
    type=isFile,
    nargs=2,
    action='store',
    help='Fastq files to screen.')

parser_screen.add_argument('-p', '--pct_threshold',
    type=float,
    default=0.2,
    action='store',
    help='Min. coverage, as %%. Default=1.0.')

parser_screen.add_argument('-n', '--num_threshold',
    type=int,
    default=1000,
    action='store',
    help='Min. coverage over the clade to score a hit, as no. of reads. Should be a positive integer. Default=1000.')

parser_screen.add_argument('-o', '--output_prefix',
    type=str,
    default="Afanc_screen",
    action='store',
    help='Output prefix for this run. Default=results.')

parser_screen.add_argument('-k', '--run_key',
    default=False,
    action='store_true',
    help='Generate and use a random key for this runs run directory. Results will be stored in {output_prefix}.{random_key}')

parser_screen.add_argument('-f', '--fetch_assemblies',
    default=False,
    action='store_true',
    help='Fetch genome assemblies for species hits using from GenBank using the ENSEMBL software suite. Default=False; get assemblies from the autodb working directory.')

parser_screen.add_argument('-t', '--threads',
    type=int,
    default=4,
    action='store',
    help='Number of threads to used for this run. Default=4.')

parser_screen.set_defaults(func=run_subtool)

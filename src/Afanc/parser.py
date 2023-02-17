"""
Afanc high-resolution Metagenomics disambiguator.
"""
import sys
import argparse
from ._version import __version__

from Afanc.utilities.generalUtils import isDir, isFile, checkDate, reformat_mapping_arg

"""
Parse arguments for Afanc
"""

def run_subtool(args):
    ## level 0 run function

    ## initialise log files to deposit stdout and stderr when necessary
    initLogFiles(args)

    if args.command == "get_dataset":
        ## get a dataset from genbank and generate a directory structure which can be used by autodatabase
        from Afanc.get_dataset.download_assemblies import runGet_dataset
        runGet_dataset(args)

    elif args.command == "autodatabase":
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


base_parser = argparse.ArgumentParser(add_help=True, description=__doc__)

base_parser.add_argument('-v', '--version', action='version',
    version='%(prog)s {version}'.format(version=__version__))

subparsers = base_parser.add_subparsers(
    title="[sub-commands]", dest="command"
)

## get_dataset args parser
parser_getDataset = subparsers.add_parser(
    "get_dataset",
    help="Download a dataset of genome assemblies from GenBank.",
)

parser_getDataset.add_argument('ID_file',
    type=isFile,
    action='store',
    help='List of line seperated IDs to download. By default, these will be assumed to be species names. If the -a flag is provided, they will be assumed as accession IDs.')

parser_getDataset.add_argument('-a', '--accessions',
    action='store_true',
    default=False,
    help='Flag to specify that ID_file contains accession IDs, rather than species names.')

parser_getDataset.add_argument('-n', '--num_assemblies',
    type=int,
    default=1,
    action='store',
    help='The number of assemblies for each species to download. If the -a flag is provided, this is ignored. Default=1.')

parser_getDataset.add_argument('-o', '--output_prefix',
    type=str,
    default='assemblies',
    action='store',
    help='Output prefix for this run. Default=assemblies.')


parser_getDataset.set_defaults(func=run_subtool)

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

parser_autodb.add_argument('-n', '--ncbi_date',
    type=checkDate,
    default="2022-05-01",
    action='store',
    help='The date of NCBI taxonomy to download. Must be of the form YYYY-05-MM. Default=2022-05-01.')

parser_autodb.add_argument('-m', '--mode_range',
    type=float,
    default=0.1,
    action='store',
    help='Range to take around the mode of the average mash distance. Default=0.1.')

parser_autodb.add_argument('-t', '--threads',
    type=int,
    default=4,
    action='store',
    help='Number of threads to used for this run. Default=4.')

parser_autodb.set_defaults(func=run_subtool)

## Afanc-screen parser
parser_screen = subparsers.add_parser(
    "screen",
    help="High-resolution metagenomic screening of short read data using a database constructed by autodatabase.",
)

parser_screen.add_argument('database',
    type=isDir,
    action='store',
    help='Path to the results directory created by running autodatabase.')

parser_screen.add_argument('fastq',
    type=isFile,
    nargs=2,
    action='store',
    help='Fastq files to screen.')

parser_screen.add_argument('-p', '--pct_threshold',
    type=float,
    default=5.0,
    action='store',
    help='Min. %% of reads to call a global hit. Default=5.0.')

parser_screen.add_argument('-n', '--num_threshold',
    type=int,
    default=1000,
    action='store',
    help='Min. coverage over the clade to score a hit, as no. of reads. Should be a positive integer. Default=1000.')

parser_screen.add_argument('-m', '--mapping_sensitivity',
    type=reformat_mapping_arg,
    default="very-sensitive",
    action='store',
    choices=["very-sensitive", "sensitive", "fast", "very-fast"],
    help='Sensitivity of mapping to use when mapping reads to suspected target genomes.')

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
    help='Fetch genome assemblies for species hits using from GenBank using the ENSEMBL software suite. Default=False; get assemblies from the autodb working directory where possible.')

parser_screen.add_argument('-t', '--threads',
    type=int,
    default=4,
    action='store',
    help='Number of threads to used for this run. Default=4.')

clean_group = parser_screen.add_mutually_exclusive_group()

clean_group.add_argument('-c', '--clean',
    default=False,
    action='store_true',
    help='Remove the bowtie2 working directory from the output directory. Default=False.')

clean_group.add_argument('-s', '--superclean',
    default=False,
    action='store_true',
    help='Delete the entire output directory, leaving only log files and the results .json. Default=False.')


parser_screen.set_defaults(func=run_subtool)

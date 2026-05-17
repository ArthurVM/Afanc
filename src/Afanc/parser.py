"""
Afanc high-resolution Metagenomics disambiguator.
"""
import sys
import argparse
from ._version import __version__

from Afanc.utilities.generalUtils import isDir, isFile, checkNcbiTaxDB, checkDate

"""
Parse arguments for Afanc
"""

def runSubtool(args):
    ## level 0 run function

    ## initialise log files to deposit stdout and stderr when necessary
    initLogFiles(args)

    if args.command == "get_dataset":
        ## get a dataset from genbank and generate a directory structure which can be used by autodatabase
        from Afanc.get_dataset.download_assemblies import runGet_dataset

        splash = f"""
        RUNNING GET_DATASET
        ===================

        """
        print(splash, file=args.stdout)
        print(splash, file=args.stderr)

        runGet_dataset(args)

    elif args.command == "autodatabase":
        ## run autodatabase from a fasta directory structure
        from Afanc.autodatabase.runFuncs import runAutoDB
        splash = f"""
        RUNNING AUTODATABASE
        ====================

        """
        print(splash, file=args.stdout)
        print(splash, file=args.stderr)

        runAutoDB(args)

    elif args.command == "screen":
        ## screen fastq using a database
        from Afanc.screen.runFuncs import runScreen
        splash = f"""
        RUNNING SCREEN
        ==============

        """
        print(splash, file=args.stdout)
        print(splash, file=args.stderr)

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


parser_getDataset.set_defaults(func=runSubtool)

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
    default="2026-05-01",
    action='store',
    help='The date of NCBI taxonomy to download. Must be of the form YYYY-05-MM. Default=2026-05-01.')

parser_autodb.add_argument('-m', '--mode_range',
    type=float,
    default=0.1,
    action='store',
    help='Range to take around the mode of the average mash distance. Default=0.1.')

parser_autodb.add_argument('-v', '--variant_index_method',
    choices=['mash'],
    default='mash',
    action='store',
    help='Method for calculating intra/inter-taxon distances for the Mash Variant Index. Default=mash.')

parser_autodb.add_argument('-x', '--ncbi_tax_db',
    type=checkNcbiTaxDB,
    default=False,
    action='store',
    help='Use a locally stored ncbi taxonomy database instead of downloading from ncbi. Default=False.')

parser_autodb.add_argument('-f', "--use_ftp",
    action="store_true",
    default=False,
    help="Flag. Use ftp instead of https/rsync when downloading, Default=False.")

# parser_autodb.add_argument('-d', '--stdev',
#     type=float,
#     default=1.0,
#     action='store',
#     help='Number of standard deviations around the mean mash distance for taxon sequence selection. Default=1.0.')

parser_autodb.add_argument('-t', '--threads',
    type=int,
    default=4,
    action='store',
    help='Number of threads to used for this run. Default=4.')

clean_group = parser_autodb.add_mutually_exclusive_group()

clean_group.add_argument('-c', '--clean',
    default=False,
    action='store_true',
    help='Remove the mash working directory from the output directory. Default=False.')

clean_group.add_argument('-s', '--superclean',
    default=False,
    action='store_true',
    help='Remove all files not required for running Afanc screen from the output directory. Default=False.')

parser_autodb.set_defaults(func=runSubtool)

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

parser_screen.add_argument('-u', '--upper_bound',
    type=float,
    default=0.75,
    action='store',
    help='The upper bound fraction of the elastic threshold to designate read commute recipients. Default=0.75.')

parser_screen.add_argument('-l', '--lower_bound',
    type=float,
    default=0.25,
    action='store',
    help='The lower bound fraction of the elastic threshold to designate read commute donors. Default=0.25.')

parser_screen.add_argument('--variant-caller', '--variant_caller',
    choices=['freebayes', 'bcftools'],
    default='freebayes',
    dest='variant_caller',
    action='store',
    help='Variant caller to use for SNP detection before lineage classification. Default=freebayes.')

parser_screen.add_argument('--lineage-profile-compound', '--compound-dataset', '--compound_dataset',
    default=False,
    dest='lineage_profile_compound',
    action='store_true',
    help='Force lineage profiling when multiple taxa are detected. Default=False because mixed datasets can confound lineage SNP profiles.')

parser_screen.add_argument('--snp-min-qual',
    type=float,
    default=30.0,
    dest='snp_min_qual',
    action='store',
    help='Minimum VCF QUAL for SNPs included in SNP JSON. Default=30.')

parser_screen.add_argument('--snp-min-dp',
    type=int,
    default=None,
    dest='snp_min_dp',
    action='store',
    help='Require INFO/DP >= this value for SNPs included in SNP JSON. Default=no INFO/DP filter.')

parser_screen.add_argument('--snp-min-missing-depth',
    type=int,
    default=10,
    dest='snp_min_missing_depth',
    action='store',
    help='Positions with samtools depth -aa coverage below this value are marked missing. Default=10.')

parser_screen.add_argument('--snp-allow-filtered',
    default=False,
    dest='snp_allow_filtered',
    action='store_true',
    help='Allow VCF records with FILTER values other than PASS/. in SNP JSON. Default=False.')

parser_screen.add_argument('--snp-accept-missing-qual',
    default=False,
    dest='snp_accept_missing_qual',
    action='store_true',
    help='Allow VCF records with missing QUAL in SNP JSON. Default=False.')

parser_screen.add_argument('--no-lineage-classify',
    default=False,
    dest='no_lineage_classify',
    action='store_true',
    help='Do not run lineage classification even when a database profile model is available. Default=False.')

parser_screen.add_argument('--profiles-dir',
    type=isDir,
    default=None,
    dest='profiles_dir',
    action='store',
    help='Directory containing profiles.json and profile references/models. Default=<database>/profiles, then sibling ../profiles if present.')

parser_screen.add_argument('--lineage-min-support',
    type=int,
    default=1,
    dest='lineage_min_support',
    action='store',
    help='Minimum matched loci required for hierarchical lineage descent. Default=1.')

parser_screen.add_argument('--lineage-min-support-fraction',
    type=float,
    default=0.75,
    dest='lineage_min_support_fraction',
    action='store',
    help='Minimum callable support fraction required for hierarchical lineage descent. Default=0.75.')

parser_screen.add_argument('--lineage-min-callable-fraction',
    type=float,
    default=0.5,
    dest='lineage_min_callable_fraction',
    action='store',
    help='Minimum callable marker fraction required for supported hierarchical lineage calls. Default=0.5.')

parser_screen.add_argument('--lineage-ambiguity-margin',
    type=float,
    default=0.1,
    dest='lineage_ambiguity_margin',
    action='store',
    help='Minimum decision-score margin used to resolve multiple compatible lineage children. Default=0.1.')

parser_screen.add_argument('--lineage-tie-delta',
    type=float,
    default=0.0,
    dest='lineage_tie_delta',
    action='store',
    help='Posterior delta used to report tied lineage calls. Default=0.0.')

parser_screen.add_argument('--lineage-disable-reference-marker-inference',
    default=False,
    dest='lineage_disable_reference_marker_inference',
    action='store_true',
    help='Disable reference-marker inference for canonical-as-empirical models. Default=False.')

parser_screen.add_argument('--lineage-disable-incomplete-descent',
    default=False,
    dest='lineage_disable_incomplete_descent',
    action='store_true',
    help='Disable incomplete hierarchical descent through partially callable lineage nodes. Default=False.')

parser_screen.add_argument('-o', '--output_prefix',
    type=str,
    default="Afanc_screen",
    action='store',
    help='Output prefix for this run. Default=results.')

parser_screen.add_argument('-k', '--run_key',
    default=False,
    action='store_true',
    help='Generate and use a random key for this runs run directory. Results will be stored in {output_prefix}.{random_key}')

# parser_screen.add_argument('-f', '--fetch_assemblies',
#     default=False,
#     action='store_true',
#     help='Fetch genome assemblies for species hits using from GenBank using the ENSEMBL software suite. Default=False; get assemblies from the autodb working directory where possible.')

parser_screen.add_argument('-t', '--threads',
    type=int,
    default=4,
    action='store',
    help='Number of threads to used for this run. Default=4.')

clean_group = parser_screen.add_mutually_exclusive_group()

clean_group.add_argument('-c', '--clean',
    default=False,
    action='store_true',
    help='Remove the mapping working directory from the output directory. Default=False.')

clean_group.add_argument('-s', '--superclean',
    default=False,
    action='store_true',
    help='Delete the entire output directory, leaving only log files and the results .json. Default=False.')

clean_group.add_argument('-a', '--no_map',
    default=False,
    action='store_true',
    help='Only perform metagenomic screening, do not proceed to mapping or variant calling. Default=False.')

parser_screen.set_defaults(func=runSubtool)

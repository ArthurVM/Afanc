from sys import stdout
from shutil import move, rmtree
from os import mkdir, chdir, path, listdir

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import vprint
from Afanc.utilities.makeWD import mkchdir

def download_genome(assembly, args):
    """ Download a genome from genbank using the entrez suite.
    """

    if not args.accessions:
        runline = f"esearch -db assembly -query \'{assembly}[organism] AND latest[filter]\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"

    else:
        runline = f"esearch -db assembly -query \'{assembly}\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"

    tmp_ftp_dir = command(runline, "DOWNLOAD_HITS").run_comm_quiet(1, args.stdout, args.stderr)

    ftp_dirs = [ f for f in tmp_ftp_dir[0].decode().split("\n") if f != '']

    if len(ftp_dirs) >= args.num_assemblies:
        ftp_dirs = ftp_dirs[:args.num_assemblies]

    # print(ftp_dirs)

    for ftp_dir in ftp_dirs:

        base = path.basename(ftp_dir.split("/")[-1])

        outfile = f"{base}_genomic.fna.gz"

        if path.exists(outfile):
            continue
        else:
            grepline = f"curl {ftp_dir}/{base}_genomic.fna.gz --output {outfile}"
            stdout, stderr = command(grepline, "DOWNLOAD_HITS").run_comm_quiet(1, args.stdout, args.stderr)


def parse_names_file(txt):
    """ Parses the line seperated text file containing species/accession IDs to download.
    """

    with open(txt, 'r') as fin:
        assemblyIDs = [ line.strip("\n") for line in fin.readlines() ]

    return assemblyIDs


def runGet_dataset(args):
    """ Get genome assemblies
    """

    subprocessID = "GET_ASSEMBLIES"
    vprint(
        subprocessID,
        f"Downloading assemblies...",
        "prYellow"
    )

    assemblyIDs = parse_names_file(args.ID_file)
    numIDs = len(assemblyIDs)

    print(f"Found {numIDs} IDs in {args.ID_file}\n")

    mkchdir(f"./{args.output_prefix}")

    for i, id in enumerate(assemblyIDs):

        stdout.write(f"\rDownloading assemblies for {id} ({i+1}/{numIDs})\r")
        stdout.flush()

        if not args.accessions:
            mkchdir(f"./{id.replace(' ', '_')}")

            download_genome(id, args)

            chdir("../")

        else:
            download_genome(id, 1, args.accessions)

    vprint(
        subprocessID,
        f"Done. Assemblies can be found in ./{args.output_prefix}.",
        "prGreen"
    )

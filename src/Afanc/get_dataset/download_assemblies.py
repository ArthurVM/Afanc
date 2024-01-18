import json

from sys import stdout, exit
from shutil import move, rmtree
from os import mkdir, chdir, path, listdir

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import vprint
from Afanc.utilities.makeWD import mkchdir

def download_genome(assembly, args):
    """ Download a genome from genbank using the entrez suite.
    """

    if not args.accessions:
        # runline = f"esearch -db assembly -query \'{assembly}[organism] AND latest[filter]\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
        runline = f"datasets summary genome taxon \"{assembly}\" --assembly-level complete \
                    --exclude-atypical \
                    --mag exclude \
                    --tax-exact-match \
                    --limit 100 \
                    --assembly-version latest \
                    --assembly-source GenBank"

    else:
        print(f"Accessions Download deprecated. Please use names.")
        exit(122)
        runline = f"esearch -db assembly -query \'{assembly}\' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"

    tmp_ftp_dir = command(runline, "DOWNLOAD_HITS").run_comm_quiet(1, args.stdout, args.stderr)

    accessions_list = []
    out = json.loads(tmp_ftp_dir[0].decode())

    ## check there are sequences on GenBank available for this species
    if out["total_count"] == 0:
        return
    
    ## check there are enough sequences to download
    ## if not, download them all
    if out["total_count"] < args.num_assemblies:
        num_assemblies = out["total_count"]
    else:
        num_assemblies = args.num_assemblies

    for k in out["reports"]:
        try:
            accessions_list.append([k["accession"], \
                                    k["organism"]["organism_name"], \
                                    k["average_nucleotide_identity"]["best_ani_match"]["ani"], \
                                    k["average_nucleotide_identity"]["best_ani_match"]["organism_name"], \
                                    k["assembly_info"]["assembly_level"], \
                                    k["assembly_stats"]["scaffold_n50"]])
        except:
            continue
        
    ## sort all accessions by N50 to select the highest quality
    sorted_accs_list = sorted(accessions_list, key=lambda x: x[5], reverse=True)
    
    accs_list = " ".join([acc[0] for acc in sorted_accs_list[:num_assemblies]])
    species_id = assembly.replace(' ', '_')

    dl_runline = f"datasets download genome accession {accs_list} --include genome --filename tmp_{species_id}.zip"
    stdout, stderr = command(dl_runline, "DOWNLOAD").run_comm_quiet(1, args.stdout, args.stderr)

    unzip_runline = f"unzip tmp_{species_id}.zip -d tmp_{species_id}; mv tmp_{species_id}/ncbi_dataset/data/*/*fna ./; rm -r ./tmp_{species_id}*"
    stdout, stderr = command(unzip_runline, "DOWNLOAD").run_comm_quiet(1, args.stdout, args.stderr)


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

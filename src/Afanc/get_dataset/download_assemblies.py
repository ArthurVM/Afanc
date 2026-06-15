import json

from pathlib import Path
from shutil import move, rmtree
from os import chdir, getcwd

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import vprint
from Afanc.utilities.makeWD import mkchdir

def download_genome(assembly, args):
    """ Download a genome from genbank using the entrez suite.
    """

    if args.accessions:
        raise ValueError("Accession downloads are deprecated. Please provide species names instead.")

    runline = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        assembly,
        "--assembly-level",
        "complete",
        "--exclude-atypical",
        "--mag",
        "exclude",
        "--tax-exact-match",
        "--limit",
        "100",
        "--assembly-version",
        "latest",
        "--assembly-source",
        "GenBank",
    ]

    tmp_ftp_dir = command(runline, "DOWNLOAD_HITS").run_comm_quiet(1, args.stdout, args.stderr)

    accessions_list = []
    out = json.loads(tmp_ftp_dir[0].decode())

    ## check there are sequences on GenBank available for this species
    if out["total_count"] == 0:
        return []
    
    mkchdir(f"./{assembly.replace(' ', '_')}")
    
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
        except (KeyError, TypeError):
            continue

    if not accessions_list:
        return []
        
    ## sort all accessions by N50 to select the highest quality
    sorted_accs_list = sorted(accessions_list, key=lambda x: x[5], reverse=True)
    
    selected_accessions = [acc[0] for acc in sorted_accs_list[:num_assemblies]]
    species_id = assembly.replace(' ', '_')

    zip_path = Path(f"tmp_{species_id}.zip")
    extract_dir = Path(f"tmp_{species_id}")
    dl_runline = [
        "datasets",
        "download",
        "genome",
        "accession",
        *selected_accessions,
        "--include",
        "genome",
        "--filename",
        str(zip_path),
    ]
    command(dl_runline, "DOWNLOAD").run_comm_quiet(1, args.stdout, args.stderr)

    command(["unzip", str(zip_path), "-d", str(extract_dir)], "DOWNLOAD").run_comm_quiet(1, args.stdout, args.stderr)

    for fasta_path in (extract_dir / "ncbi_dataset" / "data").glob("*/*.fna"):
        move(str(fasta_path), ".")

    if extract_dir.exists():
        rmtree(extract_dir)
    if zip_path.exists():
        zip_path.unlink()

    return selected_accessions


def parse_names_file(txt):
    """ Parses the line seperated text file containing species/accession IDs to download.
    """

    with open(txt, 'r') as fin:
        assemblyIDs = [line.strip() for line in fin if line.strip()]

    return assemblyIDs


def runGet_dataset(args):
    """ Get genome assemblies
    """

    subprocessID = "GET_ASSEMBLIES"
    vprint(
        subprocessID,
        f"Downloading assemblies...",
        "prYellow",
        args.stdout,
    )

    if args.accessions:
        raise ValueError("Accession downloads are deprecated. Please provide species names instead.")

    cwd = getcwd()
    assemblyIDs = parse_names_file(args.ID_file)
    numIDs = len(assemblyIDs)

    print(f"Found {numIDs} IDs in {args.ID_file}\n", file=args.stdout, flush=True)

    mkchdir(f"{cwd}/{args.output_prefix}")

    for i, assembly_id in enumerate(assemblyIDs):

        if "/" in assembly_id or "\\" in assembly_id:
            print(f"Invalid characters found in {assembly_id}. Skipping...", file=args.stdout, flush=True)
            continue

        print(f"Downloading assemblies for {assembly_id} ({i+1}/{numIDs})", file=args.stdout, flush=True)

        try:
            download_genome(assembly_id, args)
        except Exception as exc:
            print(f"Something went wrong attempting to download {assembly_id}: {exc}", file=args.stderr, flush=True)

        chdir(f"{cwd}/{args.output_prefix}")

    vprint(
        subprocessID,
        f"Done. Assemblies can be found in ./{args.output_prefix}.",
        "prGreen",
        args.stdout,
    )

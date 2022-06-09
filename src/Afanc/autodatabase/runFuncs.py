from shutil import move, rmtree
from os import mkdir, chdir, path, listdir
from collections import defaultdict

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import vprint
from Afanc.utilities.makeWD import mkchdir


def runAutoDB(args):
    """ Run autodatabase pipeline
    """
    from Afanc.utilities.makeWD import initAutoDBDirStructure

    subprocessID = "MAIN"
    vprint(
        subprocessID,
        f"Starting autodatabase pipeline using FASTA files found in {args.fastaDir}",
        "prYellow"
    )

    ## initialise autoDB directory structure
    initAutoDBDirStructure(args)

    ## construct root directory to deposit results
    mkchdir(args.autoDB_WDir)

    fasta_db_path = args.fastaDir

    ## download ncbi taxonomy and preprocess fastas
    fasta_dict = preprocessing(args, fasta_db_path)

    ## quality control assemblies and select highest quality ones from each taxa
    assemblyQC(args, fasta_dict)

    ## make the kraken2 database from the quality controlled assemblies
    makeK2db(args)

    ## make a Krona chart for pleasing visualisation
    makeKronaChart(args)

    vprint(
        "FINISHED",
        "Pipeline terminated successfully.",
        "prGreen"
    )


def preprocessing(args, fasta_db_path):
    """ Download NCBI taxonomy and add taxa to nodes and names files
    """

    from Afanc.autodatabase.prepareNewFasta import getTaxonomy
    from Afanc.autodatabase.taxadd import taxadd_Main

    subprocessID = "PREPROCESSING"
    vprint(
        subprocessID,
        f"Downloading NCBI-taxonomy files...",
        "prYellow"
    )

    mkchdir("ncbi_taxonomy")

    ## get ncbi taxonomy for given date
    names, nodes = getTaxonomy(args, "2020-05-01")

    chdir(args.autoDB_WDir)

    mkchdir(args.fasta_WDir)

    fasta_dict = taxadd_Main(args, fasta_db_path, names)

    chdir(args.autoDB_WDir)

    return fasta_dict


def assemblyQC(args, fasta_dict):
    """ Perform QC on assemblies using MASH then move high quality assemblies to the cleanFasta directory
    """
    from Afanc.autodatabase.assemblyQC import mash, buildMatrix, fastaMove
    from Afanc.autodatabase.makeFastaDirJSON import make_fasta_dir_JSON

    subprocessID = "QC"
    vprint(
        subprocessID,
        f"Performing quality control on FASTA files using MASH...",
        "prYellow"
    )

    ## construct a taxa dict where keys are the taxonID, and values a list containing paths to all fastas in that taxon
    taxa_dict = defaultdict(list)
    input_fasta_list = listdir(args.fasta_WDir)

    print(f"{len(input_fasta_list)} assemblies found in {args.fastaDir}")

    ## make the species name fasta file json
    make_fasta_dir_JSON(args.fastaDir)

    for fasta in input_fasta_list:
        taxon_id = fasta.split("_")[0]
        fasta_path = path.join(args.fasta_WDir, fasta)
        taxa_dict[taxon_id].append(fasta_path)

    ## generate directory structure
    mkchdir(args.qc_WDir, 0)
    mkchdir(args.mash_WDir, 0)
    mkchdir(args.cleanFasta_WDir, 0)

    for taxon_id, fastas in taxa_dict.items():
        ## calculate minhash distances for each taxon using Mash
        chdir(args.mash_WDir)
        mashdist_out = mash(args, taxon_id, fastas)

        chdir(args.qc_WDir)

        ## range to take around the mode of the average mash distance, default is 10%
        modeRange = args.mode_range

        calcArray, tax, modeVal = buildMatrix(args, mashdist_out)

        ## if the number of samples < 3 then continue
        if calcArray is None:
            continue

        ## select the best assembly using minhash distance
        fastaMove(args, calcArray, tax, modeVal, modeRange)

    chdir(args.autoDB_WDir)


def makeK2db(args):
    from Afanc.utilities.runCommands import command

    subprocessID = "MAKE_K2_DB"
    vprint(
        subprocessID,
        f"Making Kraken2 database using quality controlled assemblies...",
        "prYellow"
    )

    mkchdir(args.kraken2_WDir)
    mkchdir("taxonomy", 0)
    move(f"{args.autoDB_WDir}/ncbi_taxonomy/names.dmp", "./taxonomy")
    move(f"{args.autoDB_WDir}/ncbi_taxonomy/nodes.dmp", "./taxonomy")

    for fasta in listdir(args.cleanFasta_WDir):
        fasta_path = path.join(args.cleanFasta_WDir, fasta)
        k2build_line = f"kraken2-build --add-to-library {fasta_path} --db ."
        stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm_quiet(1, args.stdout, args.stderr)

    k2build_line = "kraken2-build --db . --download-taxonomy"
    stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm(1, args.stdout, args.stderr)

    k2build_line = f"kraken2-build --build --threads 10 --db ."
    stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm(1, args.stdout, args.stderr)

    k2build_line = f"kraken2-inspect --db . > database.txt"
    stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm(1, args.stdout, args.stderr)

    ## remove unnecessary directories
    rmtree("./library")
    rmtree("./taxonomy")

    chdir(args.autoDB_WDir)


def makeKronaChart(args):
    from Afanc.utilities.runCommands import command

    subprocessID = "KRONA"
    vprint(
        subprocessID,
        f"Making Krona chart...",
        "prYellow"
    )

    mkchdir(args.krona_WDir)

    kronachart_line = f"ktImportTaxonomy -m 3 -t 5 {args.kraken2_WDir}/database.txt -o database.html"
    stdout, stderr = command(kronachart_line, "KRONA").run_comm(1, args.stdout, args.stderr)

    chdir(args.autoDB_WDir)

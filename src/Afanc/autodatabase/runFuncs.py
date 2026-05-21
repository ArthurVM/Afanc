import json
from shutil import move, rmtree
from os import mkdir, chdir, path, listdir, rename, getcwd, remove
from shutil import copyfile
from collections import defaultdict

from ..utilities.runCommands import command
from ..utilities.generalUtils import vprint
from ..utilities.makeWD import mkchdir


def runAutoDB(args):
    """ Run autodatabase pipeline
    """
    from ..utilities.makeWD import initAutoDBDirStructure
    from ..utilities.getVersions import getVersionsAutodatabase

    ## check the working directory doesn't already exist
    if path.exists(args.output_prefix):
        vprint("ERROR", f"Output directory {args.output_prefix} already exists! Aborting to prevent accidental data loss...", "prRed")
        exit(1)

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

    fasta_db_path = str(args.fastaDir)

    ## capture python package and software versions for the autodatabase module in a JSON file
    getVersionsAutodatabase(args)

    ## download ncbi taxonomy and preprocess fastas
    fasta_dict, mapping_dict = preprocessing(args, fasta_db_path)

    ## quality control assemblies and select highest quality ones from each taxon
    assemblyQC(args, fasta_dict, mapping_dict)

    ## make the kraken2 database from the quality controlled assemblies
    makeK2db(args)

    ## make the variant index from quality controlled assemblies
    makeVariantIndex(args)

    ## make a Krona chart for pleasing visualisation
    makeKronaChart(args)

    ## clean the output directory
    cleanOutdir(args)

    vprint(
        "FINISHED",
        "Afanc-autodatabase finished successfully.",
        "prGreen"
    )


def preprocessing(args, fasta_db_path):
    """ Download NCBI taxonomy and add taxa to nodes and names files
    """

    from .prepareNewFasta import getTaxonomy
    from .taxadd import taxadd_Main
    
    subprocessID = "PREPROCESSING"
    mkchdir("ncbi_taxonomy")

    if args.ncbi_tax_db == False:
        vprint(
            subprocessID,
            f"Downloading NCBI-taxonomy files...",
            "prYellow"
        )

        ## get ncbi taxonomy for given date
        names, nodes = getTaxonomy(args, args.ncbi_date)
    
    else:
        vprint(
            subprocessID,
            f"Getting local NCBI-taxonomy files...",
            "prYellow"
        )

        ## copy required files
        names = f"{args.autoDB_WDir}/ncbi_taxonomy/names.dmp"
        nodes = f"{args.autoDB_WDir}/ncbi_taxonomy/nodes.dmp"
        merged = f"{args.autoDB_WDir}/ncbi_taxonomy/merged.dmp"
        copyfile(f"{args.ncbi_tax_db}/names.dmp", names)
        copyfile(f"{args.ncbi_tax_db}/nodes.dmp", nodes)
        copyfile(f"{args.ncbi_tax_db}/merged.dmp", merged)

    chdir(args.autoDB_WDir)

    mkchdir(args.fasta_WDir)

    fasta_dict, mapping_dict = taxadd_Main(args.fasta_WDir, fasta_db_path, names, nodes)

    chdir(args.autoDB_WDir)

    return fasta_dict, mapping_dict


def assemblyQC(args, fasta_dict, mapping_dict):
    """ Perform QC on assemblies using MASH then move high quality assemblies to the cleanFasta directory
    """
    from .assemblyQC import mash, buildMatrix, fastaMove
    from .makeFastaDirJSON import makeFastaDirJSON

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
    makeFastaDirJSON(str(args.fastaDir))

    for fasta in input_fasta_list:
        taxon_id = fasta.split("_")[0]
        fasta_path = path.join(args.fasta_WDir, fasta)
        taxa_dict[taxon_id].append(fasta_path)

    ## generate directory structure
    mkchdir(args.qc_WDir, 0)
    mkchdir(args.mash_WDir, 0)
    mkchdir(args.cleanFasta_WDir, 0)

    db_assemblies = defaultdict(list)

    for taxon_id, fastas in taxa_dict.items():
        ## calculate minhash distances for each taxon using Mash
        chdir(args.mash_WDir)
        mashdist_out = mash(args, taxon_id, fastas)

        chdir(args.qc_WDir)

        ## range to take around the mode of the average mash distance, default is 10%
        capture_range = args.mode_range

        calcArray, tax, centroid = buildMatrix(args, mashdist_out)

        ## if the number of samples < 3 then continue
        if calcArray is None:
            db_assemblies[taxon_id] = [path.basename(f) for f in fastas]
            continue

        ## select the best assembly using minhash distance
        clean_fastas = fastaMove(args, calcArray, tax, centroid, capture_range)

        ## remove base paths
        filenames = [path.basename(f) for f in clean_fastas]
        db_assemblies[taxon_id] = filenames

    chdir(args.autoDB_WDir)

    with open(f"taxID_mappings.json", "w") as fout:
        json.dump(mapping_dict, fout, indent=4)
    with open(f"./fastas_in_DB.json", "w") as fout:
        json.dump(db_assemblies, fout, indent=4)


def makeK2db(args):

    from ..utilities.runCommands import command

    subprocessID = "MAKE_K2_DB"
    vprint(
        subprocessID,
        f"Making Kraken2 database using quality controlled assemblies...",
        "prYellow"
    )

    mkchdir(args.kraken2_WDir)

    mkchdir("taxonomy", 0)

    ## Kraken2 only needs the prepared local taxonomy here because FASTA
    ## headers already contain explicit kraken:taxid fields. Calling
    ## --download-taxonomy can fail on restricted rsync/ftp configurations and
    ## would overwrite the taxonomy edited by taxadd.
    for taxdump_file in ["names.dmp", "nodes.dmp", "merged.dmp"]:
        source_taxdump = path.join(args.autoDB_WDir, "ncbi_taxonomy", taxdump_file)
        if path.exists(source_taxdump):
            copyfile(source_taxdump, path.join(args.kraken2_WDir, "taxonomy", taxdump_file))

    for fasta in listdir(args.cleanFasta_WDir):
        fasta_path = path.join(args.cleanFasta_WDir, fasta)
        k2build_line = f"kraken2-build --add-to-library {fasta_path} --db ."
        stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm_quiet(1, args.stdout, args.stderr)

    ## build and inspect the kraken2 database
    k2build_line = f"kraken2-build --build --threads {args.threads} --db ."
    stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm(1, args.stdout, args.stderr)

    k2build_line = f"kraken2-inspect --db . > database.txt"
    stdout, stderr = command(k2build_line, "KRAKEN2-BUILD").run_comm(1, args.stdout, args.stderr)

    ## remove large unnecessary directories
    rmtree(f"{args.kraken2_WDir}/library")
    rmtree(f"{args.kraken2_WDir}/taxonomy")

    chdir(args.autoDB_WDir)


def makeVariantIndex(args):
    """ Generate the Mash Variant Index used by screen read deconvolution.
    """
    from .makeVariantIndex import makeVariantIndex

    subprocessID = "GEN-VARIANT-INDEX"
    vprint(
        subprocessID,
        f"Calculating all-vs-all Mash distances and generating the Mash Variant Index...",
        "prYellow"
    )

    mkchdir(args.variant_index_WDir)

    makeVariantIndex(args)

    chdir(args.autoDB_WDir)


def makeKronaChart(args):

    from ..utilities.krona import run_krona_from_kraken_report

    subprocessID = "KRONA"
    vprint(
        subprocessID,
        f"Making Krona chart...",
        "prYellow"
    )

    mkchdir(args.krona_WDir)

    run_krona_from_kraken_report(
        path.join(args.kraken2_WDir, "database.txt"),
        path.join(args.krona_WDir, "database.html"),
        stdout=args.stdout,
        stderr=args.stderr,
        subprocess_id=subprocessID,
    )

    chdir(args.autoDB_WDir)


def cleanOutdir(args):
    """ Cleans the output directory according to provided arguments.

    clean : remove the mash working directory.
    superclean : remove everything not necessary for running Afanc screen
    """

    def rm_dir(dir):
        """ Remove a directory
        """
        ## check directory exists
        if path.isdir(dir):
            rmtree(dir, ignore_errors=True)
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove {dir}, directory not found.")

    def clean(args):
        rm_dir(args.qc_WDir)
        rm_dir(args.mash_WDir)

    def superclean(args):
        ## run the standard clean
        clean(args)
        rm_dir(args.fasta_WDir)
        rm_dir(f"{args.autoDB_WDir}/ncbi_taxonomy/")

    if args.clean:
        clean(args)

    elif args.superclean:
        superclean(args)

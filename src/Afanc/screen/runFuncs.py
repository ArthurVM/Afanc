import pysam
import sys
from importlib.metadata import version
from shutil import move, rmtree
from os import chdir, path, listdir, mkdir, remove
from collections import defaultdict

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import vprint


def runScreen(args):
    """ Run the Afanc-screen pipeline
    """
    from Afanc.utilities.makeWD import genScreenDirStructure, checkautodbWD

    ## generate the working directory structure for this run
    genScreenDirStructure(args)

    ## check if database is malformed
    checkautodbWD(args)

    ## screen with kraken2
    runKraken2(args)

    ## parse kraken2 report to a json
    getHits(args)

    ## map to my genomes
    map2Hits(args)

    ## generate the final report
    final_report = makeFinalReport(args)

    if args.clean or args.superclean:
        final_report = clean_outdir(args, final_report)

    vprint("FINISHED", f"Final report can be found at {final_report}\n", "prGreen")


def runKraken2(args):
    """ run kraken2 on query fqs with using the database
    """
    subprocessID="KRAKEN"

    chdir(args.k2WDir)

    vprint(
        subprocessID,
        f"Running Kraken2...",
        "prYellow",
    )

    runline = f"kraken2 \
     --db {args.k2_database} \
     --minimum-base-quality 15 \
     --confidence 0.05 \
     --threads {args.threads} \
     --output {args.output_prefix}.k2.txt \
     --report {args.output_prefix}.k2.report.txt \
     --paired {' '.join(args.fastq)}"

    stdout, stderr = command(runline, "KRAKEN").run_comm(1, args.stdout, args.stderr)

    chdir(args.runWDir)


def getHits(args):
    """ parses the kraken2 report and then gets assemblies pertaining to identified species
    """
    from Afanc.screen.report.parseK2report import parseK2reportMain
    from Afanc.screen.getGenomes import getAccessions, getGenomesbyAcc, getGenomesbyName, getLocalGenomes, get_hitIDs
    from Afanc.utilities.generalUtils import gendbdict

    ## generate the taxID database dictionary from the autodatabase cleanFasta directory
    dbdict = gendbdict(args.cleanFasta)

    ## parse kraken2 report
    chdir(args.reportsDir)
    out_json = parseK2reportMain(args, dbdict)

    ## check if there were any hits
    if out_json is None:
        vprint(
            "MAIN",
            "No hits detected exceeding the given thresholds. Exiting...",
            "prRed",
        )
        sys.exit(1)

    chdir(args.bt2WDir)

    ## GET GENOMES BRANCH
    ## branch for dealing with how hit genomes are fetched: local using the autodb results dir,
    ## or from genbank using the ensembl suite

    if args.fetch_assemblies:
        ## download hits from Genbank
        ## TODO : this branch is broken and results in a malformed results JSON
        subprocessID="DOWNLOADING_HITS"

        vprint(
            subprocessID,
            "Downloading assemblies from Genbank...",
            "prYellow",
        )
        assembly_ids = get_hitIDs(out_json)
        assembly_dict = getGenomesbyName(assembly_ids, args)

    else:
        ## get genomes from autodatabase results directory: /selectFasta_autoDatabase_cleanFasta/
        ## TODO : this branch is bugged due to inconsistent ncbi_taxID usage
        ## TEMP FIX : download genomes which cannot be found by ncbi_taxID
        subprocessID="GET_HITS"

        vprint(
            subprocessID,
            "Getting assemblies from local autodatabase directory...",
            "prYellow",
        )
        getLocalGenomes(out_json, args)

    chdir(args.runWDir)


def map2Hits(args):
    """ takes the hits from the kraken2 module
    """
    import gzip
    import json
    from numpy import sum

    from Afanc.utilities.generalUtils import parseBT2out
    from Afanc.utilities.modifyFasta import modifyFasta
    from Afanc.screen.maths.mappingMetrics import gini, genomeSize, breadthofCoverage, meanDOC, medianDOC

    subprocessID="MAP2HITS"

    vprint(
        subprocessID,
        "Mapping reads to hits...",
        "prYellow",
    )

    ## index genomes
    chdir(args.bt2WDir)
    mkdir("Hits")

    accessions = {}

    for assembly in listdir(args.bt2WDir):
        if assembly.endswith(".fna.gz") or assembly.endswith(".fna") or assembly.endswith(".fa.gz") or assembly.endswith(".fa"):
            ## store accession assembly pairs
            tmp_acc = path.splitext(assembly)[0]

            ## strip out genomic tag from the accession
            if tmp_acc.endswith("_genomic"):
                 accession = tmp_acc.split("_genomic")[0]
            else:
                accession = tmp_acc

            accessions[accession] = assembly
            modifyFasta(assembly, accession)

    ## RUN BOWTIE2

    ## generate the BT2 index from all assembly hits
    assemblies_line = ",".join([f"{acc}.tmp.fa" for acc, assembly in accessions.items()])
    bt2_index = "Hits/Hits"

    buildline = f"bowtie2-build {assemblies_line} {bt2_index}"
    stdout, stderr = command(buildline, "BUILD").run_comm(1, args.stdout, args.stderr)

    chdir(args.bt2WDir)

    ## map reads to the aggregated index
    mapline = f"bowtie2 -p {args.threads} -x {bt2_index} --{args.mapping_sensitivity} -1 {args.fastq[0]} -2 {args.fastq[1]} > Hits.sam"
    stdout, stderr = command(mapline, "MAP").run_comm(1, args.stdout, args.stderr)

    ## capture general mapping stats
    with open(f"./Hits.mapstats.txt", 'w') as fout:
        fout.write(f"{stdout.decode()}\n\n{stderr.decode()}")

    for accession, assembly in accessions.items():

        ## store data in a dict for output as a json
        datadict = defaultdict(dict)

        datadict["warnings"] = {}
        datadict["input_data"]["reference"] = f"{args.bt2WDir}{bt2_index}"
        datadict["input_data"]["fastq_1"] = f"{args.fastq[0]}"
        datadict["input_data"]["fastq_2"] = f"{args.fastq[1]}"

        ## strip out genomic tag from the accession
        # if assembly.endswith("_genomic"):
        #      accession = assembly.split("_genomic")[0]

        grepline = f"grep '{accession}' Hits.sam | samtools view -bh - | samtools sort - > {accession}.sorted.bam"
        stdout, stderr = command(grepline, "MAP").run_comm(1, args.stdout, args.stderr)

        ## calculate the depth at each covered position
        depthline = f"samtools depth {accession}.sorted.bam"

        depthstdout, depthstderr = command(depthline, "DEPTH").run_comm(1, None, args.stderr)
        covarray = [i.split('\t') for i in depthstdout.decode().split('\n')]

        ## check covarray is populated
        if len(covarray) <= 1:
            ## if not, then no reads map uniquely to this assembly
            datadict["warnings"]["no_unique_map"] = f"No reads map uniquely to assembly {assembly}."

        else:
            ## if yes then continue with read distribution calculations
            ## calculate genome size
            genomesize = genomeSize(assembly)

            ## calculate the mean and median depth of coverage, NB this is only the DOC at covered positions, and does not include non-covered positions
            meandoc = meanDOC(covarray)
            mediandoc = medianDOC(covarray)
            boc = breadthofCoverage(covarray, genomesize)

            ## calculate the Gini coefficient of distribution
            gini_co = gini(covarray)

            datadict["map_data"]["mean_DOC"] = meandoc
            datadict["map_data"]["median_DOC"] = mediandoc
            datadict["map_data"]["proportion_cov"] = boc
            datadict["map_data"]["gini"] = gini_co

        with open(f"{args.reportsDir}/{accession}.mapstats.json", 'w') as fout:
            json.dump(datadict, fout, indent = 4)

    chdir(args.runWDir)

def makeFinalReport(args):
    """ make final report from the kraken2 and bowtie2 reports
    """

    import json
    from os import listdir

    from Afanc.utilities.get_versions import get_versions_screen

    subprocessID = "REPORT"
    vprint(
        subprocessID,
        f"Generating final report...",
        "prYellow"
    )

    jsondict = defaultdict(dict)

    datadict = jsondict["results"] = {}

    ## initialise dictionary key
    if "Detection_events" not in datadict:
        datadict["Detection_events"] = []

    ## collect arguments
    jsondict["arguments"] = {
        "database" : args.k2_database,
        "fastqs" : args.fastq,
        "num_threshold" : args.num_threshold,
        "pct_threshold" : args.pct_threshold,
        "output_prefix" : args.output_prefix,
    }

    ## collect versions
    jsondict["versions"] = {
        "screen" : get_versions_screen()
    }

    ## collect k2 json reports
    bt2_reportsbox = [g for g in listdir(args.reportsDir) if g.endswith("mapstats.json")]

    ## read general report
    with open(f"{args.reportsDir}/{args.output_prefix}.k2.json", 'r') as k2fin:
        k2data = json.load(k2fin)

        for event in k2data["Detection_events"]:

            ## initialise warnings
            event["warnings"] = []

            ## block to deal with most likely variants
            if "closest_variant" in event:
                variant_flag = True
                assembly = event["closest_variant"]["assembly"]
                taxon_id = str(event["closest_variant"]["taxon_id"])
            else:
                variant_flag = False
                assembly = event["assembly"]
                taxon_id = str(event["taxon_id"])

            ## block dealing with hits which have an accompanying assembly which reads were mapped to
            if not assembly == None:

                ## find Bowtie2 report
                ## NOTE: there should only ever be 1 element in this list, since taxon_id should be unique to each taxon
                report_json = [report for report in bt2_reportsbox if report.startswith(taxon_id)][0]

                ## read Bowtie2 report
                with open(f"{args.reportsDir}/{report_json}", 'r') as bt2fin:
                    bt2data = json.load(bt2fin)

                    ## check the that the no_unique_map warning doesnt exist
                    if "no_unique_map" not in bt2data["warnings"]:

                        ## add fields to json dict
                        ## if there is a likely variant, add to that subdict
                        if variant_flag:
                            event["closest_variant"]["mean_DOC"] = bt2data["map_data"]["mean_DOC"]
                            event["closest_variant"]["median_DOC"] = bt2data["map_data"]["median_DOC"]
                            event["closest_variant"]["reference_cov"] = bt2data["map_data"]["proportion_cov"]
                            event["closest_variant"]["gini"] = bt2data["map_data"]["gini"]

                        ## otherwise add to the base dict
                        else:
                            event["mean_DOC"] = bt2data["map_data"]["mean_DOC"]
                            event["median_DOC"] = bt2data["map_data"]["median_DOC"]
                            event["reference_cov"] = bt2data["map_data"]["proportion_cov"]
                            event["gini"] = bt2data["map_data"]["gini"]

                    ## if it does, append it to the final report warnings
                    else:
                        event["warnings"].append(bt2data["warnings"])

            ## block dealing with hits which do not have an accompanying assembly
            else:
                event['warnings'].append(f"No assembly could be found for hit on {event['name']}.")

            datadict["Detection_events"].append(event)

            final_report = path.abspath(f"./{args.output_prefix}.json")

    with open(final_report, "w") as fout:
        json.dump(jsondict, fout, indent = 4)

    return final_report


def clean_outdir(args, final_report):
    """ Cleans the output directory according to provided arguments.

    clean : remove the bt2 working directory.
    superclean : move the results JSON to ./ and remove the entire output directory
    """

    if args.clean:
        ## check the bt2WDir exists
        if path.isdir(args.bt2WDir):
            rmtree(args.bt2WDir, ignore_errors=True)
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove{args.bt2WDir}, directory not found.")

        return final_report

    elif args.superclean:
        ## check the root run directory exists
        if path.isdir(args.runWDir):
            new_final_report_path = path.join(args.baseRunDir, final_report.split("/")[-1])
            move(final_report, new_final_report_path)
            rmtree(args.runWDir, ignore_errors=True)
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove {args.runWDir}, directory not found.")

        return new_final_report_path

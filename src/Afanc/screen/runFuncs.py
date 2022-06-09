import pysam
from os import chdir, path, listdir, mkdir
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
    makeFinalReport(args)

    vprint("FINISHED", "Pipeline terminated successfully.", "prGreen")


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
        ## TODO : this branch is broken due to inconsistent ncbi_taxID usage
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
    from Afanc.screen.maths.mappingMetrics import gini, genomeSize, breadthofCoverage, meanDOC

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
        if assembly.endswith(".fna.gz") or assembly.endswith(".fna"):
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
    mapline = f"bowtie2 -p {args.threads} -x {bt2_index} --very-sensitive -1 {args.fastq[0]} -2 {args.fastq[1]} > Hits.sam"
    stdout, stderr = command(mapline, "MAP").run_comm(1, args.stdout, args.stderr)

    ## capture general mapping stats
    with open(f"./Hits.mapstats.txt", 'w') as fout:
        fout.write(f"{stdout.decode()}\n\n{stderr.decode()}")

    for accession, assembly in accessions.items():

        ## store data in a dict for output as a json
        datadict = defaultdict(dict)

        ## strip out genomic tag from the accession
        # if assembly.endswith("_genomic"):
        #      accession = assembly.split("_genomic")[0]

        grepline = f"grep '{accession}' Hits.sam | samtools view -bh - | samtools sort - > {accession}.sorted.bam"
        stdout, stderr = command(grepline, "MAP").run_comm(1, args.stdout, args.stderr)

        ## calculate the depth at each covered position
        depthline = f"samtools depth {accession}.sorted.bam"
        depthstdout, depthstderr = command(depthline, "DEPTH").run_comm(1, None, args.stderr)
        covarray = [i.split('\t') for i in depthstdout.decode().split('\n')]

        ## calculate genome size
        genomesize = genomeSize(assembly)

        ## calculate the mean depth of coverage, NB this is only the DOC at covered positions, and does not include non-covered positions
        meandoc = meanDOC(covarray)
        boc = breadthofCoverage(covarray, genomesize)

        ## calculate the Gini coefficient of distribution
        gini_co = gini(covarray)

        datadict["input_data"]["reference"] = f"{args.bt2WDir}{bt2_index}"
        datadict["input_data"]["fastq_1"] = f"{args.fastq[0]}"
        datadict["input_data"]["fastq_2"] = f"{args.fastq[1]}"

        datadict["map_data"]["mean_DOC"] = meandoc
        datadict["map_data"]["proportion_cov"] = boc
        datadict["map_data"]["gini"] = gini_co

        with open(f"{args.reportsDir}/{accession}.mapstats.json", 'w') as fout:
            json.dump(datadict, fout, indent = 4)

    chdir(args.runWDir)


import pysam
from os import chdir, path, listdir, mkdir
from collections import defaultdict

from Afanc.utilities.generalUtils import vprint

def makeFinalReport(args):
    """ make final report from the kraken2 and bowtie2 reports
    """

    import json
    from os import listdir

    subprocessID = "REPORT"
    vprint(
        subprocessID,
        f"Generating final report...",
        "prYellow"
    )

    jsondict = defaultdict(dict)

    datadict = jsondict["results"] = {}
    jsondict["arguments"] = { "database" : args.k2_database,
        "fastqs" : args.fastq,
        "pct_threshold" : args.pct_threshold,
        "num_threshold" : args.num_threshold,
        "output_prefix" : args.output_prefix,
    }

    ## collect k2 json reports
    reportsbox = [g for g in listdir(args.reportsDir) if g.endswith("mapstats.json")]

    with open(f"{args.reportsDir}/{args.output_prefix}.k2.json", 'r') as k2fin:
        k2data = json.load(k2fin)

        for event in k2data["Detection_events"]:
            for report in reportsbox:
                if str(event['taxon_id']) in report:
                    report_json = report

            with open(f"{args.reportsDir}/{report_json}", 'r') as bt2fin:      ## TODO: sort out this line to make it more robust, i.e. remove _genomic
                bt2data = json.load(bt2fin)

                ## add fields to json dict
                ## if there is a likely variant, add to that subdict
                if "most_likely_variant" in event:
                    event["most_likely_variant"]["mean_DOC"] = bt2data["map_data"]["mean_DOC"]
                    event["most_likely_variant"]["reference_cov"] = bt2data["map_data"]["proportion_cov"]
                    event["most_likely_variant"]["gini"] = bt2data["map_data"]["gini"]
                ## otherwise add to the base dict
                else:
                    event["mean_DOC"] = bt2data["map_data"]["mean_DOC"]
                    event["reference_cov"] = bt2data["map_data"]["proportion_cov"]
                    event["gini"] = bt2data["map_data"]["gini"]

                ## initialise dictionary key
                if "Detection_events" not in datadict:
                    datadict["Detection_events"] = []

                datadict["Detection_events"].append(event)

    with open(f"./{args.output_prefix}.json", "w") as fout:
        json.dump(jsondict, fout, indent = 4)

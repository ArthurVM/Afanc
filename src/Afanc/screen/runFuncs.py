from sys import exit
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
    out_json = runFPScreen(args)

    ## if no_map is True then exit Afanc screen
    if args.no_map:
        final_report = makeFinalReport(args, None, None)
        vprint("FINISHED", f"no_map mode finished. Metagenomic report can be found in {final_report}\n", "prGreen")
        return 0

    ## parse kraken2 report to a json
    hits_json, variant_species = getHits(args, out_json)

    ## map to my genomes
    variant_bams, reports = map2Hits(args, variant_species)

    ## check to see if variant files were provided
    if args.variant_profile != False:
        ## perform variant profiling
        variant_profile = profileVariants(args, variant_bams)
    else:
        variant_profile = None

    ## generate the final report
    final_report = makeFinalReport(args, variant_profile, reports)

    if args.clean or args.superclean:
        final_report = cleanOutdir(args, final_report)

    vprint("FINISHED", f"Final report can be found at {final_report}\n", "prGreen")

    return 0

def runFPScreen(args):
    """ run kraken2 on query fqs with using the database
    """
    from .report.parseK2report import parseK2reportMain
    from Afanc.utilities.generalUtils import gendbdict

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

    command(runline, "KRAKEN").run_comm(0, args.stdout, args.stderr)

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
        exit(1)

    chdir(args.runWDir)

    return out_json


def getHits(args, out_json):
    """ parses the kraken2 report and then gets assemblies pertaining to identified species
    """
    from .getGenomes import getAccessions, getGenomesbyAcc, getGenomesbyName, getLocalGenomes, get_hitIDs

    subprocessID="GET_HITS"

    chdir(args.bt2WDir)

    vprint(
        subprocessID,
        "Getting assemblies from local autodatabase directory...",
        "prYellow",
    )
    variant_species = getLocalGenomes(out_json, args)

    chdir(args.runWDir)

    return out_json, variant_species


def map2Hits(args, variant_species):
    """ takes the hits from the kraken2 module
    """
    import gzip
    import json
    from numpy import sum

    from .mapHits import make_accessions_dict, modify_fastas, gen_index, run_map

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
    variant_sams = {}

    assemblies_for_mapping = make_accessions_dict(args, variant_species)

    modify_fastas(assemblies_for_mapping)

    ## RUN BOWTIE2

    bt2_index = gen_index(args, assemblies_for_mapping)

    variant_bams, reports = run_map(args, bt2_index, assemblies_for_mapping)

    chdir(args.runWDir)

    return variant_bams, reports


def profileVariants(args, variant_bams):
    """ Take a set of variant profiles stored in BED format, and screen the fastqs for polymorphisms within the BED file.

    The variants BED file must be structured as a tab seperated file with fields:

            {chr}   {start}    {end}    {variantID}    {variant_sequence}    {variant_description}
    e.g      chr1    100000     100001   strain1        A                     European
    """

    from .variant_profiler.profile import variantProfilerMain

    subprocessID="VARIANT-PROFILER"

    vprint(
        subprocessID,
        "Running variant profiling...",
        "prYellow",
    )

    chdir(args.profilerWDir)

    variant_profile = variantProfilerMain(args, variant_bams)

    chdir(args.runWDir)

    return variant_profile


def makeFinalReport(args, variant_profile, reports):
    """ make final report from the kraken2 and bowtie2 reports
    """

    import json
    from os import listdir

    from Afanc.utilities.getVersions import getVersionsScreen

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
        datadict["Detection_events"] = { "Clustering_results" : [] }

    ## if a variant profile exists then include it in the final report
    if variant_profile != None:
        datadict["Detection_events"]["Variant_profile"] = variant_profile

    ## collect arguments
    jsondict["arguments"] = {
        "database" : args.database,
        "fastqs" : args.fastq,
        "num_threshold" : args.num_threshold,
        "pct_threshold" : args.pct_threshold,
        "upper_bound_threshold" : args.upper_bound,
        "lower_bound_thresold" : args.lower_bound,
        "mapping_sensitivity" : args.mapping_sensitivity,
        "output_prefix" : args.output_prefix,
    }

    ## collect versions
    jsondict["versions"] = {
        "screen" : getVersionsScreen()
    }

    ## collect k2 json reports
    bt2_reportsbox = [g for g in listdir(args.reportsDir) if g.endswith("mapstats.json")]

    ## read general report
    with open(f"{args.reportsDir}/{args.output_prefix}.k2.json", 'r') as k2fin:
        k2data = json.load(k2fin)

        for event in k2data["Detection_events"]:

            ## initialise warnings
            event["warnings"] = []

            ## handle instances where an assembly cannot be found or no_map mode was used
            if "assembly" in event:

                ## block to deal with most likely variants
                if "closest_variant" in event:
                    variant_flag = True
                    taxon_id = str(event["closest_variant"]["taxon_id"])

                    ## in instances where this hit was subjected to variant profiling, the assembly used for mapping will
                    ## belong to the species rather than the closest variant
                    if "assembly" in event["closest_variant"]:
                        assembly = event["closest_variant"]["assembly"]
                    else:
                        assembly = event["assembly"]
                
                ## no variants
                else:
                    variant_flag = False
                    assembly = event["assembly"]
                    taxon_id = str(event["taxon_id"])

            ## block to deal with instances where there is no assembly
            ## this is to ensure that a final json is constructed when no_map mode is used   
            else:
                ## block to deal with most likely variants
                if "closest_variant" in event:
                    variant_flag = True
                    assembly = None
                    taxon_id = str(event["closest_variant"]["taxon_id"])

                ## no variants
                else:
                    variant_flag = False
                    assembly = None
                    taxon_id = str(event["taxon_id"])

            ## block dealing with hits which have an accompanying assembly which reads were mapped to
            if not assembly == None:
                assembly_prefix = path.basename(path.splitext(assembly)[0])
                if assembly_prefix.endswith("_genomic"):
                    assembly_prefix = assembly_prefix.strip("_genomic")

                ## find Bowtie2 report
                ## NOTE: there should only ever be 1 element in this list, since taxon_id should be unique to each taxon
                report_json = [report for report in bt2_reportsbox if report.endswith(assembly_prefix + ".mapstats.json")][0]

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

            datadict["Detection_events"]["Clustering_results"].append(event)

            final_report = path.abspath(f"./{args.output_prefix}.json")

    with open(final_report, "w") as fout:
        json.dump(jsondict, fout, indent = 4)

    return final_report


def cleanOutdir(args, final_report):
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
            print(f"Error: Could not remove {args.bt2WDir}, directory not found.")

        if path.isfile(f"{args.k2WDir}/{args.output_prefix}.k2.txt"):
            remove(f"{args.k2WDir}/{args.output_prefix}.k2.txt")
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove {args.k2WDir}/{args.output_prefix}.k2.txt, file not found.")

        if path.isfile(f"{args.profilerWDir}/{args.output_prefix}.tmp.bam"):
            remove(f"{args.profilerWDir}/{args.output_prefix}.tmp.bam")
            remove(f"{args.profilerWDir}/{args.output_prefix}.tmp.bai")
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove {args.profilerWDir}/{args.output_prefix}.tmp.bam, file not found.")

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

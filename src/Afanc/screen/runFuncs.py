from sys import exit
from shutil import move, rmtree
from os import chdir, path, mkdir, remove

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import vprint
from Afanc.screen.report.final_report import makeFinalReport


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
    makeKronaChart(args)

    ## if no_map is True then exit Afanc screen
    if args.no_map:
        final_report = makeFinalReport(args, None, None)
        vprint("FINISHED", f"no_map mode finished. Metagenomic report can be found in {final_report}\n", "prGreen")
        return 0

    ## parse kraken2 report to a json
    hits_json = getHits(args, out_json)

    ## map to my genomes
    reports, mapped_bams = map2Hits(args)

    ## call SNPs from the competitive per-reference BAMs before lineage classification
    snp_profile = profileSNPs(args, mapped_bams)
    lineage_profile = profileLineages(args, snp_profile, mapped_bams)

    ## generate the final report
    final_report = makeFinalReport(args, reports, snp_profile, lineage_profile)

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

    fq_string = ' '.join(str(x) for x in args.fastq)

    runline = f"kraken2 \
     --db {args.k2_database} \
     --minimum-base-quality 15 \
     --confidence 0.05 \
     --threads {args.threads} \
     --output {args.output_prefix}.k2.txt \
     --report {args.output_prefix}.k2.report.txt \
     --paired {fq_string}"

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
    from .getGenomes import getLocalGenomes

    subprocessID="GET_HITS"

    chdir(args.mappingWDir)

    vprint(
        subprocessID,
        "Getting assemblies from local autodatabase directory...",
        "prYellow",
    )
    getLocalGenomes(out_json, args)

    chdir(args.runWDir)

    return out_json


def map2Hits(args):
    """ takes the hits from the kraken2 module
    """
    from .mapHits import make_accessions_dict, gen_index, run_map

    subprocessID="MAPPING"

    vprint(
        subprocessID,
        "Mapping reads to hits...",
        "prYellow",
    )

    ## index genomes
    chdir(args.mappingWDir)
    mkdir("Hits")

    assemblies_for_mapping = make_accessions_dict(args)

    ## RUN COMPETITIVE BWA MAPPING
    combined_reference = gen_index(args, assemblies_for_mapping)

    reports, mapped_bams = run_map(args, combined_reference, assemblies_for_mapping)

    chdir(args.runWDir)

    return reports, mapped_bams


def profileSNPs(args, mapped_bams):
    """Call SNPs from per-reference competitive mapping BAMs."""
    from .variant_calling.profile import run_snp_profiling

    subprocessID="SNP-PROFILER"

    vprint(
        subprocessID,
        f"Running SNP profiling with {args.variant_caller}...",
        "prYellow",
    )

    chdir(args.profilerWDir)

    snp_profile = run_snp_profiling(args, mapped_bams)

    chdir(args.runWDir)

    return snp_profile


def profileLineages(args, snp_profile, mapped_bams):
    """Classify lineage profiles where database profile models are available."""
    from Afanc.classifier.bayesian_profile import run_lineage_classification

    if args.no_lineage_classify:
        return {
            accession: {
                "accession": accession,
                "status": "not_run",
                "reason": "disabled_by_no_lineage_classify",
                "snp_json": snp_record.get("snp_json"),
            }
            for accession, snp_record in snp_profile.items()
        }
    if len(mapped_bams) > 1 and not args.lineage_profile_compound:
        vprint(
            "LINEAGE-PROFILER",
            "Multiple mapped taxa detected; skipping lineage classification unless --lineage-profile-compound is set.",
            "prYellow",
        )
        return {
            accession: {
                "accession": accession,
                "status": "not_run",
                "reason": "mixed_taxa_requires_lineage_profile_compound",
                "snp_json": snp_record.get("snp_json"),
            }
            for accession, snp_record in snp_profile.items()
        }

    subprocessID="LINEAGE-PROFILER"

    vprint(
        subprocessID,
        "Running Bayesian lineage classification for profiled taxa...",
        "prYellow",
    )

    lineage_profile = run_lineage_classification(args, snp_profile, mapped_bams)

    return lineage_profile


def makeKronaChart(args):
    """Generate a Krona chart from the filtered Kraken-like screen report."""
    from Afanc.utilities.krona import run_krona_from_kraken_report

    subprocessID = "KRONA"
    vprint(
        subprocessID,
        "Making screen Krona chart...",
        "prYellow",
    )

    filtered_report = path.join(args.k2WDir, f"{args.output_prefix}.filtered.k2.report.txt")
    raw_report = path.join(args.k2WDir, f"{args.output_prefix}.k2.report.txt")
    krona_input = filtered_report if path.exists(filtered_report) else raw_report
    output_html = path.join(args.kronaWDir, f"{args.output_prefix}.krona.html")

    run_krona_from_kraken_report(
        krona_input,
        output_html,
        stdout=args.stdout,
        stderr=args.stderr,
        subprocess_id=subprocessID,
    )

    chdir(args.runWDir)


def cleanOutdir(args, final_report):
    """ Cleans the output directory according to provided arguments.

    clean : remove the mapping working directory.
    superclean : move the results JSON to ./ and remove the entire output directory
    """

    if args.clean:
        ## check the mappingWDir exists
        if path.isdir(args.mappingWDir):
            rmtree(args.mappingWDir, ignore_errors=True)
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove {args.mappingWDir}, directory not found.")

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
            krona_report = path.join(args.kronaWDir, f"{args.output_prefix}.krona.html")
            if path.isfile(krona_report):
                move(krona_report, path.join(args.baseRunDir, path.basename(krona_report)))
            rmtree(args.runWDir, ignore_errors=True)
        else:
            ## If it fails, inform the user.
            print(f"Error: Could not remove {args.runWDir}, directory not found.")

        return new_final_report_path

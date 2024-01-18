import sys
import gzip
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import path, listdir, remove
from collections import defaultdict

from Afanc.utilities.runCommands import command
from Afanc.utilities.generalUtils import parseBT2out
from .maths.mappingMetrics import gini, genomeSize, breadthofCoverage, meanDOC, medianDOC


def make_accessions_dict(args, variant_species):
    """ Finds fasta files for mapping, and constructs a dictionary of form

    {
        <accession> : [ <path>, <variantFlag> ]
        ...
    }
    """

    ## collect accession : assembly pairs for mapping
    assemblies_for_mapping = {}

    ## check whether to override assembly mapping with assemblies given in the variant profile
    if args.variant_profile != False:
        for species, assembly in variant_species.items():

            ## skip False species flags
            if not species:
                continue

            ## store accession assembly pairs
            tmp_acc = path.basename(path.splitext(assembly)[0])

            ## strip out genomic tag from the accession
            if tmp_acc.endswith("_genomic"):
                accession = tmp_acc.split("_genomic")[0]
            else:
                accession = tmp_acc

            ## store both the assembly and a flag this dataset for variant profiling
            assemblies_for_mapping[accession] = [assembly, species]

    for assembly in listdir(args.bt2WDir):
        if assembly.endswith(".fna.gz") or assembly.endswith(".fna") or assembly.endswith(".fa.gz") or assembly.endswith(".fa"):
            ## store accession assembly pairs
            tmp_acc = path.basename(path.splitext(assembly)[0])

            ## strip out genomic tag from the accession
            if tmp_acc.endswith("_genomic"):
                accession = tmp_acc.split("_genomic")[0]
            else:
                accession = tmp_acc

            ## store both the assembly and a False flag to skip variant profiling on this dataset
            if accession not in assemblies_for_mapping:
                assemblies_for_mapping[accession] = [assembly, False]

    return assemblies_for_mapping


def modify_fastas(assemblies_for_mapping):

    for accession, [ assembly, variant_flag ] in assemblies_for_mapping.items():
        modify_fasta(assembly, accession)


def modify_fasta(fasta_file, accession):
    """ Takes a fasta file and modifies each chromosome header to include
    the file name for easy grepping out of a SAM file.
    """

    records = []

    if fasta_file.endswith(".gz"):
        fin = gzip.open(fasta_file, "rt")
    else:
        fin = open(fasta_file, "r")

    for rec in SeqIO.parse(fin, "fasta"):
        record = SeqRecord(
            Seq(rec.seq),
            id=f"{rec.id}___{accession}",
            name=rec.id,
            description=rec.description,
        )
        records.append(record)

    SeqIO.write(records, f"{accession}.tmp.fa", "fasta")
    fin.close()


def gen_index(args, assemblies_for_mapping):
    """ Generate the BT2 index from all assembly hits
    """

    assemblies_line = ",".join([f"{acc}.tmp.fa" for acc in assemblies_for_mapping])
    bt2_index = "Hits/Hits"

    buildline = f"bowtie2-build {assemblies_line} {bt2_index}"
    stdout, stderr = command(buildline, "BUILD").run_comm(1, args.stdout, args.stderr)

    return bt2_index


def parse_sam_lines(args, tmp_sam, accession):
    """ takes a set of mapped reads grepped out of a sam file using a grep flag and
    strip out the grep flag. Creates a bam file from mapped reads

    ___{accession} was appended to the end of each chromosome ID prior to mapping to facilitate
    partitioning of reads by mapped assembly.
    """

    read_box = []
    sorted_bam = path.abspath(f"{accession}.sorted.bam")

    ## read in the temporary sam and reformat the chromosome ID
    ## write the new sam to disk
    with open(tmp_sam, 'r') as fin, open(f"{accession}.sam", 'w') as fout:
        for read in fin.readlines():
            new_read = read.replace(f"___{accession}", "")
            fout.write(new_read)

    sortline = f"samtools view -bh {accession}.sam | samtools sort - > {sorted_bam}"
    command(sortline, "MAP").run_comm(0, args.stdout, args.stderr)

    ## remove large sam files
    remove(tmp_sam)
    remove(f"{accession}.sam")

    return sorted_bam


def run_map(args, bt2_index, assemblies_for_mapping):
    """ Map reads to the aggregated index and capture mapping statistics
    """

    ## map reads to the aggregated index
    mapline = f"bowtie2 -p {args.threads} -x {bt2_index} --{args.mapping_sensitivity} -1 {args.fastq[0]} -2 {args.fastq[1]} > Hits.sam"
    stdout, stderr = command(mapline, "MAP").run_comm(1, args.stdout, args.stderr)

    ## capture general mapping stats
    with open(f"./Hits.mapstats.txt", 'w') as fout:
        fout.write(f"{stdout.decode()}\n\n{stderr.decode()}")

    ## store sam files used for variant profiling
    variant_bams = {}
    ## store accession : report_json pairs in a dictionary
    reports = {}

    for accession, [ assembly, variant_flag ] in assemblies_for_mapping.items():

        ## store data in a dict for output as a json
        datadict = defaultdict(dict)

        datadict["warnings"] = {}
        datadict["input_data"]["reference"] = f"{args.bt2WDir}{bt2_index}"
        datadict["input_data"]["fastq_1"] = f"{args.fastq[0]}"
        datadict["input_data"]["fastq_2"] = f"{args.fastq[1]}"

        ## grep out SAM lines for each accession
        # grepline = f"grep '{accession}' Hits.sam | samtools view -bh - | samtools sort - > {sorted_bam}"
        grepline = f"grep '{accession}' Hits.sam > '{accession}'.tmp.sam"
        command(grepline, "MAP").run_comm_quiet(0, args.stdout, args.stderr)

        ## strip out the accession grep flags from the SAM lines and convert to sorted BAM
        mapped_reads = stdout.decode()
        sorted_bam = parse_sam_lines(args, f"{accession}.tmp.sam", accession)

        ## check whether to pass this bam to the variant profiling module and ensure the species name is correctly formatted
        if variant_flag != False:
            variant_bams[accession] = [ sorted_bam, variant_flag.replace(" ", "_") ]

        ## calculate the depth at each covered position
        depthline = f"samtools depth {sorted_bam}"

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

        report_json_out = f"{args.reportsDir}/{accession}.mapstats.json"
        with open(report_json_out, 'w') as fout:
            json.dump(datadict, fout, indent = 4)

        reports[accession] = report_json_out

    return variant_bams, reports

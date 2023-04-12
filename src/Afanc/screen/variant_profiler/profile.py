import json

from os import path, listdir, mkdir

from .bam import Bam
from Afanc.utilities.runCommands import command


def get_species_hits(hits_json):
    """ Reads the hits JSON and finds the species which were identified in the sample.
    """

    with open(hits_json, 'r') as fin:
        json_dict = json.load(fin)

        return [ f["name"].replace(" ", "_") for f in json_dict["Detection_events"] ]


def run_map(args, bam_path, species, ref_fasta, variants_bed):
    """ Index the reference fasta of the species hits and map reads
    """

    ## make reference fasta BT2 index
    mkdir(f"{species}_BT2index")
    bt2_index = f"{species}_BT2index/{species}_BT2index"
    buildline = f"bowtie2-build {ref_fasta} {bt2_index}"
    command(buildline, "MAP").run_comm_quiet(0, args.stdout, args.stderr)

    ## map reads and convert to sorted bam file
    mapline = f"bowtie2 -p 4 -x {bt2_index} --very-sensitive -1 {args.fastq[0]} -2 {args.fastq[1]} | samtools view -bh - | samtools sort - > {args.output_prefix}.sorted.bam"
    command(mapline, "MAP").run_comm_quiet(0, args.stdout, args.stderr)

    ## index bam file
    command(f"samtools index {args.output_prefix}.sorted.bam", "MAP").run_comm_quiet(0, args.stdout, args.stderr)

    return f"{args.output_prefix}.sorted.bam"


def variantProfilerMain(args, variant_bams):
    """ Main function for Afanc variant-profiler
    """
    # species_box = get_species_hits(hits_json)

    variant_profile = {}

    for accession, [ bam_path, species ] in variant_bams.items():
        if species in args.variant_profile:
            print(f"Carrying out variant profiling for detected species : {species} using {bam_path}")

            ## index bam file
            command(f"samtools index {bam_path}", "MAP").run_comm_quiet(0, args.stdout, args.stderr)

            # sorted_bam = run_map(args, bam_path, species, *args.variant_profile[species])

            prefix = path.basename(path.splitext(bam_path)[0])
            bam = Bam(bam_path, prefix, "illumina")
            top_variants = bam.bam2profile(*args.variant_profile[species])

            variant_profile[species] = top_variants

    return variant_profile

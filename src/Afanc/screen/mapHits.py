import gzip
import json
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os import path, listdir, remove
from collections import defaultdict

from Afanc.utilities.runCommands import command
from .mapping.bwa import map_reads_to_bam
from .maths.mappingMetrics import gini, genomeSize, breadthofCoverage, meanDOC, medianDOC


FASTA_EXTENSIONS = (".fna.gz", ".fna", ".fa.gz", ".fa", ".fasta.gz", ".fasta")


def make_accessions_dict(args):
    """ Finds fasta files for mapping, and constructs a dictionary of form

    {
        <accession> : [ <path>, <variantFlag> ]
        ...
    }
    """

    ## collect accession : assembly pairs for mapping
    assemblies_for_mapping = {}

    for assembly in listdir(args.mappingWDir):
        if assembly.endswith(FASTA_EXTENSIONS):
            ## store accession assembly pairs
            tmp_acc = path.basename(path.splitext(assembly)[0])

            ## strip out genomic tag from the accession
            if tmp_acc.endswith("_genomic"):
                accession = tmp_acc.split("_genomic")[0]
            else:
                accession = tmp_acc

            ## store accession assembly pairs
            if accession not in assemblies_for_mapping:
                assemblies_for_mapping[accession] = assembly

    return assemblies_for_mapping


def build_combined_reference(assemblies_for_mapping, output_fasta="Hits/Hits.combined.fa"):
    """Build a single competitive reference FASTA for BWA mapping."""
    with open(output_fasta, "w") as fout:
        for accession, assembly in assemblies_for_mapping.items():
            if assembly.endswith(".gz"):
                fin = gzip.open(assembly, "rt")
            else:
                fin = open(assembly, "r")

            with fin:
                for rec in SeqIO.parse(fin, "fasta"):
                    record = SeqRecord(
                        rec.seq,
                        id=f"{rec.id}___{accession}",
                        name=f"{rec.id}___{accession}",
                        description=f"{rec.description}___{accession}",
                    )
                    SeqIO.write(record, fout, "fasta")

    return path.abspath(output_fasta)


def gen_index(args, assemblies_for_mapping):
    """Build the combined BWA reference used for competitive mapping."""
    return build_combined_reference(assemblies_for_mapping)


def _strip_accession_suffix(reference_name, accession):
    suffix = f"___{accession}"
    if reference_name.endswith(suffix):
        return reference_name[:-len(suffix)]
    return reference_name


def partition_bam_by_accession(args, combined_bam, assemblies_for_mapping):
    """Split a combined competitive BAM into per-accession sorted BAMs."""
    output_bams = {}

    with pysam.AlignmentFile(combined_bam, "rb") as bam_in:
        source_header = bam_in.header.to_dict()

        for accession in assemblies_for_mapping:
            suffix = f"___{accession}"
            source_sq = source_header.get("SQ", [])
            selected_refs = [
                (idx, ref)
                for idx, ref in enumerate(source_sq)
                if ref["SN"].endswith(suffix)
            ]

            output_sam = path.abspath(f"{accession}.sam")
            sorted_bam = path.abspath(f"{accession}.sorted.bam")
            ref_id_map = {}
            target_sq = []

            for target_idx, (source_idx, ref) in enumerate(selected_refs):
                ref_id_map[source_idx] = target_idx
                new_ref = dict(ref)
                new_ref["SN"] = _strip_accession_suffix(ref["SN"], accession)
                target_sq.append(new_ref)

            target_header = dict(source_header)
            target_header["SQ"] = target_sq

            with pysam.AlignmentFile(output_sam, "w", header=target_header) as bam_out:
                for read in bam_in.fetch(until_eof=True):
                    if read.reference_id not in ref_id_map:
                        continue

                    read.reference_id = ref_id_map[read.reference_id]
                    if read.next_reference_id in ref_id_map:
                        read.next_reference_id = ref_id_map[read.next_reference_id]
                    else:
                        read.next_reference_id = -1
                        read.next_reference_start = -1

                    bam_out.write(read)

                bam_in.reset()

            sortline = f"samtools view -bh {output_sam} | samtools sort - > {sorted_bam}"
            command(sortline, "MAP").run_comm(0, args.stdout, args.stderr)
            command(f"samtools index {sorted_bam}", "MAP").run_comm(0, args.stdout, args.stderr)
            remove(output_sam)

            output_bams[accession] = sorted_bam

    return output_bams


def _write_mapping_stats(args, accession, assembly, sorted_bam, combined_reference):
    datadict = defaultdict(dict)

    datadict["warnings"] = {}
    datadict["input_data"]["reference"] = combined_reference
    datadict["input_data"]["fastq_1"] = f"{args.fastq[0]}"
    datadict["input_data"]["fastq_2"] = f"{args.fastq[1]}"

    depthline = f"samtools depth {sorted_bam}"

    depthstdout, depthstderr = command(depthline, "DEPTH").run_comm(1, None, args.stderr)
    covarray = [i.split('\t') for i in depthstdout.decode().split('\n')]

    if len(covarray) <= 1:
        datadict["warnings"]["no_unique_map"] = f"No reads map uniquely to assembly {assembly}."

    else:
        genomesize = genomeSize(assembly)

        meandoc = meanDOC(covarray)
        mediandoc = medianDOC(covarray)
        boc = breadthofCoverage(covarray, genomesize)

        gini_co = gini(covarray)

        datadict["map_data"]["mean_DOC"] = meandoc
        datadict["map_data"]["median_DOC"] = mediandoc
        datadict["map_data"]["proportion_cov"] = boc
        datadict["map_data"]["gini"] = gini_co

        if boc < 0.05:
            datadict["warnings"]["FP-warning"] = f"Coverage across {assembly} low (<5%). Result could be false-positive."
        if boc < 0.01:
            datadict["warnings"]["FP-warning"] = f"Coverage across {assembly} very low (<1%). Result likely to be false-positive."

    report_json_out = f"{args.reportsDir}/{accession}.mapstats.json"
    with open(report_json_out, 'w') as fout:
        json.dump(datadict, fout, indent = 4, default=str)

    return report_json_out


def run_map(args, combined_reference, assemblies_for_mapping):
    """Map reads competitively with BWA and capture per-reference statistics."""

    combined_bam = path.abspath("Hits/Hits.combined.sorted.bam")
    map_reads_to_bam(
        ref_fasta=combined_reference,
        r1_fastq=args.fastq[0],
        r2_fastq=args.fastq[1],
        output_bam=combined_bam,
        sample_name=args.output_prefix,
        cpus=args.threads,
        tmpdir=path.abspath("Hits"),
    )

    ## capture general mapping stats
    with open(f"./Hits.mapstats.txt", 'w') as fout:
        fout.write(f"competitive_mapper\tbwa\ncombined_reference\t{combined_reference}\ncombined_bam\t{combined_bam}\n")

    ## store accession : report_json pairs in a dictionary
    reports = {}
    mapped_bams = {}
    accession_bams = partition_bam_by_accession(args, combined_bam, assemblies_for_mapping)

    for accession, assembly in assemblies_for_mapping.items():

        sorted_bam = accession_bams[accession]
        mapped_bams[accession] = {
            "bam": sorted_bam,
            "assembly": path.abspath(assembly),
            "lineage_profile": getattr(args, "lineage_profiles_by_accession", {}).get(accession),
        }

        report_json_out = _write_mapping_stats(args, accession, assembly, sorted_bam, combined_reference)
        reports[accession] = report_json_out

    return reports, mapped_bams

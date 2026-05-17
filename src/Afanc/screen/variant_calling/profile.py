from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Set

from Bio import SeqIO

from Afanc.utilities.runCommands import command
from .callers import run_variant_caller


DEFAULT_ALLELE_ID_FORMAT = "{chrom}.{start}.{ref}.{alt}"
DNA4 = {"A", "C", "G", "T"}


def run_snp_profiling(
    args: Any,
    mapped_bams: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Call SNPs for each competitively mapped reference and write SNP JSON."""
    snp_profiles = {}

    for accession, mapping_record in mapped_bams.items():
        bam_path = Path(mapping_record["bam"])
        source_reference = Path(mapping_record["assembly"])
        reference_fasta = prepare_variant_reference(
            source_reference=source_reference,
            accession=accession,
            output_dir=Path(args.profilerWDir),
            stdout=args.stdout,
            stderr=args.stderr,
        )
        output_prefix = Path(args.profilerWDir) / accession

        caller_kwargs = {}
        if args.variant_caller == "freebayes":
            caller_kwargs["cpus"] = args.threads

        caller_outputs = run_variant_caller(
            bam_path=bam_path,
            reference_fasta=reference_fasta,
            output_prefix=output_prefix,
            caller_name=args.variant_caller,
            **caller_kwargs,
        )

        vcf_for_json = caller_outputs.get("filtered_vcf") or caller_outputs["raw_vcf"]
        depth_path = output_prefix.with_suffix(".depth.bed")
        write_depth_bed(
            bam_path=bam_path,
            output_depth=depth_path,
            stdout=args.stdout,
            stderr=args.stderr,
        )
        snp_json_path = Path(args.reportsDir) / f"{accession}.snps.json"
        snp_json = vcf_to_snp_json(
            vcf_path=vcf_for_json,
            depth_path=depth_path,
            output_json=snp_json_path,
            min_qual=args.snp_min_qual,
            allow_filtered=args.snp_allow_filtered,
            accept_missing_qual=args.snp_accept_missing_qual,
            min_dp=args.snp_min_dp,
            min_missing_depth=args.snp_min_missing_depth,
        )

        snp_profiles[accession] = {
            "accession": accession,
            "assembly": str(source_reference),
            "bam": str(bam_path),
            "reference_fasta": str(reference_fasta),
            "caller": args.variant_caller,
            "raw_vcf": str(caller_outputs["raw_vcf"]),
            "filtered_vcf": str(caller_outputs.get("filtered_vcf")) if caller_outputs.get("filtered_vcf") else None,
            "depth_bed": str(depth_path),
            "snp_json": str(snp_json_path),
            "snp_count": len(snp_json["alleles"]),
            "missing_count": len(snp_json["missing"]),
            "allele_id_format": DEFAULT_ALLELE_ID_FORMAT,
        }

    return snp_profiles


def prepare_variant_reference(
    source_reference: Path,
    accession: str,
    output_dir: Path,
    stdout: Optional[Any] = None,
    stderr: Optional[Any] = None,
) -> Path:
    """Write an uncompressed per-accession reference and ensure it is indexed."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_reference = output_dir / f"{accession}.variant_ref.fa"

    opener = gzip.open if str(source_reference).endswith(".gz") else open
    with opener(source_reference, "rt") as fin, output_reference.open("w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            SeqIO.write(record, fout, "fasta")

    fai_path = Path(str(output_reference) + ".fai")
    if not fai_path.is_file():
        command(f"samtools faidx {output_reference}", "VARCALL").run_comm_quiet(0, stdout, stderr)

    return output_reference


def write_depth_bed(
    bam_path: Path,
    output_depth: Path,
    stdout: Optional[Any] = None,
    stderr: Optional[Any] = None,
) -> Path:
    command(f"samtools depth -aa {bam_path} > {output_depth}", "VARCALL").run_comm_quiet(0, stdout, stderr)
    return output_depth


def vcf_to_snp_json(
    vcf_path: Path,
    depth_path: Path,
    output_json: Path,
    min_qual: float = 30.0,
    allow_filtered: bool = False,
    accept_missing_qual: bool = False,
    min_dp: Optional[int] = None,
    min_missing_depth: int = 10,
) -> Dict[str, Any]:
    """Convert a VCF and depth BED into Ardal-compatible SNP JSON."""
    allele_map = sample_alleles_from_vcf(
        vcf_path=vcf_path,
        min_qual=min_qual,
        allow_filtered=allow_filtered,
        accept_missing_qual=accept_missing_qual,
        min_dp=min_dp,
    )
    alleles = set()
    for sample_alleles in allele_map.values():
        alleles.update(sample_alleles)

    snp_json = {
        "alleles": sorted(alleles),
        "missing": missing_positions_from_bed(depth_path, min_missing_depth),
    }

    with output_json.open("w") as fout:
        json.dump(snp_json, fout, indent=2)

    return snp_json


def sample_alleles_from_vcf(
    vcf_path: Path,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: Optional[int],
) -> Dict[str, Set[str]]:
    sample_ids = None
    allele_map = {}

    with open_maybe_gzip(vcf_path) as fin:
        for line in fin:
            if not line:
                continue
            if line[0] == "#":
                if line.startswith("#CHROM"):
                    header_parts = line.rstrip("\n").split("\t")
                    sample_ids = header_parts[9:] if len(header_parts) > 9 else ["sample"]
                    allele_map = {sample: set() for sample in sample_ids}
                continue

            if sample_ids is None:
                raise ValueError(f"Missing #CHROM header in {vcf_path}")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _id, ref, alt_field, qual, filt, info = parts[:8]
            if len(ref) != 1 or ref.upper() not in DNA4:
                continue
            if not passes_filters(
                qual_field=qual,
                filter_field=filt,
                min_qual=min_qual,
                allow_filtered=allow_filtered,
                accept_missing_qual=accept_missing_qual,
                min_dp=min_dp,
                info_field=info,
            ):
                continue

            alts = alt_field.split(",")
            if len(parts) < 10:
                for alt in alts:
                    if is_simple_snp(ref, alt):
                        allele_map.setdefault("sample", set()).add(f"{chrom}.{pos}.{ref.upper()}.{alt.upper()}")
                continue

            fmt = parts[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                continue

            for sample_offset, sample_id in enumerate(sample_ids):
                sample_col_index = 9 + sample_offset
                if sample_col_index >= len(parts):
                    continue
                sample_fields = parts[sample_col_index].split(":")
                if len(sample_fields) <= gt_index:
                    continue
                gt_value = sample_fields[gt_index]
                if not gt_value or gt_value == ".":
                    continue

                for allele_token in re.split(r"[|/]", gt_value):
                    if not allele_token or allele_token == ".":
                        continue
                    try:
                        allele_index = int(allele_token)
                    except ValueError:
                        continue
                    if allele_index <= 0:
                        continue
                    alt_index = allele_index - 1
                    if alt_index >= len(alts):
                        continue
                    alt = alts[alt_index]
                    if is_simple_snp(ref, alt):
                        allele_map[sample_id].add(f"{chrom}.{pos}.{ref.upper()}.{alt.upper()}")

    return allele_map


def missing_positions_from_bed(depth_path: Path, min_depth: int) -> List[List[str]]:
    missing = []
    seen = set()

    with open_maybe_gzip(depth_path) as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, pos, depth_str = parts[:3]
            try:
                depth = int(float(depth_str))
            except ValueError:
                continue
            if depth >= min_depth:
                continue
            key = f"{chrom}.{pos}"
            if key in seen:
                continue
            seen.add(key)
            missing.append([chrom, pos])

    return missing


def open_maybe_gzip(input_path: Path):
    return gzip.open(input_path, "rt") if input_path.suffix == ".gz" else input_path.open("r")


def is_simple_snp(ref: str, alt: str) -> bool:
    ref = ref.upper()
    alt = alt.upper()
    return len(ref) == 1 and len(alt) == 1 and ref in DNA4 and alt in DNA4


def passes_filters(
    qual_field: str,
    filter_field: str,
    min_qual: float,
    allow_filtered: bool,
    accept_missing_qual: bool,
    min_dp: Optional[int],
    info_field: str,
) -> bool:
    if qual_field == ".":
        if not accept_missing_qual:
            return False
    else:
        try:
            if float(qual_field) < min_qual:
                return False
        except ValueError:
            return False

    if not allow_filtered and filter_field not in ("PASS", "."):
        return False

    if min_dp is not None:
        dp = parse_info_dp(info_field)
        if dp is None or dp < min_dp:
            return False

    return True


def parse_info_dp(info_field: str) -> Optional[int]:
    for item in info_field.split(";"):
        if item.startswith("DP="):
            try:
                return int(item[3:])
            except ValueError:
                return None
    return None

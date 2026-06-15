from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple

from Bio import SeqIO


DNA4 = {"A", "C", "G", "T"}


def read_reference_sequences(reference_fasta: Path) -> Dict[str, str]:
    opener = gzip.open if str(reference_fasta).endswith(".gz") else open
    sequences = {}
    with opener(reference_fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[str(record.id)] = str(record.seq).upper()
    return sequences


def reference_lengths(reference_fasta: Path) -> Dict[str, int]:
    return {
        chrom: len(sequence)
        for chrom, sequence in read_reference_sequences(reference_fasta).items()
    }


def parse_vcf_contigs(vcf_path: Path) -> Dict[str, int]:
    contigs = {}
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                break
            if not line.startswith("##contig=<"):
                continue
            inner = line.strip()[len("##contig=<"):-1]
            fields = {}
            for item in inner.split(","):
                if "=" not in item:
                    continue
                key, value = item.split("=", 1)
                fields[key] = value
            if "ID" not in fields or "length" not in fields:
                continue
            contigs[fields["ID"]] = int(fields["length"])
    return contigs


def validate_vcf_reference(vcf_path: Path, reference_fasta: Path) -> Dict[str, object]:
    reference_sequences = read_reference_sequences(reference_fasta)
    reference = {
        chrom: len(sequence)
        for chrom, sequence in reference_sequences.items()
    }
    contigs = parse_vcf_contigs(vcf_path)
    validation = {
        "valid": True,
        "errors": [],
        "checked_contigs": 0,
        "checked_records": 0,
        "vcf": str(vcf_path),
        "reference": str(reference_fasta),
    }

    if not contigs:
        validation["valid"] = False
        validation["errors"].append(
            "VCF has no ##contig headers with lengths; reference compatibility cannot be checked."
        )
        return validation

    for chrom, length in contigs.items():
        if chrom not in reference:
            validation["valid"] = False
            validation["errors"].append(f"VCF contig {chrom!r} is not present in the profile reference.")
            continue
        validation["checked_contigs"] += 1
        if reference[chrom] != length:
            validation["valid"] = False
            validation["errors"].append(
                f"VCF contig {chrom!r} length is {length}, but profile reference length is {reference[chrom]}."
            )

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos_str, _id, ref, _alt = parts[:5]
            ref = ref.upper()
            if len(ref) != 1 or ref not in DNA4:
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                validation["valid"] = False
                validation["errors"].append(f"VCF record has invalid position {pos_str!r}.")
                continue
            if chrom not in reference_sequences:
                validation["valid"] = False
                validation["errors"].append(f"VCF record uses contig {chrom!r}, which is not in the profile reference.")
                continue
            if pos < 1 or pos > len(reference_sequences[chrom]):
                validation["valid"] = False
                validation["errors"].append(f"VCF record position {chrom}:{pos} is outside the profile reference.")
                continue
            validation["checked_records"] += 1
            observed = reference_sequences[chrom][pos - 1].upper()
            if observed != ref:
                validation["valid"] = False
                validation["errors"].append(
                    f"VCF REF allele at {chrom}:{pos} is {ref}, but profile reference has {observed}."
                )

    return validation


def parse_allele_id_format(allele_id_format: str) -> list[Tuple[str | None, str]]:
    tokens = []
    pattern = re.compile(r"\{([^{}]+)\}")
    cursor = 0
    for match in pattern.finditer(allele_id_format):
        literal = allele_id_format[cursor:match.start()]
        if literal:
            tokens.append((None, re.escape(literal)))
        tokens.append((match.group(1), r"(?P<%s>.+?)" % match.group(1)))
        cursor = match.end()
    literal = allele_id_format[cursor:]
    if literal:
        tokens.append((None, re.escape(literal)))
    return tokens


def allele_format_fields(allele_id_format: str) -> set[str]:
    return {
        field
        for field, _pattern in parse_allele_id_format(allele_id_format)
        if field is not None
    }


def parse_allele_id(allele_id: str, allele_id_format: str) -> Dict[str, str]:
    pattern = "".join(pattern for _field, pattern in parse_allele_id_format(allele_id_format))
    match = re.fullmatch(pattern, str(allele_id))
    if not match:
        raise ValueError(f"Allele ID {allele_id!r} does not match format {allele_id_format!r}.")
    return match.groupdict()


def normalise_position(parsed: Mapping[str, str]) -> int:
    if "pos" in parsed:
        return int(parsed["pos"])
    if "start" in parsed:
        return int(parsed["start"])
    raise ValueError("Allele ID format must include 'pos' or 'start'.")


def validate_allele_reference(
    allele_ids: Iterable[str],
    allele_id_format: str,
    reference_fasta: Path,
) -> Dict[str, object]:
    fields = allele_format_fields(allele_id_format)
    validation = {
        "valid": True,
        "errors": [],
        "checked_alleles": 0,
        "reference": str(reference_fasta),
        "allele_id_format": allele_id_format,
    }

    required = {"chrom", "ref", "alt"}
    missing_required = sorted(required - fields)
    if missing_required:
        validation["valid"] = False
        validation["errors"].append(
            "Allele ID format is missing required fields: " + ", ".join(missing_required)
        )
    if "pos" not in fields and "start" not in fields:
        validation["valid"] = False
        validation["errors"].append("Allele ID format must include 'pos' or 'start'.")
    if "ref" not in fields:
        validation["valid"] = False
        validation["errors"].append(
            "Allele ID format must include 'ref' so profile reference correctness can be checked."
        )
    if not validation["valid"]:
        return validation

    reference = read_reference_sequences(reference_fasta)
    for allele_id in allele_ids:
        try:
            parsed = parse_allele_id(str(allele_id), allele_id_format)
            chrom = str(parsed["chrom"])
            pos = normalise_position(parsed)
            ref = str(parsed["ref"]).upper()
        except Exception as exc:
            validation["valid"] = False
            validation["errors"].append(str(exc))
            continue

        if len(ref) != 1 or ref not in DNA4:
            validation["valid"] = False
            validation["errors"].append(f"Allele ID {allele_id!r} has unsupported ref allele {ref!r}.")
            continue
        if chrom not in reference:
            validation["valid"] = False
            validation["errors"].append(f"Allele ID {allele_id!r} uses missing contig {chrom!r}.")
            continue
        if pos < 1 or pos > len(reference[chrom]):
            validation["valid"] = False
            validation["errors"].append(f"Allele ID {allele_id!r} position {pos} is outside contig {chrom}.")
            continue
        validation["checked_alleles"] += 1
        observed = reference[chrom][pos - 1].upper()
        if observed != ref:
            validation["valid"] = False
            validation["errors"].append(
                f"Allele ID {allele_id!r} expects ref {ref} at {chrom}:{pos}, but profile reference has {observed}."
            )

    return validation

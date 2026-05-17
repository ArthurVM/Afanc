#!/usr/bin/env python3
"""Convert a canonical SNP BED/barcode file into a JSON lineage model.

The input is expected to be a BED-like tab-delimited file where the
first 5 columns are required:

    chrom  start  end  lineage_ID  alt_allele

Optional legacy annotation columns 6-8 are also accepted and preserved
when present:

    clade_label  spoligotype  rd

The resulting JSON model is intentionally species-level and sparse:
each locus names the lineage it is intended to define plus a shared
background frequency for all non-target lineages.

If a reference FASTA is supplied, reference alleles are resolved and
added to the model. If not, the model still carries chrom/pos/alt and
sets ref to null.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

DEFAULT_ALLELE_ID_FORMAT = "{chrom}.{start}.{ref}.{alt}"
HIERARCHY_ROOT_NODE = "__root__"

try:
    from ardal.core.ArdalHeaderUtils import ArdalHeaderUtils
except ImportError:  # pragma: no cover
    ArdalHeaderUtils = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_bed",
        type=Path,
        help="Canonical SNP BED/barcode file.",
    )
    parser.add_argument(
        "-o",
        "--output-json",
        type=Path,
        default=None,
        help="Output JSON path. Defaults to <input>.canonical_model.json.",
    )
    parser.add_argument(
        "--species-id",
        default="Mycobacterium_tuberculosis",
        help="Species identifier to store in the model.",
    )
    parser.add_argument(
        "--model-id",
        default=None,
        help="Model identifier. Defaults to the output stem.",
    )
    parser.add_argument(
        "--reference-name",
        default=None,
        help="Optional reference name to store in the model.",
    )
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        default=None,
        help="Optional FASTA used to resolve reference alleles.",
    )
    parser.add_argument(
        "--target-frequency",
        type=float,
        default=0.99,
        help="Target-lineage allele frequency for canonical SNPs.",
    )
    parser.add_argument(
        "--background-frequency",
        type=float,
        default=0.01,
        help="Background allele frequency for non-target lineages.",
    )
    parser.add_argument(
        "--allele-id-format",
        default=DEFAULT_ALLELE_ID_FORMAT,
        help="Format string used when the first column already contains an encoded SNP ID.",
    )
    return parser.parse_args()


def natural_sort_key(value: str) -> List[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def load_reference_sequences(reference_fasta: Optional[Path]) -> Dict[str, str]:
    if reference_fasta is None:
        return {}

    sequences: Dict[str, List[str]] = {}
    current_name: Optional[str] = None

    with open_text(reference_fasta) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                if current_name in sequences:
                    raise ValueError(f"Duplicate contig in FASTA: {current_name}")
                sequences[current_name] = []
                continue
            if current_name is None:
                raise ValueError("Invalid FASTA: sequence data encountered before a header.")
            sequences[current_name].append(line.upper())

    return {name: "".join(parts) for name, parts in sequences.items()}


def lineage_parent(lineage_id: str, known_lineages: Iterable[str]) -> Optional[str]:
    known = set(known_lineages)
    parts = lineage_id.split(".")
    for size in range(len(parts) - 1, 0, -1):
        candidate = ".".join(parts[:size])
        if candidate in known:
            return candidate
    return None


def lineage_ancestors(lineage_id: str, known_lineages: Iterable[str]) -> List[str]:
    known = set(known_lineages)
    parts = lineage_id.split(".")
    ancestors = []
    for size in range(1, len(parts) + 1):
        candidate = ".".join(parts[:size])
        if candidate in known:
            ancestors.append(candidate)
    return ancestors


def validate_frequency(name: str, value: float) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be between 0.0 and 1.0, got {value}.")


def read_barcode_rows(
    input_bed: Path,
    reference_sequences: Dict[str, str],
    allele_id_format: str = DEFAULT_ALLELE_ID_FORMAT,
) -> List[dict]:
    rows = []
    with input_bed.open("r") as handle:
        raw_rows = []
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            raw_rows.append((line_number, line.split("\t")))

    encoded_position_map = decode_allele_positions_with_ardal(
        [fields[0] for _, fields in raw_rows],
        allele_id_format=allele_id_format,
    )
    allele_id_pattern = compile_allele_id_pattern(allele_id_format)

    for line_number, fields in raw_rows:
            if len(fields) < 4:
                raise ValueError(
                    f"{input_bed}:{line_number} expected at least 4 tab-separated columns, got {len(fields)}."
                )

            decoded_allele = try_decode_allele_id(
                fields[0],
                pattern=allele_id_pattern,
                ardal_position_map=encoded_position_map,
            )
            parsed_bed_row = None
            bed_parse_error = None
            if len(fields) >= 5:
                try:
                    parsed_bed_row = parse_bed_barcode_fields(fields, input_bed, line_number)
                except ValueError as exc:
                    bed_parse_error = exc

            if decoded_allele is not None and (
                parsed_bed_row is None or encoded_allele_disagrees_with_bed(decoded_allele, parsed_bed_row)
            ):
                chrom, pos_1based, ref, alt = decoded_allele
                lineage_id = fields[3]
                parsed_row = {
                    "chrom": chrom,
                    "start": pos_1based - 1,
                    "end": pos_1based,
                    "pos": pos_1based,
                    "lineage_id": lineage_id,
                    "ref": ref,
                    "alt": alt.upper(),
                    "clade_label": None,
                    "spoligotype": None,
                    "region_of_difference": None,
                    "source_line_number": line_number,
                }
            elif parsed_bed_row is not None:
                parsed_row = parsed_bed_row
            elif bed_parse_error is not None:
                raise bed_parse_error
            else:
                raise ValueError(
                    f"{input_bed}:{line_number} expected either BED-style canonical SNP columns or an encoded SNP ID in column 1."
                )

            parsed_row["ref"] = resolve_reference_allele(
                input_bed=input_bed,
                line_number=line_number,
                reference_sequences=reference_sequences,
                chrom=parsed_row["chrom"],
                pos_1based=parsed_row["pos"],
                expected_ref=parsed_row["ref"],
            )
            rows.append(parsed_row)

    return rows


def parse_bed_barcode_fields(fields: List[str], input_bed: Path, line_number: int) -> dict:
    if len(fields) < 5:
        raise ValueError(f"{input_bed}:{line_number} expected at least 5 tab-separated columns, got {len(fields)}.")

    chrom, start, end, lineage_id, alt = fields[:5]
    start_i = int(start)
    end_i = int(end)
    if end_i != start_i + 1:
        raise ValueError(f"{input_bed}:{line_number} expected a 1-bp interval, got start={start_i}, end={end_i}.")

    clade_label = fields[5] if len(fields) > 5 else None
    spoligotype = fields[6] if len(fields) > 6 else None
    rd = fields[7] if len(fields) > 7 else None
    return {
        "chrom": chrom,
        "start": start_i,
        "end": end_i,
        "pos": end_i,
        "lineage_id": lineage_id,
        "ref": None,
        "alt": alt.upper(),
        "clade_label": clade_label,
        "spoligotype": spoligotype,
        "region_of_difference": rd,
        "source_line_number": line_number,
    }


def compile_allele_id_pattern(allele_id_format: str) -> re.Pattern:
    placeholders = re.findall(r"\{([^}]+)\}", allele_id_format)
    accepted_placeholders = {"chr", "chrom", "start", "end", "ref", "alt"}
    invalid_placeholders = sorted(set(placeholders) - accepted_placeholders)
    if invalid_placeholders:
        raise ValueError(
            f"Invalid allele_id_format placeholder(s) {invalid_placeholders}. Accepted placeholders are {sorted(accepted_placeholders)}."
        )

    if "chr" in placeholders and "chrom" in placeholders:
        raise ValueError("allele_id_format cannot contain both {chr} and {chrom}.")
    if "start" not in placeholders or "alt" not in placeholders:
        raise ValueError("allele_id_format must contain at least {start} and {alt}.")
    if "chr" not in placeholders and "chrom" not in placeholders:
        raise ValueError("allele_id_format must contain either {chrom} or {chr}.")

    pattern = re.escape(allele_id_format)
    replacements = {
        "chr": r"(?P<chrom>.+)",
        "chrom": r"(?P<chrom>.+)",
        "start": r"(?P<start>\d+)",
        "end": r"(?P<end>\d+)",
        "ref": r"(?P<ref>.+)",
        "alt": r"(?P<alt>.+)",
    }
    for key, replacement in replacements.items():
        pattern = pattern.replace(re.escape(f"{{{key}}}"), replacement)
    return re.compile(f"^{pattern}$")


def decode_allele_positions_with_ardal(
    allele_ids: List[str],
    allele_id_format: str,
) -> Dict[str, tuple[str, int]]:
    if ArdalHeaderUtils is None:
        return {}

    unique_alleles = list(dict.fromkeys(allele_ids))
    header_utils = ArdalHeaderUtils(
        headers={"guids": [], "alleles": unique_alleles},
        meta={},
        allele_id_format=normalise_ardal_allele_id_format(allele_id_format),
    )
    allele_positions = header_utils.get_allele_positions()
    return {
        allele_id: (str(chrom), int(pos))
        for allele_id, (chrom, pos) in allele_positions.items()
    }


def normalise_ardal_allele_id_format(allele_id_format: str) -> str:
    return allele_id_format.replace("{chrom}", "{chr}")


def try_decode_allele_id(
    allele_id: str,
    pattern: re.Pattern,
    ardal_position_map: Dict[str, tuple[str, int]],
) -> Optional[tuple[str, int, Optional[str], str]]:
    match = pattern.match(allele_id)
    if match is None:
        return None

    parts = match.groupdict()
    alt = parts.get("alt")
    if alt is None:
        return None

    chrom_pos = ardal_position_map.get(allele_id)
    if chrom_pos is not None:
        chrom, pos = chrom_pos
    else:
        chrom = parts.get("chrom")
        if chrom is None:
            return None
        pos = int(parts["start"])

    return chrom, pos, parts.get("ref"), alt


def encoded_allele_disagrees_with_bed(decoded_allele: tuple[str, int, Optional[str], str], parsed_bed_row: dict) -> bool:
    chrom, pos_1based, _, alt = decoded_allele
    return (
        chrom != parsed_bed_row["chrom"]
        or pos_1based != parsed_bed_row["pos"]
        or alt.upper() != parsed_bed_row["alt"]
    )


def resolve_reference_allele(
    input_bed: Path,
    line_number: int,
    reference_sequences: Dict[str, str],
    chrom: str,
    pos_1based: int,
    expected_ref: Optional[str],
) -> Optional[str]:
    if not reference_sequences:
        return expected_ref.upper() if expected_ref is not None else None

    if chrom not in reference_sequences:
        raise ValueError(f"{input_bed}:{line_number} contig {chrom!r} not found in reference FASTA.")

    sequence = reference_sequences[chrom]
    if pos_1based < 1 or pos_1based > len(sequence):
        raise ValueError(f"{input_bed}:{line_number} position {pos_1based} out of bounds for contig {chrom}.")

    ref = sequence[pos_1based - 1]
    if expected_ref is not None and expected_ref.upper() != ref:
        raise ValueError(
            f"{input_bed}:{line_number} encoded reference allele {expected_ref!r} does not match FASTA reference {ref!r}."
        )
    return ref


def build_model(args: argparse.Namespace) -> dict:
    validate_frequency("target_frequency", args.target_frequency)
    validate_frequency("background_frequency", args.background_frequency)

    if args.background_frequency >= args.target_frequency:
        raise ValueError("background_frequency should be lower than target_frequency for canonical SNPs.")

    reference_sequences = load_reference_sequences(args.reference_fasta)
    rows = read_barcode_rows(
        args.input_bed,
        reference_sequences,
        allele_id_format=getattr(args, "allele_id_format", DEFAULT_ALLELE_ID_FORMAT),
    )

    ## collect the direct marker counts before expanding any inherited hierarchy
    direct_marker_counts = Counter(row["lineage_id"] for row in rows)
    lineage_ids = sorted(direct_marker_counts, key=natural_sort_key)
    lineage_set = set(lineage_ids)

    children = defaultdict(list)
    for lineage_id in lineage_ids:
        parent = lineage_parent(lineage_id, lineage_set)
        if parent is not None:
            children[parent].append(lineage_id)

    direct_locus_ids_by_lineage = defaultdict(list)
    loci = []
    for row in rows:
        ref_token = row["ref"] if row["ref"] is not None else "N"
        locus_id = f"{row['chrom']}:{row['pos']}:{ref_token}>{row['alt']}:{row['lineage_id']}"
        loci.append(
            {
                "locus_id": locus_id,
                "chrom": row["chrom"],
                "pos": row["pos"],
                "bed_start": row["start"],
                "bed_end": row["end"],
                "ref": row["ref"],
                "alt": row["alt"],
                "target_lineage": row["lineage_id"],
                "empirical_bayes": {
                    "target_frequency": args.target_frequency,
                    "background_frequency": args.background_frequency,
                },
                "annotations": {
                    "clade_label": row["clade_label"],
                    "spoligotype": row["spoligotype"],
                    "region_of_difference": row["region_of_difference"],
                    "source_line_number": row["source_line_number"],
                },
            }
        )
        direct_locus_ids_by_lineage[row["lineage_id"]].append(locus_id)

    ## materialise the hierarchy within the model so later code does not need to
    ## infer lineage relationships from string prefixes
    lineage_records = []
    hierarchy_edges = []
    root_lineages = []
    parent_map = {}
    for lineage_id in lineage_ids:
        parent = lineage_parent(lineage_id, lineage_set)
        ancestors = lineage_ancestors(lineage_id, lineage_set)
        direct_children = sorted(children.get(lineage_id, []), key=natural_sort_key)
        direct_locus_ids = sorted(direct_locus_ids_by_lineage.get(lineage_id, []), key=natural_sort_key)
        inherited_locus_ids = []
        for ancestor in ancestors:
            inherited_locus_ids.extend(direct_locus_ids_by_lineage.get(ancestor, []))
        inherited_locus_ids = sorted(inherited_locus_ids, key=natural_sort_key)

        parent_map[lineage_id] = parent
        if parent is None:
            root_lineages.append(lineage_id)
        else:
            hierarchy_edges.append({"parent": parent, "child": lineage_id})

        lineage_records.append(
            {
                "lineage_id": lineage_id,
                "prior": round(1.0 / len(lineage_ids), 12),
                "parent_lineage": parent,
                "is_root": parent is None,
                "ancestor_lineages": ancestors[:-1],
                "direct_children": direct_children,
                "direct_locus_ids": direct_locus_ids,
                "inherited_locus_ids": inherited_locus_ids,
                "direct_marker_count": len(direct_locus_ids),
                "inherited_marker_count": len(inherited_locus_ids),
            }
        )

    children_by_parent = {
        lineage_id: sorted(children.get(lineage_id, []), key=natural_sort_key)
        for lineage_id in lineage_ids
        if children.get(lineage_id)
    }
    leaf_lineages = sorted(
        [lineage_id for lineage_id in lineage_ids if not children.get(lineage_id)],
        key=natural_sort_key,
    )
    hierarchy_is_hierarchical = any(parent is not None for parent in parent_map.values())
    serialised_parent_map = {lineage_id: parent_map[lineage_id] for lineage_id in lineage_ids}
    serialised_children = dict(children_by_parent)
    serialised_roots = sorted(root_lineages, key=natural_sort_key)
    top_classifier_node = None
    if hierarchy_is_hierarchical:
        if len(serialised_roots) > 1:
            serialised_parent_map[HIERARCHY_ROOT_NODE] = None
            serialised_children[HIERARCHY_ROOT_NODE] = list(serialised_roots)
            for root in serialised_roots:
                serialised_parent_map[root] = HIERARCHY_ROOT_NODE
            serialised_roots = [HIERARCHY_ROOT_NODE]
            top_classifier_node = HIERARCHY_ROOT_NODE
        elif len(serialised_roots) == 1:
            top_classifier_node = serialised_roots[0]

    output_json = args.output_json
    if output_json is None:
        output_json = args.input_bed.with_suffix("").with_suffix(".canonical_model.json")

    model_id = args.model_id or output_json.stem
    created_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    model = {
        "schema_version": "1.2",
        "model_id": model_id,
        "species_id": args.species_id,
        "reference": {
            "name": args.reference_name,
            "path": str(args.reference_fasta) if args.reference_fasta else None,
            "contigs": sorted(reference_sequences) if reference_sequences else None,
            "ref_alleles_resolved": bool(reference_sequences),
        },
        "architecture": {
            "feature_family": "canonical",
            "topology": "hierarchical" if hierarchy_is_hierarchical else "flat",
            "classification_workflow": (
                "canonical_hierarchical" if hierarchy_is_hierarchical else "canonical_flat"
            ),
            "emission_model": "none",
        },
        "model_type": "canonical_snp",
        "hierarchy": {
            "is_hierarchical": hierarchy_is_hierarchical,
            "parent_map": serialised_parent_map if hierarchy_is_hierarchical else None,
            "children_by_parent": serialised_children if hierarchy_is_hierarchical else {},
            "root_lineages": serialised_roots if hierarchy_is_hierarchical else sorted(root_lineages, key=natural_sort_key),
            "leaf_lineages": leaf_lineages,
            "top_classifier_node": top_classifier_node,
            "edges": sorted(hierarchy_edges, key=lambda edge: (natural_sort_key(edge["parent"]), natural_sort_key(edge["child"]))),
            "lineage_naming_rule": "dot-delimited prefix, when parent exists in the model",
            "note": (
                "Lineage relationships and direct/inherited locus memberships are materialized in the model. "
                "Classifiers should use these declared relationships instead of inferring hierarchy ad hoc."
            ),
        },
        "lineages": lineage_records,
        "loci": loci,
        "empirical_bayes": {
            "emission_family": "binomial",
            "default_target_frequency": args.target_frequency,
            "default_background_frequency": args.background_frequency,
            "snp_encoding": "sparse_target_plus_shared_background",
        },
        "full_bayes": {
            "status": "not_fit",
            "model_family": "beta_binomial_hierarchical",
            "priors": {
                "mu_prior": None,
                "kappa_prior": None,
            },
            "posterior_summary": None,
            "draw_store": None,
        },
        "reporting_support": {
            "status": "not_fit",
            "method": None,
            "evidence_metric": "usable_loci",
            "base_thresholds": None,
            "collision_map": None,
            "thresholds": None,
            "fallback_policy": None,
            "provenance": None,
        },
        "summary": {
            "lineage_count": len(lineage_records),
            "locus_count": len(loci),
            "has_resolved_reference_alleles": bool(reference_sequences),
        },
        "provenance": {
            "created_at": created_at,
            "created_by": Path(__file__).name,
            "input_file": str(args.input_bed),
            "input_format": "canonical_barcode_bed_5plus",
        },
    }

    return model, output_json


def main() -> None:
    args = parse_args()
    model, output_json = build_model(args)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as handle:
        json.dump(model, handle, indent=2)
    print(f"Wrote canonical lineage model to {output_json}")


if __name__ == "__main__":
    main()

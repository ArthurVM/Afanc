from __future__ import annotations

import json
import math
import re
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

import pysam


ModelInput = Union[str, Path, Mapping[str, Any]]
SnpJsonInput = Union[str, Path, Mapping[str, Any], Sequence[Any]]
DEFAULT_ALLELE_ID_FORMAT = "{chrom}.{start}.{ref}.{alt}"
SNP_JSON_LIST_KEYS = ("snps", "SNPs", "alleles", "variants")
SNP_JSON_ITEM_KEYS = ("snp", "name", "allele", "allele_id", "variant")
HIERARCHICAL_MODEL_TYPE = "hierarchical_empirical_geolineage"
CANONICAL_MODEL_TYPE = "canonical_snp"
HIERARCHY_ROOT_NODE = "__root__"


@dataclass(frozen=True)
class VariantEvidence:
    locus_id: str
    chrom: str
    pos: int
    ref: Optional[str]
    alt: str
    target_lineage: str
    observed: bool
    usable: bool
    quality: Optional[float]
    filters: Tuple[str, ...]
    depth: Optional[int]
    ref_depth: Optional[int]
    alt_depth: Optional[int]
    allele_frequency: Optional[float]
    genotype: Optional[Tuple[Optional[int], ...]]
    excluded_reasons: Tuple[str, ...]


@dataclass(frozen=True)
class ParsedVcfEvidence:
    sample_name: Optional[str]
    evidence_by_locus: Dict[str, VariantEvidence]
    observed_locus_ids: Tuple[str, ...]
    usable_locus_ids: Tuple[str, ...]


ParsedEvidence = ParsedVcfEvidence


def validate_vcf_qc(
    vcf_path: Union[str, Path],
    expected_caller: Optional[str] = None,
    require_filtered_vcf: bool = False,
    expected_min_quality: Optional[int] = None,
    expected_min_depth: Optional[int] = None,
    expected_min_alt_fraction: Optional[float] = None,
    require_hom_alt: Optional[bool] = None,
) -> Dict[str, Any]:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        header_text = str(vcf.header)
        format_fields = set(vcf.header.formats.keys())

    header_lower = header_text.lower()

    detected_caller = _detect_caller_from_header(header_lower)
    if expected_caller is not None and detected_caller != expected_caller:
        raise ValueError(
            f"VCF caller validation failed for {vcf_path}: expected {expected_caller!r}, got {detected_caller!r}."
        )

    ## ensure the basic depth/allele fields required for parsing are present
    if detected_caller == "freebayes":
        required_formats = {"GT", "DP"}
        if not ({"AD"} <= format_fields or {"AO", "RO"} <= format_fields):
            raise ValueError(
                f"VCF QC validation failed for {vcf_path}: expected AD or AO/RO format fields for FreeBayes-style parsing."
            )
        missing_formats = required_formats - format_fields
    elif detected_caller == "bcftools":
        required_formats = {"GT", "DP", "AD"}
        missing_formats = required_formats - format_fields
    else:
        required_formats = {"GT", "DP"}
        missing_formats = required_formats - format_fields

    if missing_formats:
        raise ValueError(
            f"VCF QC validation failed for {vcf_path}: missing required FORMAT fields {sorted(missing_formats)}."
        )

    ## check that the expected bcftools filter stuff is still present in the header
    if require_filtered_vcf:
        if "bcftools_viewcommand" not in header_lower or "bcftools_normcommand" not in header_lower:
            raise ValueError(
                f"VCF QC validation failed for {vcf_path}: expected bcftools view/norm command lines in the VCF header."
            )

        if require_hom_alt:
            _require_header_substring(header_text, 'FMT/GT="1/1"', vcf_path)
        if expected_min_quality is not None:
            _require_header_substring(header_text, f"QUAL>={int(expected_min_quality)}", vcf_path)
        if expected_min_depth is not None:
            _require_header_substring(header_text, f"FMT/DP>={int(expected_min_depth)}", vcf_path)
        if expected_min_alt_fraction is not None:
            _require_any_header_substring(
                header_text,
                [
                    f"(FMT/AO)/(FMT/DP)>={expected_min_alt_fraction}",
                    f"(FMT/AD[1])/(FMT/DP)>={expected_min_alt_fraction}",
                    f"FORMAT/AD[1]/FORMAT/DP>={expected_min_alt_fraction}",
                ],
                vcf_path,
            )

    return {
        "vcf_path": Path(vcf_path),
        "detected_caller": detected_caller,
        "format_fields": sorted(format_fields),
        "is_filtered_vcf": "bcftools_viewcommand" in header_lower and "bcftools_normcommand" in header_lower,
    }


def load_lineage_model(model: ModelInput) -> Dict[str, Any]:
    ## allow either a loaded model dictionary or a path to a json model
    if isinstance(model, Mapping):
        return dict(model)

    model_path = Path(model)
    with model_path.open("r") as handle:
        return json.load(handle)


def write_model_targets_bed(model: ModelInput, output_bed: Union[str, Path]) -> Path:
    model_dict = _prepare_model_for_classification(load_lineage_model(model))
    output_path = Path(output_bed)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ## deduplicate intervals before writing the targeting bed
    with output_path.open("w") as handle:
        for locus in _iter_unique_bed_loci(model_dict):
            handle.write(
                "\t".join(
                    [
                        str(locus["chrom"]),
                        str(locus["bed_start"]),
                        str(locus["bed_end"]),
                        str(locus["locus_id"]),
                    ]
                )
            )
            handle.write("\n")

    return output_path


def parse_vcf_for_classification(
    vcf_path: Union[str, Path],
    model: ModelInput,
    sample_name: Optional[str] = None,
    min_depth: int = 0,
    min_quality: Optional[float] = None,
    require_pass: bool = False,
) -> ParsedVcfEvidence:
    model_dict = _prepare_model_for_classification(load_lineage_model(model))
    loci_by_key = _build_model_locus_lookup(model_dict)
    evidence_by_locus = _initialise_evidence(model_dict, missing_reason="not_in_vcf")

    with pysam.VariantFile(str(vcf_path)) as vcf:
        resolved_sample = _resolve_sample_name(vcf, sample_name)

        ## capture only loci which appear in the model
        ## everything else missing
        for record in vcf:
            if not record.alts:
                continue

            filter_keys = tuple(record.filter.keys())
            for alt_index, alt in enumerate(record.alts):
                loci = _match_model_loci(
                    loci_by_key,
                    record.chrom,
                    record.pos,
                    record.ref,
                    alt,
                )
                if not loci:
                    continue

                for locus in loci:
                    evidence = _record_to_evidence(
                        record=record,
                        locus=locus,
                        alt_index=alt_index,
                        sample_name=resolved_sample,
                        min_depth=min_depth,
                        min_quality=min_quality,
                        require_pass=require_pass,
                        filter_keys=filter_keys,
                    )
                    evidence_by_locus[locus["locus_id"]] = evidence

    observed_locus_ids = tuple(
        locus_id for locus_id, evidence in evidence_by_locus.items() if evidence.observed
    )
    usable_locus_ids = tuple(
        locus_id for locus_id, evidence in evidence_by_locus.items() if evidence.usable
    )

    return ParsedVcfEvidence(
        sample_name=resolved_sample,
        evidence_by_locus=evidence_by_locus,
        observed_locus_ids=observed_locus_ids,
        usable_locus_ids=usable_locus_ids,
    )


def parse_snp_json_for_classification(
    snp_json: SnpJsonInput,
    model: ModelInput,
    allele_id_format: str = DEFAULT_ALLELE_ID_FORMAT,
    sample_name: Optional[str] = None,
) -> ParsedEvidence:
    model_dict = _prepare_model_for_classification(load_lineage_model(model))
    loci_by_key = _build_model_locus_lookup(model_dict)
    evidence_by_locus = _initialise_evidence(model_dict, missing_reason="not_in_json")

    loaded_sample_name, snp_names = _load_snp_json_snps(snp_json)
    resolved_sample = sample_name if sample_name is not None else loaded_sample_name

    observed_locus_ids = []
    usable_locus_ids = []
    seen_locus_ids = set()

    pattern = _compile_allele_id_pattern(allele_id_format)

    for snp_name in snp_names:
        chrom, pos, _, ref, alt = _decode_allele_id(
            allele_id=snp_name,
            pattern=pattern,
            allele_id_format=allele_id_format,
        )
        loci = _match_model_loci(
            loci_by_key,
            chrom,
            pos,
            ref,
            alt,
            allow_refless_input=True,
        )
        if not loci:
            continue

        for locus in loci:
            locus_id = locus["locus_id"]
            if locus_id in seen_locus_ids:
                continue
            seen_locus_ids.add(locus_id)

            evidence_by_locus[locus_id] = VariantEvidence(
                locus_id=locus_id,
                chrom=locus["chrom"],
                pos=int(locus["pos"]),
                ref=locus.get("ref"),
                alt=locus["alt"],
                target_lineage=locus["target_lineage"],
                observed=True,
                usable=True,
                quality=None,
                filters=(),
                depth=1,
                ref_depth=0,
                alt_depth=1,
                allele_frequency=1.0,
                genotype=None,
                excluded_reasons=(),
            )
            observed_locus_ids.append(locus_id)
            usable_locus_ids.append(locus_id)

    return ParsedEvidence(
        sample_name=resolved_sample,
        evidence_by_locus=evidence_by_locus,
        observed_locus_ids=tuple(observed_locus_ids),
        usable_locus_ids=tuple(usable_locus_ids),
    )


def classify_evidence_against_model(
    parsed_evidence: ParsedVcfEvidence,
    model: ModelInput,
    mode: str = "empirical_bayes",
    lineage_scope: str = "inherited",
    kappa_weighting: bool = False,
    min_kappa_weight: float = 0.25,
    full_weight_kappa: Optional[float] = 5.0,
) -> Dict[str, Any]:
    model_dict = _prepare_model_for_classification(load_lineage_model(model))
    architecture = model_dict["architecture"]
    feature_family = architecture["feature_family"]
    topology = architecture["topology"]

    if feature_family == "empirical" and topology == "hierarchical":
        return _classify_hierarchical_evidence_against_model(
            parsed_evidence=parsed_evidence,
            model=model_dict,
            mode=mode,
            lineage_scope=lineage_scope,
            kappa_weighting=kappa_weighting,
            min_kappa_weight=min_kappa_weight,
            full_weight_kappa=full_weight_kappa,
        )

    if feature_family == "empirical" and topology == "flat":
        return _classify_flat_evidence_against_model(
            parsed_evidence=parsed_evidence,
            model=model_dict,
            mode=mode,
            lineage_scope=lineage_scope,
            kappa_weighting=kappa_weighting,
            min_kappa_weight=min_kappa_weight,
            full_weight_kappa=full_weight_kappa,
        )

    if feature_family == "canonical" and topology == "hierarchical":
        return _classify_hierarchical_canonical_evidence_against_model(
            parsed_evidence=parsed_evidence,
            model=model_dict,
            lineage_scope=lineage_scope,
        )

    if feature_family == "canonical" and topology == "flat":
        return _classify_flat_canonical_evidence_against_model(
            parsed_evidence=parsed_evidence,
            model=model_dict,
            lineage_scope=lineage_scope,
        )

    raise ValueError(
        f"Unsupported model architecture: feature_family={feature_family!r}, topology={topology!r}."
    )


def classify_vcf_against_model(
    vcf_path: Union[str, Path],
    model: ModelInput,
    mode: str = "empirical_bayes",
    lineage_scope: str = "inherited",
    sample_name: Optional[str] = None,
    min_depth: int = 0,
    min_quality: Optional[float] = None,
    require_pass: bool = False,
    validate_qc: bool = False,
    expected_caller: Optional[str] = None,
    require_filtered_vcf: bool = False,
    expected_min_quality: Optional[int] = None,
    expected_min_depth: Optional[int] = None,
    expected_min_alt_fraction: Optional[float] = None,
    require_hom_alt: Optional[bool] = None,
    kappa_weighting: bool = False,
    min_kappa_weight: float = 0.25,
    full_weight_kappa: Optional[float] = 5.0,
) -> Dict[str, Any]:
    prepared_model = _prepare_model_for_classification(load_lineage_model(model))
    qc_summary = None
    if validate_qc:
        qc_summary = validate_vcf_qc(
            vcf_path=vcf_path,
            expected_caller=expected_caller,
            require_filtered_vcf=require_filtered_vcf,
            expected_min_quality=expected_min_quality,
            expected_min_depth=expected_min_depth,
            expected_min_alt_fraction=expected_min_alt_fraction,
            require_hom_alt=require_hom_alt,
        )

    parsed_evidence = parse_vcf_for_classification(
        vcf_path=vcf_path,
        model=prepared_model,
        sample_name=sample_name,
        min_depth=min_depth,
        min_quality=min_quality,
        require_pass=require_pass,
    )
    result = classify_evidence_against_model(
        parsed_evidence=parsed_evidence,
        model=prepared_model,
        mode=mode,
        lineage_scope=lineage_scope,
        kappa_weighting=kappa_weighting,
        min_kappa_weight=min_kappa_weight,
        full_weight_kappa=full_weight_kappa,
    )
    if qc_summary is not None:
        result["qc_validation"] = qc_summary
    return result


def classify_snp_json_against_model(
    snp_json: SnpJsonInput,
    model: ModelInput,
    mode: str = "empirical_bayes",
    lineage_scope: str = "inherited",
    allele_id_format: str = DEFAULT_ALLELE_ID_FORMAT,
    sample_name: Optional[str] = None,
    kappa_weighting: bool = False,
    min_kappa_weight: float = 0.25,
    full_weight_kappa: Optional[float] = 5.0,
) -> Dict[str, Any]:
    prepared_model = _prepare_model_for_classification(load_lineage_model(model))
    parsed_evidence = parse_snp_json_for_classification(
        snp_json=snp_json,
        model=prepared_model,
        allele_id_format=allele_id_format,
        sample_name=sample_name,
    )
    
    # print(f"""mode={mode} ; kappa_weighting={kappa_weighting} ; min_kappa_weight={min_kappa_weight} ; full_weight_kappa={full_weight_kappa}""")
    
    result = classify_evidence_against_model(
        parsed_evidence=parsed_evidence,
        model=prepared_model,
        mode=mode,
        lineage_scope=lineage_scope,
        kappa_weighting=kappa_weighting,
        min_kappa_weight=min_kappa_weight,
        full_weight_kappa=full_weight_kappa,
    )
    result["allele_id_format"] = allele_id_format
    return result


def _build_model_locus_lookup(model: Mapping[str, Any]) -> Dict[str, Dict[Tuple[Any, ...], list[Dict[str, Any]]]]:
    exact_lookup: Dict[Tuple[str, int, Optional[str], str], list[Dict[str, Any]]] = defaultdict(list)
    model_refless_lookup: Dict[Tuple[str, int, str], list[Dict[str, Any]]] = defaultdict(list)
    any_ref_grouped: Dict[Tuple[str, int, str], list[Dict[str, Any]]] = defaultdict(list)
    any_ref_lookup: Dict[Tuple[str, int, str], list[Dict[str, Any]]] = {}

    for locus in model["loci"]:
        chrom = str(locus["chrom"])
        pos = int(locus["pos"])
        ref = locus.get("ref")
        alt = str(locus["alt"])
        exact_lookup[(chrom, pos, ref, alt)].append(locus)

        refless_key = (chrom, pos, alt)
        if ref is None:
            model_refless_lookup[refless_key].append(locus)
        any_ref_grouped[refless_key].append(locus)

    for refless_key, loci in any_ref_grouped.items():
        non_null_refs = {locus.get("ref") for locus in loci if locus.get("ref") is not None}
        if len(non_null_refs) <= 1:
            any_ref_lookup[refless_key] = loci

    return {
        "exact": dict(exact_lookup),
        "model_refless": dict(model_refless_lookup),
        "any_ref_unique": any_ref_lookup,
    }


def _initialise_evidence(model: Mapping[str, Any], missing_reason: str) -> Dict[str, VariantEvidence]:
    evidence = {}
    for locus in model["loci"]:
        ## initialise every model locus as absent from the current input
        ## overwrite only when a matching record is observed
        evidence[locus["locus_id"]] = VariantEvidence(
            locus_id=locus["locus_id"],
            chrom=locus["chrom"],
            pos=int(locus["pos"]),
            ref=locus.get("ref"),
            alt=locus["alt"],
            target_lineage=locus["target_lineage"],
            observed=False,
            usable=False,
            quality=None,
            filters=(),
            depth=None,
            ref_depth=None,
            alt_depth=None,
            allele_frequency=None,
            genotype=None,
            excluded_reasons=(missing_reason,),
        )
    return evidence


def _resolve_sample_name(vcf: pysam.VariantFile, sample_name: Optional[str]) -> Optional[str]:
    samples = list(vcf.header.samples)
    if sample_name is not None:
        if sample_name not in samples:
            raise ValueError(f"Sample {sample_name!r} was not found in VCF.")
        return sample_name
    if not samples:
        return None
    return samples[0]


def _match_model_loci(
    lookup: Mapping[str, Mapping[Tuple[Any, ...], Sequence[Dict[str, Any]]]],
    chrom: str,
    pos: int,
    ref: Optional[str],
    alt: str,
    allow_refless_input: bool = False,
) -> list[Dict[str, Any]]:
    exact_key = (chrom, pos, ref, alt)
    if exact_key in lookup["exact"]:
        return list(lookup["exact"][exact_key])

    fallback_key = (chrom, pos, alt)
    if ref is None and allow_refless_input:
        return list(lookup["any_ref_unique"].get(fallback_key, ()))

    return list(lookup["model_refless"].get(fallback_key, ()))


def _record_to_evidence(
    record: pysam.VariantRecord,
    locus: Mapping[str, Any],
    alt_index: int,
    sample_name: Optional[str],
    min_depth: int,
    min_quality: Optional[float],
    require_pass: bool,
    filter_keys: Tuple[str, ...],
) -> VariantEvidence:
    sample_data = record.samples[sample_name] if sample_name is not None else None
    depth = _extract_depth(sample_data, record)
    ref_depth, alt_depth = _extract_allele_depths(sample_data, alt_index)
    allele_frequency = _extract_allele_frequency(sample_data, record, alt_index, depth, alt_depth)
    genotype = tuple(sample_data.get("GT")) if sample_data is not None and sample_data.get("GT") else None

    excluded_reasons = []
    if depth is None:
        excluded_reasons.append("missing_depth")
    elif depth < min_depth:
        excluded_reasons.append("below_min_depth")

    if min_quality is not None and record.qual is not None and record.qual < min_quality:
        excluded_reasons.append("below_min_quality")

    if require_pass and filter_keys and filter_keys != ("PASS",):
        excluded_reasons.append("filter_not_pass")

    ## classifier expects depth aware evidence so loci missing DP/AD are treated as observed but unusable
    usable = not excluded_reasons and depth is not None and alt_depth is not None

    return VariantEvidence(
        locus_id=locus["locus_id"],
        chrom=locus["chrom"],
        pos=int(locus["pos"]),
        ref=locus.get("ref"),
        alt=locus["alt"],
        target_lineage=locus["target_lineage"],
        observed=True,
        usable=usable,
        quality=record.qual,
        filters=filter_keys,
        depth=depth,
        ref_depth=ref_depth,
        alt_depth=alt_depth,
        allele_frequency=allele_frequency,
        genotype=genotype,
        excluded_reasons=tuple(excluded_reasons),
    )


def _load_snp_json_snps(snp_json: SnpJsonInput) -> Tuple[Optional[str], Tuple[str, ...]]:
    if isinstance(snp_json, Mapping):
        data = dict(snp_json)
    elif isinstance(snp_json, Sequence) and not isinstance(snp_json, (str, bytes, Path)):
        data = snp_json
    else:
        json_path = Path(snp_json)
        with json_path.open("r") as handle:
            data = json.load(handle)

    sample_name = data.get("sample_name") if isinstance(data, Mapping) else None
    raw_snps = _extract_raw_snp_list(data)
    snp_names = []
    for item in raw_snps:
        snp_names.append(_coerce_snp_name(item))
    return sample_name, tuple(snp_names)


def _extract_raw_snp_list(data: Any) -> Sequence[Any]:
    if isinstance(data, Sequence) and not isinstance(data, (str, bytes, Path)):
        return data

    if not isinstance(data, Mapping):
        raise ValueError("SNP JSON must be either a list of SNP names or an object containing one.")

    for key in SNP_JSON_LIST_KEYS:
        if key in data:
            value = data[key]
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes, Path)):
                return value
            raise ValueError(f"SNP JSON field {key!r} must contain a list of SNP names.")

    list_fields = [
        key
        for key, value in data.items()
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, Path))
    ]
    if len(list_fields) == 1:
        return data[list_fields[0]]

    raise ValueError(
        "SNP JSON must be a list or contain exactly one SNP list field such as 'snps', 'SNPs', 'alleles', or 'variants'."
    )


def _coerce_snp_name(item: Any) -> str:
    if isinstance(item, str):
        return item

    if isinstance(item, Mapping):
        for key in SNP_JSON_ITEM_KEYS:
            value = item.get(key)
            if isinstance(value, str):
                return value

    raise ValueError("Each SNP JSON entry must be a string or an object with a SNP name field.")


def _compile_allele_id_pattern(allele_id_format: str) -> re.Pattern[str]:
    if not allele_id_format:
        raise ValueError("allele_id_format cannot be empty.")

    placeholders = re.findall(r"\{([^}]+)\}", allele_id_format)
    accepted_placeholders = {"chr", "chrom", "start", "end", "ref", "alt"}
    invalid_placeholders = sorted(set(placeholders) - accepted_placeholders)
    if invalid_placeholders:
        raise ValueError(
            f"Invalid allele_id_format placeholder(s) {invalid_placeholders}. "
            f"Accepted placeholders are {sorted(accepted_placeholders)}."
        )

    if "chr" in placeholders and "chrom" in placeholders:
        raise ValueError("allele_id_format cannot contain both {chr} and {chrom}; use one or the other.")
    if "start" not in placeholders:
        raise ValueError("allele_id_format must contain a {start} placeholder.")
    if "alt" not in placeholders:
        raise ValueError("allele_id_format must contain an {alt} placeholder.")
    if "chr" not in placeholders and "chrom" not in placeholders:
        raise ValueError("allele_id_format must contain either a {chrom} or {chr} placeholder.")

    pattern = re.escape(allele_id_format)
    placeholder_patterns = {
        "chr": r"(?P<chrom>.+)",
        "chrom": r"(?P<chrom>.+)",
        "start": r"(?P<start>\d+)",
        "end": r"(?P<end>\d+)",
        "ref": r"(?P<ref>.+)",
        "alt": r"(?P<alt>.+)",
    }

    found_placeholders = 0
    for key, regex_pattern in placeholder_patterns.items():
        escaped_placeholder = re.escape(f"{{{key}}}")
        if escaped_placeholder in pattern:
            pattern = pattern.replace(escaped_placeholder, regex_pattern)
            found_placeholders += 1

    if found_placeholders == 0:
        raise ValueError("allele_id_format does not contain any valid placeholders.")

    return re.compile(f"^{pattern}$")


def _decode_allele_id(
    allele_id: str,
    pattern: re.Pattern[str],
    allele_id_format: str,
) -> Tuple[str, int, Optional[int], Optional[str], str]:
    match = pattern.match(allele_id)
    if not match:
        raise ValueError(f"Allele ID {allele_id!r} does not match the format {allele_id_format!r}.")

    parts = match.groupdict()
    chrom = parts.get("chrom")
    alt = parts.get("alt")
    if chrom is None or alt is None:
        raise ValueError(
            f"Allele ID {allele_id!r} could not be decoded into chrom/start/alt fields using {allele_id_format!r}."
        )

    end = int(parts["end"]) if parts.get("end") is not None else None
    return chrom, int(parts["start"]), end, parts.get("ref"), alt


def _extract_depth(sample_data: Any, record: pysam.VariantRecord) -> Optional[int]:
    if sample_data is not None and sample_data.get("DP") is not None:
        return int(sample_data.get("DP"))
    try:
        info_dp = record.info.get("DP")
    except ValueError:
        info_dp = None
    if info_dp is not None:
        return int(info_dp)
    return None


def _extract_allele_depths(sample_data: Any, alt_index: int) -> Tuple[Optional[int], Optional[int]]:
    if sample_data is None:
        return None, None

    ## prefer AD where available
    ## accept freebayes style RO/AO fields too
    ad = sample_data.get("AD")
    if ad is None or len(ad) <= alt_index + 1:
        ro = sample_data.get("RO")
        ao = sample_data.get("AO")
        if ro is None or ao is None:
            return None, None

        if isinstance(ao, Iterable) and not isinstance(ao, (str, bytes)):
            ao = tuple(ao)
            if len(ao) <= alt_index:
                return None, None
            return int(ro), int(ao[alt_index])

        if alt_index != 0:
            return None, None
        return int(ro), int(ao)

    return int(ad[0]), int(ad[alt_index + 1])


def _extract_allele_frequency(
    sample_data: Any,
    record: pysam.VariantRecord,
    alt_index: int,
    depth: Optional[int],
    alt_depth: Optional[int],
) -> Optional[float]:
    if sample_data is not None and sample_data.get("AF") is not None:
        af_value = sample_data.get("AF")
        if isinstance(af_value, Iterable) and not isinstance(af_value, (str, bytes)):
            return float(tuple(af_value)[alt_index])
        return float(af_value)

    try:
        info_af = record.info.get("AF")
    except ValueError:
        info_af = None
    if info_af is not None:
        if isinstance(info_af, Iterable) and not isinstance(info_af, (str, bytes)):
            return float(tuple(info_af)[alt_index])
        return float(info_af)

    if depth and alt_depth is not None:
        return alt_depth / depth

    return None


def _detect_caller_from_header(header_lower: str) -> Optional[str]:
    if "freebayes" in header_lower or "source=freebayes" in header_lower:
        return "freebayes"
    if (
        "source=bcftools" in header_lower
        or "bcftools_callcommand" in header_lower
        or "bcftools_mpileupcommand" in header_lower
        or "bcftools call" in header_lower
        or "bcftools mpileup" in header_lower
    ):
        return "bcftools"
    return None


def _require_header_substring(header_text: str, expected_substring: str, vcf_path: Union[str, Path]) -> None:
    if expected_substring not in header_text:
        raise ValueError(
            f"VCF QC validation failed for {vcf_path}: expected to find {expected_substring!r} in the VCF header."
        )


def _require_any_header_substring(
    header_text: str,
    expected_substrings: Sequence[str],
    vcf_path: Union[str, Path],
) -> None:
    for expected_substring in expected_substrings:
        if expected_substring in header_text:
            return

    raise ValueError(
        f"VCF QC validation failed for {vcf_path}: expected one of {list(expected_substrings)!r} in the VCF header."
    )


def _lineage_scope_field(lineage_scope: str) -> str:
    ## lineage membership is driven from the model hierarchy rather than inferred from lineage name prefixes at classification time
    ## TODO: this has caused issues in the past
    if lineage_scope == "direct":
        return "direct_locus_ids"
    if lineage_scope == "inherited":
        return "inherited_locus_ids"
    raise ValueError("lineage_scope must be either 'direct' or 'inherited'.")


def _normalise_model_architecture(model: Mapping[str, Any]) -> Dict[str, Any]:
    architecture = dict(model.get("architecture") or {})

    if "feature_family" not in architecture:
        if model.get("model_type") == CANONICAL_MODEL_TYPE:
            architecture["feature_family"] = "canonical"
        else:
            architecture["feature_family"] = "empirical"

    if "topology" not in architecture:
        if architecture["feature_family"] == "empirical":
            is_hierarchical = _is_legacy_empirical_hierarchical_model(model)
        else:
            hierarchy = model.get("hierarchy") or {}
            is_hierarchical = bool(hierarchy.get("is_hierarchical"))
            if not is_hierarchical:
                is_hierarchical = any(lineage.get("parent_lineage") is not None for lineage in model.get("lineages", []))
        architecture["topology"] = "hierarchical" if is_hierarchical else "flat"

    if "classification_workflow" not in architecture:
        architecture["classification_workflow"] = (
            f"{architecture['feature_family']}_{architecture['topology']}"
        )

    if "emission_model" not in architecture:
        if architecture["feature_family"] == "canonical":
            architecture["emission_model"] = "none"
        else:
            architecture["emission_model"] = (
                (model.get("provenance") or {}).get("loso", {}).get("emission_model")
                or "empirical_bayes"
            )

    return architecture


def _build_parent_map_from_lineages(model: Mapping[str, Any]) -> Dict[str, Optional[str]]:
    parent_map: Dict[str, Optional[str]] = {}
    for lineage in model.get("lineages", []):
        lineage_id = str(lineage["lineage_id"])
        parent = lineage.get("parent_lineage")
        parent_map[lineage_id] = None if parent is None else str(parent)
    return parent_map


def _build_parent_map_from_edges(hierarchy: Mapping[str, Any]) -> Dict[str, Optional[str]]:
    parent_map: Dict[str, Optional[str]] = {}
    for edge in hierarchy.get("edges", []):
        parent = str(edge["parent"])
        child = str(edge["child"])
        parent_map.setdefault(parent, None)
        parent_map[child] = parent
    return parent_map


def _build_children_by_parent(parent_map: Mapping[str, Optional[str]]) -> Dict[str, list[str]]:
    children_by_parent: Dict[str, list[str]] = defaultdict(list)
    for child, parent in parent_map.items():
        if parent is None:
            continue
        children_by_parent[str(parent)].append(str(child))
    return {
        str(parent): sorted(children)
        for parent, children in children_by_parent.items()
    }


def _normalise_hierarchy(model: Mapping[str, Any], architecture: Mapping[str, Any]) -> Dict[str, Any]:
    hierarchy = dict(model.get("hierarchy") or {})
    topology = architecture["topology"]

    if topology == "flat":
        root_lineages = sorted(
            [
                str(lineage["lineage_id"])
                for lineage in model.get("lineages", [])
                if lineage.get("parent_lineage") is None
            ]
            or list(hierarchy.get("root_lineages") or ()),
        )
        return {
            **hierarchy,
            "is_hierarchical": False,
            "parent_map": None,
            "children_by_parent": {},
            "root_lineages": root_lineages,
            "leaf_lineages": sorted(str(lineage["lineage_id"]) for lineage in model.get("lineages", [])),
            "top_classifier_node": None,
        }

    parent_map = hierarchy.get("parent_map")
    if parent_map:
        normalised_parent_map = {
            str(child): (None if parent is None else str(parent))
            for child, parent in parent_map.items()
        }
    else:
        normalised_parent_map = _build_parent_map_from_lineages(model)
        if not normalised_parent_map:
            normalised_parent_map = _build_parent_map_from_edges(hierarchy)

    children_by_parent = hierarchy.get("children_by_parent")
    if children_by_parent:
        normalised_children = {
            str(parent): sorted(str(child) for child in children)
            for parent, children in children_by_parent.items()
        }
    else:
        normalised_children = _build_children_by_parent(normalised_parent_map)

    roots = sorted(
        hierarchy.get("root_lineages")
        or [node for node, parent in normalised_parent_map.items() if parent is None]
    )
    leaf_lineages = sorted(
        hierarchy.get("leaf_lineages")
        or [
            str(lineage["lineage_id"])
            for lineage in model.get("lineages", [])
            if not lineage.get("direct_children")
        ]
    )
    top_classifier_node = hierarchy.get("top_classifier_node")

    if len(roots) > 1 and top_classifier_node is None:
        normalised_parent_map = dict(normalised_parent_map)
        normalised_parent_map[HIERARCHY_ROOT_NODE] = None
        normalised_children = dict(normalised_children)
        normalised_children[HIERARCHY_ROOT_NODE] = list(roots)
        for root in roots:
            normalised_parent_map[str(root)] = HIERARCHY_ROOT_NODE
        roots = [HIERARCHY_ROOT_NODE]
        top_classifier_node = HIERARCHY_ROOT_NODE
    elif len(roots) == 1 and top_classifier_node is None:
        top_classifier_node = roots[0]

    return {
        **hierarchy,
        "is_hierarchical": True,
        "parent_map": normalised_parent_map,
        "children_by_parent": normalised_children,
        "root_lineages": roots,
        "leaf_lineages": leaf_lineages,
        "top_classifier_node": top_classifier_node,
    }


def _prepare_model_for_classification(model: Mapping[str, Any]) -> Dict[str, Any]:
    if model.get("_prepared_for_classification"):
        return dict(model)

    prepared = dict(model)
    prepared["architecture"] = _normalise_model_architecture(prepared)
    prepared["hierarchy"] = _normalise_hierarchy(prepared, prepared["architecture"])

    if (
        prepared["architecture"]["feature_family"] == "empirical"
        and prepared["architecture"]["topology"] == "hierarchical"
        and "node_models" in prepared
    ):
        return _prepare_hierarchical_model(prepared)

    prepared["is_hierarchical"] = False
    prepared["_prepared_for_classification"] = True
    prepared["_loci_by_id"] = {locus["locus_id"]: locus for locus in prepared.get("loci", [])}
    prepared["is_hierarchical"] = prepared["architecture"]["topology"] == "hierarchical"
    return prepared


def _is_legacy_empirical_hierarchical_model(model: Mapping[str, Any]) -> bool:
    if model.get("model_type") == HIERARCHICAL_MODEL_TYPE:
        return True
    hierarchy = model.get("hierarchy") or {}
    return bool(hierarchy.get("is_hierarchical")) and "node_models" in model


def _prepare_hierarchical_model(model: Mapping[str, Any]) -> Dict[str, Any]:
    prepared_node_models = {}
    flattened_loci = []

    for node_name, node_model in (model.get("node_models") or {}).items():
        locus_id_map = {}
        prepared_loci = []
        for locus in node_model.get("loci", []):
            prepared_locus = dict(locus)
            prepared_locus["raw_locus_id"] = locus["locus_id"]
            prepared_locus["node_name"] = node_name
            prepared_locus["locus_id"] = _qualify_hierarchical_locus_id(node_name, locus["locus_id"])
            locus_id_map[locus["locus_id"]] = prepared_locus["locus_id"]
            prepared_loci.append(prepared_locus)
            flattened_loci.append(prepared_locus)

        prepared_lineages = []
        for lineage in node_model.get("lineages", []):
            prepared_lineage = dict(lineage)
            for field in ("direct_locus_ids", "inherited_locus_ids"):
                prepared_lineage[field] = [
                    locus_id_map[locus_id]
                    for locus_id in lineage.get(field, [])
                    if locus_id in locus_id_map
                ]
            prepared_lineages.append(prepared_lineage)

        prepared_node_model = dict(node_model)
        prepared_node_model["loci"] = prepared_loci
        prepared_node_model["lineages"] = prepared_lineages
        prepared_node_model["_loci_by_id"] = {
            locus["locus_id"]: locus for locus in prepared_loci
        }
        prepared_node_model["architecture"] = {
            "feature_family": "empirical",
            "topology": "flat",
            "classification_workflow": "empirical_flat",
            "emission_model": model["architecture"].get("emission_model", "empirical_bayes"),
        }
        prepared_node_model["_prepared_for_classification"] = True
        prepared_node_model["is_hierarchical"] = False
        prepared_node_models[node_name] = prepared_node_model

    prepared = dict(model)
    prepared["node_models"] = prepared_node_models
    prepared["loci"] = flattened_loci
    prepared["_loci_by_id"] = {locus["locus_id"]: locus for locus in flattened_loci}
    prepared["_prepared_for_classification"] = True
    prepared["is_hierarchical"] = True
    return prepared


def _qualify_hierarchical_locus_id(node_name: str, locus_id: str) -> str:
    return f"{node_name}::{locus_id}"


def _iter_unique_bed_loci(model: Mapping[str, Any]) -> Sequence[Mapping[str, Any]]:
    unique_loci = []
    seen_keys = set()
    for locus in model.get("loci", []):
        key = (
            str(locus["chrom"]),
            int(locus["bed_start"]),
            int(locus["bed_end"]),
            locus.get("ref"),
            str(locus["alt"]) if locus.get("alt") is not None else None,
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)
        unique_loci.append(locus)
    return unique_loci


def _model_locus_by_id(model: Mapping[str, Any], locus_id: str) -> Mapping[str, Any]:
    loci_by_id = model.get("_loci_by_id")
    if loci_by_id is not None:
        try:
            return loci_by_id[locus_id]
        except KeyError as exc:
            raise KeyError(f"Locus {locus_id!r} not found in model.") from exc

    for locus in model["loci"]:
        if locus["locus_id"] == locus_id:
            return locus
    raise KeyError(f"Locus {locus_id!r} not found in model.")


def _classify_flat_evidence_against_model(
    parsed_evidence: ParsedVcfEvidence,
    model: Mapping[str, Any],
    mode: str,
    lineage_scope: str,
    kappa_weighting: bool,
    min_kappa_weight: float,
    full_weight_kappa: Optional[float],
) -> Dict[str, Any]:
    lineage_field = _lineage_scope_field(lineage_scope)
    model_locus_ids = set(model.get("_loci_by_id", {}))

    lineage_scores = []
    usable_evidence = [
        evidence
        for evidence in parsed_evidence.evidence_by_locus.values()
        if evidence.usable and evidence.locus_id in model_locus_ids
    ]

    if not usable_evidence:
        raise ValueError("No usable loci were available for classification.")

    for lineage in model["lineages"]:
        lineage_id = lineage["lineage_id"]
        active_locus_ids = set(lineage.get(lineage_field, []))
        log_score = math.log(float(lineage.get("prior", 0.0) or 1e-300))
        informative_loci = 0
        matched_loci = 0

        for evidence in usable_evidence:
            informative_loci += 1
            is_active = evidence.locus_id in active_locus_ids
            if is_active:
                matched_loci += 1

            log_score += _locus_log_likelihood(
                evidence=evidence,
                locus=_model_locus_by_id(model, evidence.locus_id),
                is_active_for_lineage=is_active,
                mode=mode,
                model=model,
                kappa_weighting=kappa_weighting,
                min_kappa_weight=min_kappa_weight,
                full_weight_kappa=full_weight_kappa,
            )

        lineage_scores.append(
            {
                "lineage_id": lineage_id,
                "log_score": log_score,
                "matched_loci": matched_loci,
                "informative_loci": informative_loci,
            }
        )

    posterior_scores = _normalise_log_scores(lineage_scores)
    posterior_scores.sort(key=lambda item: item["posterior"], reverse=True)

    return {
        "mode": mode,
        "lineage_scope": lineage_scope,
        "kappa_weighting": bool(kappa_weighting),
        "min_kappa_weight": float(min_kappa_weight),
        "full_weight_kappa": None if full_weight_kappa is None else float(full_weight_kappa),
        "sample_name": parsed_evidence.sample_name,
        "best_lineage": posterior_scores[0]["lineage_id"],
        "best_posterior": posterior_scores[0]["posterior"],
        "lineage_posteriors": posterior_scores,
        "evidence_summary": {
            "observed_loci": sum(
                1
                for evidence in parsed_evidence.evidence_by_locus.values()
                if evidence.observed and evidence.locus_id in model_locus_ids
            ),
            "usable_loci": len(usable_evidence),
            "total_model_loci": len(model["loci"]),
        },
    }


def _classify_hierarchical_evidence_against_model(
    parsed_evidence: ParsedVcfEvidence,
    model: Mapping[str, Any],
    mode: str,
    lineage_scope: str,
    kappa_weighting: bool,
    min_kappa_weight: float,
    full_weight_kappa: Optional[float],
) -> Dict[str, Any]:
    hierarchy = model["hierarchy"]
    node_models = model["node_models"]
    children_by_parent = {
        str(parent): list(children)
        for parent, children in (hierarchy.get("children_by_parent") or {}).items()
    }
    leaf_lineages = set(hierarchy.get("leaf_lineages", []))
    top_classifier_node = hierarchy.get("top_classifier_node")
    node_results: Dict[str, Dict[str, Any]] = {}

    def descend(node_name: str, incoming_posterior: float) -> Dict[str, float]:
        if node_name in leaf_lineages and node_name not in node_models:
            return {node_name: incoming_posterior}

        children = list(children_by_parent.get(node_name, []))
        if node_name not in node_models:
            if not children:
                return {node_name: incoming_posterior}
            if len(children) == 1:
                return descend(children[0], incoming_posterior)
            raise ValueError(
                f"Hierarchical node {node_name!r} has multiple children but no node model for classification."
            )

        local_result = _classify_flat_evidence_against_model(
            parsed_evidence=parsed_evidence,
            model=node_models[node_name],
            mode=mode,
            lineage_scope=lineage_scope,
            kappa_weighting=kappa_weighting,
            min_kappa_weight=min_kappa_weight,
            full_weight_kappa=full_weight_kappa,
        )
        node_results[node_name] = local_result
        child_posteriors = {
            item["lineage_id"]: item["posterior"]
            for item in local_result["lineage_posteriors"]
        }

        leaf_posteriors: Dict[str, float] = {}
        for child in children:
            child_posterior = incoming_posterior * child_posteriors.get(child, 0.0)
            for leaf, leaf_posterior in descend(child, child_posterior).items():
                leaf_posteriors[leaf] = leaf_posteriors.get(leaf, 0.0) + leaf_posterior
        return leaf_posteriors

    if top_classifier_node is None:
        root_lineages = list(hierarchy.get("root_lineages") or ())
        if len(root_lineages) == 1:
            top_classifier_node = root_lineages[0]
        else:
            raise ValueError("Hierarchical model is missing hierarchy.top_classifier_node.")

    leaf_posteriors = descend(str(top_classifier_node), 1.0)
    posterior_total = sum(leaf_posteriors.values())
    if posterior_total <= 0.0:
        raise ValueError("Hierarchical classification did not produce any leaf posterior mass.")

    lineage_posteriors = [
        {"lineage_id": lineage_id, "posterior": posterior / posterior_total}
        for lineage_id, posterior in leaf_posteriors.items()
    ]
    lineage_posteriors.sort(key=lambda item: item["posterior"], reverse=True)

    return {
        "mode": mode,
        "lineage_scope": lineage_scope,
        "kappa_weighting": bool(kappa_weighting),
        "min_kappa_weight": float(min_kappa_weight),
        "full_weight_kappa": None if full_weight_kappa is None else float(full_weight_kappa),
        "sample_name": parsed_evidence.sample_name,
        "best_lineage": lineage_posteriors[0]["lineage_id"],
        "best_posterior": lineage_posteriors[0]["posterior"],
        "lineage_posteriors": lineage_posteriors,
        "hierarchical": True,
        "hierarchy": {
            "top_classifier_node": top_classifier_node,
            "root_lineages": list(hierarchy.get("root_lineages") or ()),
            "leaf_lineages": sorted(leaf_lineages),
        },
        "node_classifications": node_results,
        "evidence_summary": {
            "observed_loci": len(parsed_evidence.observed_locus_ids),
            "usable_loci": len(parsed_evidence.usable_locus_ids),
            "total_model_loci": len(model["loci"]),
        },
    }


def _canonical_locus_state(evidence: VariantEvidence) -> str:
    if not evidence.usable:
        return "no_call"
    if (evidence.alt_depth or 0) > 0:
        return "alt_present"
    return "ref_callable"


def _classify_flat_canonical_evidence_against_model(
    parsed_evidence: ParsedVcfEvidence,
    model: Mapping[str, Any],
    lineage_scope: str,
) -> Dict[str, Any]:
    lineage_scores, evidence_summary = _score_canonical_lineages(
        parsed_evidence=parsed_evidence,
        model=model,
        candidate_lineages=model["lineages"],
        lineage_scope=lineage_scope,
    )

    return {
        "mode": "canonical",
        "lineage_scope": lineage_scope,
        "sample_name": parsed_evidence.sample_name,
        "best_lineage": lineage_scores[0]["lineage_id"],
        "best_posterior": lineage_scores[0]["posterior"],
        "lineage_posteriors": lineage_scores,
        "evidence_summary": evidence_summary,
    }


def _score_canonical_lineages(
    parsed_evidence: ParsedVcfEvidence,
    model: Mapping[str, Any],
    candidate_lineages: Sequence[Mapping[str, Any]],
    lineage_scope: str,
) -> Tuple[list[Dict[str, Any]], Dict[str, Any]]:
    lineage_field = _lineage_scope_field(lineage_scope)
    model_locus_ids = set(model.get("_loci_by_id", {}))
    present_locus_ids = {
        evidence.locus_id
        for evidence in parsed_evidence.evidence_by_locus.values()
        if evidence.locus_id in model_locus_ids and _canonical_locus_state(evidence) == "alt_present"
    }
    callable_locus_ids = {
        evidence.locus_id
        for evidence in parsed_evidence.evidence_by_locus.values()
        if evidence.locus_id in model_locus_ids and _canonical_locus_state(evidence) != "no_call"
    }

    lineage_scores = []
    for lineage in candidate_lineages:
        active_locus_ids = set(lineage.get(lineage_field, ())) & model_locus_ids
        matched_loci = len(active_locus_ids & present_locus_ids)
        callable_active_loci = len(active_locus_ids & callable_locus_ids)
        missing_active_loci = len(active_locus_ids) - callable_active_loci
        log_score = float(matched_loci) + math.log(float(lineage.get("prior", 0.0) or 1e-300))

        lineage_scores.append(
            {
                "lineage_id": lineage["lineage_id"],
                "log_score": log_score,
                "matched_loci": matched_loci,
                "active_loci": len(active_locus_ids),
                "callable_active_loci": callable_active_loci,
                "missing_active_loci": missing_active_loci,
            }
        )

    posterior_scores = _normalise_log_scores(lineage_scores)
    posterior_scores.sort(
        key=lambda item: (
            item["posterior"],
            item["matched_loci"],
            item["callable_active_loci"],
            -item["missing_active_loci"],
            str(item["lineage_id"]),
        ),
        reverse=True,
    )

    evidence_summary = {
        "observed_loci": sum(
            1
            for evidence in parsed_evidence.evidence_by_locus.values()
            if evidence.observed and evidence.locus_id in model_locus_ids
        ),
        "usable_loci": sum(
            1
            for evidence in parsed_evidence.evidence_by_locus.values()
            if evidence.usable and evidence.locus_id in model_locus_ids
        ),
        "total_model_loci": len(model.get("loci", ())),
    }
    return posterior_scores, evidence_summary


def _classify_hierarchical_canonical_evidence_against_model(
    parsed_evidence: ParsedVcfEvidence,
    model: Mapping[str, Any],
    lineage_scope: str,
) -> Dict[str, Any]:
    hierarchy = model["hierarchy"]
    children_by_parent = {
        str(parent): list(children)
        for parent, children in (hierarchy.get("children_by_parent") or {}).items()
    }
    lineages_by_id = {
        str(lineage["lineage_id"]): lineage
        for lineage in model.get("lineages", [])
    }
    top_classifier_node = hierarchy.get("top_classifier_node")
    if top_classifier_node is None:
        root_lineages = list(hierarchy.get("root_lineages") or ())
        if len(root_lineages) == 1:
            top_classifier_node = root_lineages[0]
        else:
            raise ValueError("Hierarchical canonical model is missing hierarchy.top_classifier_node.")

    node_results: Dict[str, Dict[str, Any]] = {}
    resolved_path = []
    current_node = str(top_classifier_node)

    while True:
        children = [str(child) for child in children_by_parent.get(current_node, ())]
        if not children:
            best_lineage = current_node
            break

        candidate_lineages = [
            lineages_by_id[child]
            for child in children
            if child in lineages_by_id
        ]
        if not candidate_lineages:
            if len(children) == 1:
                current_node = children[0]
                resolved_path.append(current_node)
                continue
            raise ValueError(
                f"Hierarchical canonical node {current_node!r} has multiple children but no lineage records."
            )

        lineage_scores, evidence_summary = _score_canonical_lineages(
            parsed_evidence=parsed_evidence,
            model=model,
            candidate_lineages=candidate_lineages,
            lineage_scope="direct",
        )
        local_result = {
            "mode": "canonical",
            "lineage_scope": "direct",
            "sample_name": parsed_evidence.sample_name,
            "best_lineage": lineage_scores[0]["lineage_id"],
            "best_posterior": lineage_scores[0]["posterior"],
            "lineage_posteriors": lineage_scores,
            "evidence_summary": evidence_summary,
            "node_name": current_node,
            "children": children,
        }
        node_results[current_node] = local_result

        current_is_lineage = current_node in lineages_by_id
        if lineage_scores[0]["matched_loci"] <= 0 and current_is_lineage:
            best_lineage = current_node
            break

        current_node = str(local_result["best_lineage"])
        resolved_path.append(current_node)

    model_locus_ids = set(model.get("_loci_by_id", {}))
    return {
        "mode": "canonical",
        "lineage_scope": lineage_scope,
        "sample_name": parsed_evidence.sample_name,
        "best_lineage": best_lineage,
        "best_posterior": 1.0,
        "lineage_posteriors": [{"lineage_id": best_lineage, "posterior": 1.0}],
        "hierarchical": True,
        "hierarchy": {
            "top_classifier_node": top_classifier_node,
            "root_lineages": list(hierarchy.get("root_lineages") or ()),
            "leaf_lineages": sorted(str(lineage["lineage_id"]) for lineage in model.get("lineages", [])),
        },
        "resolved_path": resolved_path,
        "node_classifications": node_results,
        "evidence_summary": {
            "observed_loci": sum(
                1
                for evidence in parsed_evidence.evidence_by_locus.values()
                if evidence.observed and evidence.locus_id in model_locus_ids
            ),
            "usable_loci": sum(
                1
                for evidence in parsed_evidence.evidence_by_locus.values()
                if evidence.usable and evidence.locus_id in model_locus_ids
            ),
            "total_model_loci": len(model.get("loci", ())),
        },
    }


def _locus_log_likelihood(
    evidence: VariantEvidence,
    locus: Mapping[str, Any],
    is_active_for_lineage: bool,
    mode: str,
    model: Mapping[str, Any],
    kappa_weighting: bool,
    min_kappa_weight: float,
    full_weight_kappa: Optional[float],
) -> float:
    if evidence.depth is None or evidence.alt_depth is None:
        raise ValueError(f"Locus {evidence.locus_id} is missing depth information required for classification.")

    if mode == "empirical_bayes":
        ## use the sparse target/background frequencies stored directly on the locus
        frequencies = locus["empirical_bayes"]
        p = (
            float(frequencies["target_frequency"])
            if is_active_for_lineage
            else float(frequencies["background_frequency"])
        )
        log_likelihood = _log_binomial_pmf(evidence.alt_depth, evidence.depth, p)
        return log_likelihood * _kappa_weight_for_locus(
            locus=locus,
            model=model,
            is_active_for_lineage=is_active_for_lineage,
            kappa_weighting=kappa_weighting,
            min_kappa_weight=min_kappa_weight,
            full_weight_kappa=full_weight_kappa,
        )

    if mode == "full_bayes":
        ## use prefit BB parameters stored in the model
        alpha, beta = _resolve_full_bayes_beta_params(locus, model, is_active_for_lineage)
        log_likelihood = _log_beta_binomial_pmf(evidence.alt_depth, evidence.depth, alpha, beta)
        return log_likelihood * _kappa_weight_for_locus(
            locus=locus,
            model=model,
            is_active_for_lineage=is_active_for_lineage,
            kappa_weighting=kappa_weighting,
            min_kappa_weight=min_kappa_weight,
            full_weight_kappa=full_weight_kappa,
        )

    raise ValueError("mode must be either 'empirical_bayes' or 'full_bayes'.")


def _resolve_full_bayes_beta_params(
    locus: Mapping[str, Any],
    model: Mapping[str, Any],
    is_active_for_lineage: bool,
) -> Tuple[float, float]:
    full_bayes = model.get("full_bayes") or {}
    if full_bayes.get("status") != "fit":
        raise ValueError("full_bayes was requested, but model['full_bayes']['status'] is not 'fit'.")

    locus_params = locus.get("full_bayes", {})
    param_prefix = "target" if is_active_for_lineage else "background"

    alpha = locus_params.get(f"{param_prefix}_alpha")
    beta = locus_params.get(f"{param_prefix}_beta")
    if alpha is not None and beta is not None:
        return float(alpha), float(beta)

    mu = locus_params.get(f"{param_prefix}_mu")
    kappa = locus_params.get(f"{param_prefix}_kappa")
    if mu is not None and kappa is not None:
        return _mu_kappa_to_alpha_beta(float(mu), float(kappa))

    defaults = full_bayes.get("defaults", {})
    alpha = defaults.get(f"{param_prefix}_alpha")
    beta = defaults.get(f"{param_prefix}_beta")
    if alpha is not None and beta is not None:
        return float(alpha), float(beta)

    mu = defaults.get(f"{param_prefix}_mu")
    kappa = defaults.get(f"{param_prefix}_kappa")
    if mu is not None and kappa is not None:
        return _mu_kappa_to_alpha_beta(float(mu), float(kappa))

    raise ValueError(
        "full_bayes parameters were not found for this locus. "
        "Expected alpha/beta or mu/kappa values on the locus or in model['full_bayes']['defaults']."
    )


def _mu_kappa_to_alpha_beta(mu: float, kappa: float) -> Tuple[float, float]:
    if not 0.0 < mu < 1.0:
        raise ValueError(f"mu must be strictly between 0 and 1, got {mu}.")
    if kappa <= 0.0:
        raise ValueError(f"kappa must be positive, got {kappa}.")
    return mu * kappa, (1.0 - mu) * kappa


def _resolve_min_groups_for_full_bayes(model: Mapping[str, Any]) -> Optional[int]:
    full_bayes = model.get("full_bayes") or {}
    min_groups = full_bayes.get("min_groups_for_full_bayes")
    if min_groups is None:
        min_groups = ((model.get("provenance") or {}).get("loso") or {}).get("min_groups_for_full_bayes")
    if min_groups is None:
        return None
    return int(min_groups)


def _fallback_kappa_scale(
    fit_status: Optional[str],
    group_count: Optional[int],
    min_groups_for_full_bayes: Optional[int],
) -> float:
    if fit_status in (None, "fit"):
        return 1.0

    if min_groups_for_full_bayes is not None and min_groups_for_full_bayes > 0 and group_count is not None:
        return min(max(float(group_count) / float(min_groups_for_full_bayes), 0.0), 1.0)

    return 0.5


def _resolve_full_bayes_kappa(
    locus: Mapping[str, Any],
    model: Mapping[str, Any],
    is_active_for_lineage: bool,
) -> Optional[float]:
    full_bayes = model.get("full_bayes") or {}
    if full_bayes.get("status") != "fit":
        return None

    locus_params = locus.get("full_bayes", {})
    param_prefix = "target" if is_active_for_lineage else "background"
    kappa = locus_params.get(f"{param_prefix}_kappa")
    fit_status = locus_params.get(f"{param_prefix}_fit_status")
    group_count = locus_params.get(f"{param_prefix}_group_count")
    if kappa is not None:
        scale = _fallback_kappa_scale(
            fit_status=fit_status,
            group_count=None if group_count is None else int(group_count),
            min_groups_for_full_bayes=_resolve_min_groups_for_full_bayes(model),
        )
        return float(kappa) * float(scale)

    defaults = full_bayes.get("defaults", {})
    kappa = defaults.get(f"{param_prefix}_kappa")
    if kappa is not None:
        return float(kappa)

    return None


def _kappa_weight_for_locus(
    locus: Mapping[str, Any],
    model: Mapping[str, Any],
    is_active_for_lineage: bool,
    kappa_weighting: bool,
    min_kappa_weight: float,
    full_weight_kappa: Optional[float],
) -> float:
    if not kappa_weighting:
        return 1.0

    if not 0.0 < float(min_kappa_weight) <= 1.0:
        raise ValueError(f"min_kappa_weight must be in (0, 1], got {min_kappa_weight}.")

    kappa = _resolve_full_bayes_kappa(
        locus=locus,
        model=model,
        is_active_for_lineage=is_active_for_lineage,
    )
    if kappa is None:
        return 1.0

    if full_weight_kappa is None:
        raw_weight = kappa / (kappa + 1.0)
    else:
        raw_weight = min(kappa / float(full_weight_kappa), 1.0)

    return float(min(max(raw_weight, float(min_kappa_weight)), 1.0))


def _log_binomial_pmf(k: int, n: int, p: float) -> float:
    p = min(max(p, 1e-12), 1.0 - 1e-12)
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log(1.0 - p)
    )


def _log_beta_binomial_pmf(k: int, n: int, alpha: float, beta: float) -> float:
    if alpha <= 0.0 or beta <= 0.0:
        raise ValueError(f"alpha and beta must be positive, got alpha={alpha}, beta={beta}.")
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + math.lgamma(k + alpha)
        + math.lgamma(n - k + beta)
        - math.lgamma(n + alpha + beta)
        + math.lgamma(alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
    )


def _normalise_log_scores(lineage_scores: Sequence[Dict[str, Any]]) -> list[Dict[str, Any]]:
    max_log_score = max(item["log_score"] for item in lineage_scores)
    normaliser = sum(math.exp(item["log_score"] - max_log_score) for item in lineage_scores)

    normalised = []
    for item in lineage_scores:
        posterior = math.exp(item["log_score"] - max_log_score) / normaliser
        enriched = dict(item)
        enriched["posterior"] = posterior
        normalised.append(enriched)

    return normalised

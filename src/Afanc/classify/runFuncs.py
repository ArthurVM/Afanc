from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from Afanc.classify.allele_json import load_allele_json
from Afanc.classify.reference_validation import validate_allele_reference, validate_vcf_reference
from Afanc.utilities.generalUtils import vprint
from Afanc.screen.profiles import (
    load_profiles_manifest,
    profile_paths,
    resolve_profile_for_species,
    validate_profile_reference,
)
from Afanc.screen.variant_calling.profile import (
    missing_positions_from_bed,
    sample_alleles_from_vcf,
)
from Afanc.classifier.bayesian_profile import classify_snp_payload


def runClassify(args: Any) -> int:
    """Classify lineage from externally supplied VCF or allele JSON evidence."""
    subprocessID = "CLASSIFY"
    vprint(
        subprocessID,
        f"Resolving profile model for {args.species}...",
        "prYellow",
        getattr(args, "stdout", None),
    )
    profiles_dir, manifest = load_profiles_manifest(
        getattr(args, "database", "."),
        profiles_dir=getattr(args, "profiles_dir", None),
    )
    if profiles_dir is None:
        raise ValueError("No profiles.json found. Provide --profiles-dir or a database with profiles.")

    profile = resolve_profile_for_species(args.species, manifest)
    if profile is None:
        supported = supported_reference_species(manifest)
        supported_text = "\n".join(f"  - {name}" for name in supported) if supported else "  (none)"
        vprint(
            "ERROR",
            f"No enabled profile model found for {args.species!r}.\nSupported reference species/profile names:\n{supported_text}",
            "prRed",
            getattr(args, "stderr", None),
        )
        raise ValueError(
            f"No enabled profile model found for species/profile name {args.species!r}. "
            f"Supported reference species/profile names: {', '.join(supported) if supported else '(none)'}"
        )

    vprint(
        subprocessID,
        "Validating profile model against its reference FASTA...",
        "prYellow",
        getattr(args, "stdout", None),
    )
    profile_validation = validate_profile_reference(profile, profiles_dir)
    if not profile_validation["valid"]:
        raise ValueError("Profile reference validation failed: " + "; ".join(profile_validation["errors"]))

    paths = profile_paths(profile, profiles_dir)
    reference_fasta = Path(paths["reference"])
    model_path = Path(paths["model"])

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    sample_id = args.sample_id or output_prefix.name
    snp_json_path = output_prefix.with_suffix(".snps.json")
    classification_json = output_prefix.with_suffix(".lineage_classification.json")
    summary_json = output_prefix.with_suffix(".classify.json")

    if args.vcf:
        input_type = "vcf"
        input_path = Path(args.vcf)
        vprint(
            subprocessID,
            f"Validating VCF reference compatibility for {input_path}...",
            "prYellow",
            getattr(args, "stdout", None),
        )
        reference_validation = validate_vcf_reference(input_path, reference_fasta)
        if not reference_validation["valid"]:
            raise ValueError("VCF reference validation failed: " + "; ".join(reference_validation["errors"]))
        vprint(
            subprocessID,
            "Converting VCF evidence to Afanc SNP JSON payload...",
            "prYellow",
            getattr(args, "stdout", None),
        )
        sample_payload = _payload_from_vcf(args, input_path)
    else:
        input_type = "allele_json"
        input_path = Path(args.allele_json)
        if not args.allele_id_format:
            raise ValueError("--allele-id-format is required with --allele-json.")
        vprint(
            subprocessID,
            f"Loading allele JSON evidence from {input_path}...",
            "prYellow",
            getattr(args, "stdout", None),
        )
        sample_payload = load_allele_json(input_path)
        vprint(
            subprocessID,
            "Validating allele IDs against the profile reference FASTA...",
            "prYellow",
            getattr(args, "stdout", None),
        )
        reference_validation = validate_allele_reference(
            sample_payload["alleles"],
            args.allele_id_format,
            reference_fasta,
        )
        if not reference_validation["valid"]:
            raise ValueError("Allele JSON reference validation failed: " + "; ".join(reference_validation["errors"]))

    with snp_json_path.open("w") as handle:
        json.dump(sample_payload, handle, indent=2)

    vprint(
        subprocessID,
        f"Running Bayesian lineage classification for {sample_id}...",
        "prYellow",
        getattr(args, "stdout", None),
    )
    classification_result = classify_snp_payload(
        sample_payload,
        model_path,
        classification_json,
        sample_id,
        min_support=args.lineage_min_support,
        min_support_fraction=args.lineage_min_support_fraction,
        min_callable_fraction=args.lineage_min_callable_fraction,
        ambiguity_score_margin=args.lineage_ambiguity_margin,
        allow_incomplete_descent=not args.lineage_disable_incomplete_descent,
        infer_reference_allele_markers=not args.lineage_disable_reference_marker_inference,
        tie_delta=args.lineage_tie_delta,
    )

    summary = {
        "sample_id": sample_id,
        "species": args.species,
        "profile": {
            "taxon_id": profile.get("taxon_id"),
            "name": profile.get("name"),
            "model": str(model_path),
            "reference": str(reference_fasta),
            "profile_match": profile.get("profile_match"),
        },
        "input": {
            "type": input_type,
            "path": str(input_path),
            "reference_validation": reference_validation,
        },
        "snp_json": str(snp_json_path),
        "snp_count": len(sample_payload["alleles"]),
        "missing_count": len(sample_payload["missing"]),
        "classification_json": classification_result["classification_json"],
        "classification": {
            "status": classification_result["status"],
            "best_lineage": classification_result["best_lineage"],
            "best_posterior": classification_result["best_posterior"],
            "call_status": classification_result["call_status"],
            "classifier": classification_result["classifier"],
            "model_family": classification_result["model_family"],
        },
    }
    with summary_json.open("w") as handle:
        json.dump(summary, handle, indent=2)

    vprint(
        "FINISHED",
        f"Classification summary written to {summary_json}",
        "prGreen",
        getattr(args, "stdout", None),
    )
    return 0


def _payload_from_vcf(args: Any, vcf_path: Path) -> Dict[str, Any]:
    allele_map = sample_alleles_from_vcf(
        vcf_path=vcf_path,
        min_qual=args.snp_min_qual,
        allow_filtered=args.snp_allow_filtered,
        accept_missing_qual=args.snp_accept_missing_qual,
        min_dp=args.snp_min_dp,
    )
    alleles = set()
    for sample_alleles in allele_map.values():
        alleles.update(sample_alleles)

    missing = []
    if args.depth_bed:
        missing = missing_positions_from_bed(Path(args.depth_bed), args.snp_min_missing_depth)

    return {
        "alleles": sorted(alleles),
        "missing": missing,
    }


def supported_reference_species(manifest: Dict[str, Any]) -> list[str]:
    profiles = manifest.get("profiles") or {}
    names = []
    for profile in profiles.values():
        if not profile or not profile.get("enabled", True):
            continue
        if profile.get("name"):
            names.append(str(profile["name"]))
        names.extend(str(alias) for alias in profile.get("aliases", []) if alias)
    return sorted(set(names), key=str.lower)

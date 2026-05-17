from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping

from .bayesian_classifier import classify_samples, load_model_json, load_sample_alleles_json, write_classifications_json


def run_lineage_classification(
    args: Any,
    snp_profile: Mapping[str, Mapping[str, Any]],
    mapped_bams: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    lineage_profiles = {}

    for accession, snp_record in snp_profile.items():
        mapping_record = mapped_bams.get(accession, {})
        profile = mapping_record.get("lineage_profile")
        if not profile:
            lineage_profiles[accession] = {
                "accession": accession,
                "status": "not_run",
                "reason": "no_profile_model",
                "snp_json": str(snp_record.get("snp_json")),
                "snp_count": snp_record.get("snp_count"),
                "missing_count": snp_record.get("missing_count"),
            }
            continue

        model_path = Path(profile["model"])
        snp_json_path = Path(snp_record["snp_json"])
        output_json = Path(args.reportsDir) / f"{accession}.lineage_classification.json"

        model = load_model_json(model_path)
        sample_payload = load_sample_alleles_json(snp_json_path)
        classifications = classify_samples(
            sample_payload,
            model,
            sample_id=args.output_prefix,
            min_support=args.lineage_min_support,
            min_support_fraction=args.lineage_min_support_fraction,
            min_callable_fraction=args.lineage_min_callable_fraction,
            ambiguity_score_margin=args.lineage_ambiguity_margin,
            allow_incomplete_descent=not args.lineage_disable_incomplete_descent,
            infer_reference_allele_markers=not args.lineage_disable_reference_marker_inference,
            tie_delta=args.lineage_tie_delta,
        )
        write_classifications_json(classifications, output_json)

        best = classifications[0] if classifications else {}
        lineage_profiles[accession] = {
            "accession": accession,
            "status": "classified" if classifications else "not_determined",
            "taxon_id": profile.get("taxon_id"),
            "name": profile.get("name"),
            "model": str(model_path),
            "reference": profile.get("reference"),
            "profile_match": profile.get("profile_match"),
            "snp_json": str(snp_json_path),
            "classification_json": str(output_json),
            "best_lineage": best.get("best_lineage"),
            "best_posterior": best.get("best_posterior"),
            "call_status": best.get("call_status"),
            "classifier": best.get("classifier"),
            "model_family": best.get("model_family"),
        }

    return lineage_profiles

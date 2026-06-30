import json
from os import listdir, path

from Afanc.utilities.generalUtils import vprint
from Afanc.utilities.getVersions import getVersionsScreen


def makeFinalReport(args, reports, snp_profile=None, lineage_profile=None):
    """Make final report from the kraken2 and mapping reports."""
    subprocessID = "REPORT"
    vprint(
        subprocessID,
        "Generating final report...",
        "prYellow"
    )

    jsondict = {
        "schema_version": "2.0",
        "results": {
            "detections": [],
        },
        "run": {
            "inputs": {
                "database": str(args.database),
                "fastqs": [str(fastq) for fastq in args.fastq],
            },
            "settings": _report_settings(args),
            "software_versions": {
                "screen": getVersionsScreen(),
            },
        },
    }
    clustering_results = []

    mapping_reports_box = [g for g in listdir(args.reportsDir) if g.endswith("mapstats.json")]
    final_report = path.abspath(f"./{args.output_prefix}.json")

    with open(f"{args.reportsDir}/{args.output_prefix}.k2.json", 'r') as k2fin:
        k2data = json.load(k2fin)

        for event in k2data["Detection_events"]:

            event["warnings"] = event.get("warnings", [])
            if "lineage_profile_validation" in event and not event["lineage_profile_validation"].get("valid", False):
                event["warnings"].append({
                    "code": "lineage_profile_validation",
                    "details": event["lineage_profile_validation"].get("errors", []),
                })

            mapping_target = (
                event["closest_variant"]
                if isinstance(event.get("closest_variant"), dict)
                else event
            )
            assembly = mapping_target.get("assembly") or event.get("assembly")
            if event.get("lineage_profile") and event.get("assembly"):
                assembly = event["assembly"]

            if assembly is not None:
                assembly_prefix = _assembly_to_accession(assembly)
                report_path = reports.get(assembly_prefix) if reports else None
                if report_path is None:
                    report_matches = [
                        report
                        for report in mapping_reports_box
                        if report.endswith(assembly_prefix + ".mapstats.json")
                    ]
                    report_path = f"{args.reportsDir}/{report_matches[0]}" if report_matches else None

                if report_path is None:
                    event["warnings"].append({
                        "code": "mapping_statistics_unavailable",
                        "message": (
                            f"No mapping stats report could be found for assembly {assembly} "
                            f"using key {assembly_prefix}."
                        ),
                    })
                else:
                    with open(report_path, 'r') as mapping_report_handle:
                        mapping_data = json.load(mapping_report_handle)
                    mapping_warnings = mapping_data.get("warnings", {})
                    if "no_unique_map" not in mapping_warnings:
                        mapping_target.update({
                            "mean_DOC": mapping_data["map_data"]["mean_DOC"],
                            "median_DOC": mapping_data["map_data"]["median_DOC"],
                            "reference_cov": mapping_data["map_data"]["proportion_cov"],
                            "gini": mapping_data["map_data"]["gini"],
                        })
                    _append_mapping_warnings(event["warnings"], mapping_warnings)

            elif not getattr(args, "no_map", False):
                event["warnings"].append({
                    "code": "no_reference_assembly",
                    "message": f"No assembly could be found for hit on {event['name']}.",
                })

            clustering_results.append(event)
            jsondict["results"]["detections"].append(
                _build_detection_summary(
                    event,
                    no_map=getattr(args, "no_map", False),
                    snp_profile=snp_profile,
                    lineage_profile=lineage_profile,
                )
            )

    jsondict["run"]["artifacts"] = _write_report_subreports(
        args,
        clustering_results=clustering_results,
        snp_profile=snp_profile,
        lineage_profile=lineage_profile,
    )

    with open(final_report, "w") as fout:
        json.dump(jsondict, fout, indent=2, default=str)

    return final_report


def _report_settings(args):
    """Group report settings by pipeline stage."""
    return {
        "output_prefix": args.output_prefix,
        "screening": {
            "minimum_reads": args.num_threshold,
            "minimum_percentage": args.pct_threshold,
            "redistribution_bounds": {
                "lower": args.lower_bound,
                "upper": args.upper_bound,
            },
        },
        "mapping": {
            "enabled": not getattr(args, "no_map", False),
        },
        "variant_calling": {
            "enabled": not getattr(args, "no_map", False),
            "caller": args.variant_caller,
            "filters": {
                "minimum_quality": args.snp_min_qual,
                "minimum_depth": args.snp_min_dp,
                "minimum_missing_depth": args.snp_min_missing_depth,
                "allow_filtered": args.snp_allow_filtered,
                "accept_missing_quality": args.snp_accept_missing_qual,
            },
        },
        "lineage_classification": {
            "enabled": not args.no_lineage_classify,
            "profiles_directory": args.profiles_dir,
            "allow_compound_sample": args.lineage_profile_compound,
            "skip_mixed_taxa_by_default": True,
            "minimum_support": args.lineage_min_support,
            "minimum_support_fraction": args.lineage_min_support_fraction,
            "minimum_callable_fraction": args.lineage_min_callable_fraction,
            "ambiguity_margin": args.lineage_ambiguity_margin,
            "tie_delta": args.lineage_tie_delta,
            "reference_marker_inference": not args.lineage_disable_reference_marker_inference,
            "incomplete_descent": not args.lineage_disable_incomplete_descent,
        },
    }


def _append_mapping_warnings(event_warnings, mapping_warnings):
    """Copy map-stat warnings into the final detection."""
    for warning_code, warning_payload in mapping_warnings.items():
        if isinstance(warning_payload, dict):
            warning = dict(warning_payload)
            warning.setdefault("code", warning_code)
        else:
            warning = {"code": warning_code, "message": str(warning_payload)}
        event_warnings.append(warning)


def _assembly_to_accession(assembly):
    assembly_prefix = path.basename(path.splitext(assembly)[0])
    if assembly_prefix.endswith("_genomic"):
        assembly_prefix = assembly_prefix[:-len("_genomic")]
    return assembly_prefix


def _build_detection_summary(event, no_map=False, snp_profile=None, lineage_profile=None):
    """Build one compact, human-readable detection result."""
    target = event.get("closest_variant") if isinstance(event.get("closest_variant"), dict) else event
    taxon = {
        "taxon_id": str(event["taxon_id"]),
        "name": event["name"],
    }
    if target is not event:
        taxon["closest_variant"] = {
            "taxon_id": str(target["taxon_id"]),
            "name": target["name"],
        }

    snp_record = _record_for_event(snp_profile, event)
    lineage_record = _record_for_event(lineage_profile, event)

    return {
        "taxon": taxon,
        "read_support": {
            "count": event.get("reads"),
            "percentage": event.get("percentage"),
        },
        "mapping": _mapping_summary(event, target, no_map=no_map),
        "variant_calling": _variant_calling_summary(snp_record, no_map=no_map),
        "lineage": _lineage_summary(event, lineage_record, no_map=no_map),
        "warnings": _normalise_warnings(event.get("warnings", [])),
    }


def _mapping_summary(event, target, no_map=False):
    assembly = target.get("assembly") or event.get("assembly")
    accession = target.get("accession") or event.get("accession")
    if accession in (None, "", "NA") and assembly:
        accession = _assembly_to_accession(assembly)

    reference = {}
    if accession not in (None, "", "NA"):
        reference["accession"] = str(accession)
    if assembly:
        reference["assembly"] = path.basename(str(assembly))

    if no_map:
        summary = {
            "status": "not_run",
            "reason": "disabled_by_no_map",
        }
    elif "mean_DOC" in target:
        summary = {
            "status": "completed",
            "mean_depth": target.get("mean_DOC"),
            "median_depth": target.get("median_DOC"),
            "genome_coverage_fraction": target.get("reference_cov"),
            "depth_gini": target.get("gini"),
        }
    elif _has_warning_code(event, "no_unique_map"):
        summary = {
            "status": "not_determined",
            "reason": "no_unique_mappings",
        }
    else:
        summary = {
            "status": "not_determined",
            "reason": "mapping_statistics_unavailable",
        }

    if reference:
        summary["reference"] = reference
    return summary


def _variant_calling_summary(snp_record, no_map=False):
    if snp_record is None:
        return {
            "status": "not_run",
            "reason": "mapping_disabled" if no_map else "no_variant_calls",
        }

    return {
        "status": "completed",
        "caller": snp_record.get("caller"),
        "snp_count": snp_record.get("snp_count"),
        "missing_position_count": snp_record.get("missing_count"),
        "allele_id_format": snp_record.get("allele_id_format"),
    }


def _lineage_summary(event, lineage_record, no_map=False):
    if lineage_record is None:
        return {
            "status": "not_run",
            "reason": "mapping_disabled" if no_map else "no_lineage_result",
        }

    summary = {
        "status": lineage_record.get("status", "not_determined"),
    }
    if lineage_record.get("reason"):
        summary["reason"] = lineage_record["reason"]

    profile_record = lineage_record
    if not profile_record.get("taxon_id") and isinstance(event.get("lineage_profile"), dict):
        profile_record = event["lineage_profile"]

    profile = {}
    if profile_record.get("taxon_id") is not None:
        profile["taxon_id"] = str(profile_record["taxon_id"])
    if profile_record.get("name"):
        profile["name"] = profile_record["name"]
    if profile_record.get("profile_match"):
        profile["match"] = profile_record["profile_match"]
    if profile:
        summary["profile"] = profile

    if lineage_record.get("best_lineage") is not None:
        summary["call"] = {
            "name": lineage_record["best_lineage"],
            "posterior": lineage_record.get("best_posterior"),
            "status": lineage_record.get("call_status"),
        }

    method = {}
    if lineage_record.get("classifier"):
        method["classifier"] = lineage_record["classifier"]
    if lineage_record.get("model_family"):
        method["model_family"] = lineage_record["model_family"]
    if method:
        summary["method"] = method

    return summary


def _record_for_event(records, event):
    if not records:
        return None

    candidates = []
    target = event.get("closest_variant") if isinstance(event.get("closest_variant"), dict) else event
    for node in (target, event):
        if node.get("assembly"):
            candidates.append(_assembly_to_accession(node["assembly"]))
        if node.get("accession") not in (None, "", "NA"):
            candidates.append(str(node["accession"]))

    for candidate in candidates:
        if candidate in records:
            return records[candidate]
    return None


def _normalise_warnings(warnings):
    return [
        warning
        if isinstance(warning, dict) and "code" in warning
        else {"code": "general_warning", "message": str(warning)}
        for warning in warnings or []
    ]


def _has_warning_code(event, warning_code):
    return any(
        isinstance(warning, dict)
        and warning.get("code") == warning_code
        for warning in event.get("warnings", [])
    )


def _write_report_subreports(args, clustering_results, snp_profile=None, lineage_profile=None):
    artifacts = {}
    detail_reports = {
        "detection_details": ("detections.details", clustering_results),
        "variant_calling_details": ("variant_calling.details", snp_profile or {}),
        "lineage_classification_details": ("lineage_classification.details", lineage_profile or {}),
    }

    for artifact_name, (file_label, payload) in detail_reports.items():
        report_path = path.abspath(
            path.join(args.reportsDir, f"{args.output_prefix}.{file_label}.json")
        )
        with open(report_path, "w") as fout:
            json.dump(payload, fout, indent=4, default=str)
        artifacts[artifact_name] = {
            "path": report_path,
            "record_count": len(payload),
        }

    if hasattr(args, "kronaWDir"):
        artifacts["krona_report"] = {
            "path": path.abspath(
                path.join(args.kronaWDir, f"{args.output_prefix}.krona.html")
            ),
        }
    return artifacts

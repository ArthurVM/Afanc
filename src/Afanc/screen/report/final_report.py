import json
from collections import defaultdict
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

    jsondict = defaultdict(dict)
    datadict = jsondict["results"] = {}

    ## index detection events by taxid
    datadict["Detection_events"] = {}
    clustering_results = []

    jsondict["arguments"] = {
        "database" : str(args.database),
        "fastqs" : [str(fastq) for fastq in args.fastq],
        "num_threshold" : args.num_threshold,
        "pct_threshold" : args.pct_threshold,
        "upper_bound_threshold" : args.upper_bound,
        "lower_bound_thresold" : args.lower_bound,
        "variant_caller" : args.variant_caller,
        "lineage_profile_compound" : args.lineage_profile_compound,
        "snp_min_qual" : args.snp_min_qual,
        "snp_min_dp" : args.snp_min_dp,
        "snp_min_missing_depth" : args.snp_min_missing_depth,
        "snp_allow_filtered" : args.snp_allow_filtered,
        "snp_accept_missing_qual" : args.snp_accept_missing_qual,
        "lineage_classification" : "disabled" if args.no_lineage_classify else "profile_models_if_available",
        "lineage_skips_mixed_taxa_by_default" : True,
        "profiles_dir" : args.profiles_dir,
        "lineage_min_support" : args.lineage_min_support,
        "lineage_min_support_fraction" : args.lineage_min_support_fraction,
        "lineage_min_callable_fraction" : args.lineage_min_callable_fraction,
        "lineage_ambiguity_margin" : args.lineage_ambiguity_margin,
        "lineage_tie_delta" : args.lineage_tie_delta,
        "lineage_reference_marker_inference" : not args.lineage_disable_reference_marker_inference,
        "lineage_incomplete_descent" : not args.lineage_disable_incomplete_descent,
        "output_prefix" : args.output_prefix,
    }

    jsondict["versions"] = {
        "screen" : getVersionsScreen()
    }

    mapping_reports_box = [g for g in listdir(args.reportsDir) if g.endswith("mapstats.json")]
    final_report = path.abspath(f"./{args.output_prefix}.json")

    with open(f"{args.reportsDir}/{args.output_prefix}.k2.json", 'r') as k2fin:
        k2data = json.load(k2fin)

        for event in k2data["Detection_events"]:

            event["warnings"] = event.get("warnings", [])
            if "lineage_profile_validation" in event and not event["lineage_profile_validation"].get("valid", False):
                event["warnings"].append({"lineage_profile_validation": event["lineage_profile_validation"].get("errors", [])})

            if "closest_variant" in event:
                closest_variant_flag = True

                ## profile overrides store the assembly at event level
                if "lineage_profile" in event and "assembly" in event:
                    assembly = event["assembly"]
                elif "assembly" in event["closest_variant"]:
                    assembly = event["closest_variant"]["assembly"]
                elif "assembly" in event:
                    assembly = event["assembly"]
                else:
                    assembly = None

            elif "assembly" in event:

                closest_variant_flag = False
                assembly = event["assembly"]

            else:
                closest_variant_flag = False
                assembly = None

            if assembly is not None:
                assembly_prefix = _assembly_to_accession(assembly)

                ## fall back to legacy filename matching
                report_path = reports.get(assembly_prefix) if reports else None
                if report_path is None:
                    report_matches = [
                        report
                        for report in mapping_reports_box
                        if report.endswith(assembly_prefix + ".mapstats.json")
                    ]
                    report_path = f"{args.reportsDir}/{report_matches[0]}" if report_matches else None

                if report_path is None:
                    event["warnings"].append(
                        f"No mapping stats report could be found for assembly {assembly} using key {assembly_prefix}."
                    )
                    clustering_results.append(event)
                    _index_detection_event_by_taxid(
                        datadict["Detection_events"],
                        event,
                        snp_profile=snp_profile,
                        lineage_profile=lineage_profile,
                    )
                    continue

                with open(report_path, 'r') as mapping_report_handle:
                    mapping_data = json.load(mapping_report_handle)

                    if "no_unique_map" not in mapping_data["warnings"]:

                        if closest_variant_flag:
                            event["closest_variant"]["mean_DOC"] = mapping_data["map_data"]["mean_DOC"]
                            event["closest_variant"]["median_DOC"] = mapping_data["map_data"]["median_DOC"]
                            event["closest_variant"]["reference_cov"] = mapping_data["map_data"]["proportion_cov"]
                            event["closest_variant"]["gini"] = mapping_data["map_data"]["gini"]

                        else:
                            event["mean_DOC"] = mapping_data["map_data"]["mean_DOC"]
                            event["median_DOC"] = mapping_data["map_data"]["median_DOC"]
                            event["reference_cov"] = mapping_data["map_data"]["proportion_cov"]
                            event["gini"] = mapping_data["map_data"]["gini"]

                    else:
                        event["warnings"].append(mapping_data["warnings"])

            else:
                event['warnings'].append(f"No assembly could be found for hit on {event['name']}.")

            clustering_results.append(event)
            _index_detection_event_by_taxid(
                datadict["Detection_events"],
                event,
                snp_profile=snp_profile,
                lineage_profile=lineage_profile,
            )

    jsondict["sub_reports"] = _write_report_subreports(
        args,
        clustering_results=clustering_results,
        snp_profile=snp_profile,
        lineage_profile=lineage_profile,
    )

    with open(final_report, "w") as fout:
        json.dump(jsondict, fout, indent = 4, default=str)

    return final_report


def _assembly_to_accession(assembly):
    assembly_prefix = path.basename(path.splitext(assembly)[0])
    if assembly_prefix.endswith("_genomic"):
        assembly_prefix = assembly_prefix[:-len("_genomic")]
    return assembly_prefix


def _write_report_subreports(args, clustering_results, snp_profile=None, lineage_profile=None):
    subreports = {
        "clustering": {
            "path": path.abspath(path.join(args.reportsDir, f"{args.output_prefix}.clustering_results.json")),
            "record_count": len(clustering_results),
        },
        "snp_profile": {
            "path": path.abspath(path.join(args.reportsDir, f"{args.output_prefix}.snp_profile.json")),
            "record_count": len(snp_profile or {}),
        },
        "lineage_profile": {
            "path": path.abspath(path.join(args.reportsDir, f"{args.output_prefix}.lineage_profile.json")),
            "record_count": len(lineage_profile or {}),
        },
    }
    if hasattr(args, "kronaWDir"):
        krona_report = path.join(args.kronaWDir, f"{args.output_prefix}.krona.html")
        subreports["krona"] = {
            "path": path.abspath(krona_report),
            "record_count": None,
        }

    payloads = {
        "clustering": clustering_results,
        "snp_profile": snp_profile or {},
        "lineage_profile": lineage_profile or {},
    }

    for key, payload in payloads.items():
        with open(subreports[key]["path"], "w") as fout:
            json.dump(payload, fout, indent=4, default=str)

    return subreports


def _index_detection_event_by_taxid(by_taxid, event, snp_profile=None, lineage_profile=None):
    taxid = str(event["taxon_id"])
    record = by_taxid.setdefault(
        taxid,
        {
            "taxon_id": taxid,
            "name": event["name"],
            "taxonomic_assignment": None,
            "snp_profile": {},
            "lineage_profile": {},
        },
    )
    record["name"] = event["name"]
    record["taxonomic_assignment"] = event
    _attach_profiles_for_event(record, event, snp_profile=snp_profile, lineage_profile=lineage_profile)


def _attach_profiles_for_event(record, event, snp_profile=None, lineage_profile=None):
    accessions = _event_accession_keys(event)
    if isinstance(event.get("closest_variant"), dict):
        accessions.update(_event_accession_keys(event["closest_variant"]))
    for accession in accessions:
        if snp_profile and accession in snp_profile:
            record["snp_profile"][accession] = snp_profile[accession]
        if lineage_profile and accession in lineage_profile:
            record["lineage_profile"][accession] = lineage_profile[accession]


def _event_accession_keys(event):
    keys = set()
    for assembly_key in ("assembly",):
        if event.get(assembly_key):
            keys.add(_assembly_to_accession(event[assembly_key]))
    if event.get("accession"):
        keys.add(str(event["accession"]))
    return keys

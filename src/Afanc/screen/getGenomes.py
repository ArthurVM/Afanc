import shutil
import json
from os import path


def getLocalGenomes(out_json, args):
    """ if args.fetch_assemblies is False, get genomes from the autodatabase results directory:
            autodb_results/selectFasta_autoDatabase_cleanFasta/
    """
    with open(args.db_fastas_json, "r") as fin:
        assembly_dict = json.load(fin)

    ## read in species report json to get hit genomes for mapping to
    with open(out_json, 'r') as fin:
        json_data = json.load(fin)

        for subdict in json_data["Detection_events"]:

            if _select_profile_reference(subdict, args):
                continue

            mapping_event = _mapping_target_event(subdict)
            taxID = str(mapping_event["taxon_id"])
            species_ID = mapping_event["name"].replace(" ", "_")

            ## copy assembly to the mapping working directory for index construction
            assembly_path = _resolve_database_reference(taxID, species_ID, assembly_dict, args.cleanFasta)
            if assembly_path is None:
                mapping_event.setdefault("warnings", []).append(
                    f"No local reference assembly found for taxID {taxID} or species key {species_ID}."
                )
                continue

            print(f"Copying {assembly_path}...")
            shutil.copy(assembly_path, "./")

            mapping_event["assembly"] = path.basename(assembly_path)

    with open(out_json, "w") as fout:
        json.dump(json_data, fout, indent=4, default=str)


def _select_profile_reference(event, args):
    from .profiles import profile_paths, resolve_profile_for_event, validate_profile_reference

    if not getattr(args, "profiles_dir", None) or not getattr(args, "profiles_manifest", None):
        return False

    lineage_profile = resolve_profile_for_event(event, args.profiles_manifest)
    if lineage_profile is None:
        return False

    validation = validate_profile_reference(lineage_profile, args.profiles_dir)
    event["lineage_profile_validation"] = validation
    if not validation["valid"]:
        event["lineage_profile"] = {
            "taxon_id": lineage_profile["taxon_id"],
            "name": lineage_profile["name"],
            "enabled": False,
            "validation_errors": validation["errors"],
        }
        return False

    resolved_paths = profile_paths(lineage_profile, args.profiles_dir)
    assembly_path = resolved_paths["reference"]
    print(f"Copying lineage profile reference {assembly_path}...")
    shutil.copy(assembly_path, "./")

    assembly_basename = path.basename(assembly_path)
    accession = _accession_from_assembly(assembly_basename)
    event["assembly"] = assembly_basename
    if isinstance(event.get("closest_variant"), dict):
        event["closest_variant"]["assembly"] = assembly_basename
    event["lineage_profile"] = {
        "taxon_id": lineage_profile["taxon_id"],
        "name": lineage_profile["name"],
        "reference": resolved_paths["reference"],
        "model": resolved_paths["model"],
        "allele_id_format": lineage_profile.get("allele_id_format", "{chrom}.{start}.{ref}.{alt}"),
        "profile_match": lineage_profile.get("profile_match"),
    }
    args.lineage_profiles_by_accession[accession] = event["lineage_profile"]
    return True


def _mapping_target_event(event):
    if "closest_variant" in event:
        return event["closest_variant"]
    return event


def _resolve_database_reference(taxID, species_ID, assembly_dict, cleanFasta):
    if taxID in assembly_dict:
        return path.join(cleanFasta, assembly_dict[taxID][0])
    if species_ID in assembly_dict:
        return path.join(cleanFasta, assembly_dict[species_ID][0])
    return None

def _accession_from_assembly(assembly):
    tmp_acc = path.basename(path.splitext(assembly)[0])
    if tmp_acc.endswith("_genomic"):
        return tmp_acc.split("_genomic")[0]
    return tmp_acc

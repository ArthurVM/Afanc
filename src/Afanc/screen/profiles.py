from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from Bio import SeqIO


def load_profiles_manifest(database_dir: str, profiles_dir: Optional[str] = None) -> tuple[Optional[str], Dict[str, Any]]:
    candidate_dirs = []
    if profiles_dir is not None:
        candidate_dirs.append(Path(profiles_dir))

    database_path = Path(database_dir)
    candidate_dirs.extend(
        [
            database_path / "profiles",
            database_path.parent / "profiles",
        ]
    )

    selected_profiles_dir = None
    manifest_path = None
    for candidate_dir in candidate_dirs:
        candidate_manifest = candidate_dir / "profiles.json"
        if candidate_manifest.is_file():
            selected_profiles_dir = candidate_dir
            manifest_path = candidate_manifest
            break

    if manifest_path is None:
        return None, {"profiles": {}}

    with manifest_path.open("r") as handle:
        manifest = json.load(handle)
    return str(selected_profiles_dir), manifest


def resolve_profile_for_event(event: Mapping[str, Any], manifest: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    profiles = manifest.get("profiles") or {}
    candidate_taxids = _event_candidate_taxids(event)

    ## prefer exact taxid matches to aliases
    for candidate_taxid in candidate_taxids:
        profile = profiles.get(candidate_taxid)
        if profile and profile.get("enabled", True):
            return _resolved_profile(profile, profile_taxid=candidate_taxid, match_type="taxid", matched_value=candidate_taxid)

    name_candidates = _event_candidate_names(event)
    exact_profile = _resolve_profile_by_name(profiles, name_candidates, allow_contains=False)
    if exact_profile is not None:
        return exact_profile

    return _resolve_profile_by_name(profiles, name_candidates, allow_contains=True)


def resolve_profile_for_species(species_name: str, manifest: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    """Resolve an enabled profile from a user-provided species/profile name."""
    profiles = manifest.get("profiles") or {}
    name_candidates = [str(species_name)]
    exact_profile = _resolve_profile_by_name(profiles, name_candidates, allow_contains=False)
    if exact_profile is not None:
        return exact_profile
    return _resolve_profile_by_name(profiles, name_candidates, allow_contains=True)


def _event_candidate_taxids(event: Mapping[str, Any]) -> list[str]:
    candidate_taxids = []

    if "taxon_id" in event:
        candidate_taxids.append(str(event["taxon_id"]))
    if isinstance(event.get("closest_variant"), Mapping) and "taxon_id" in event["closest_variant"]:
        candidate_taxids.append(str(event["closest_variant"]["taxon_id"]))

    return candidate_taxids


def _event_candidate_names(event: Mapping[str, Any]) -> list[str]:
    candidate_names = []

    if event.get("name"):
        candidate_names.append(str(event["name"]))
    if isinstance(event.get("closest_variant"), Mapping) and event["closest_variant"].get("name"):
        candidate_names.append(str(event["closest_variant"]["name"]))

    return candidate_names


def _resolve_profile_by_name(
    profiles: Mapping[str, Mapping[str, Any]],
    name_candidates: list[str],
    allow_contains: bool,
) -> Optional[Dict[str, Any]]:
    normalised_candidates = {
        name: _normalise_profile_name(name)
        for name in name_candidates
        if name
    }

    for profile_taxid, profile in profiles.items():
        if not profile or not profile.get("enabled", True):
            continue

        profile_names = _profile_match_names(profile)
        for raw_candidate, normalised_candidate in normalised_candidates.items():
            for profile_name in profile_names:
                if normalised_candidate == profile_name:
                    return _resolved_profile(
                        profile,
                        profile_taxid=profile_taxid,
                        match_type="name",
                        matched_value=raw_candidate,
                    )

                if allow_contains and _safe_contains_match(normalised_candidate, profile_name):
                    return _resolved_profile(
                        profile,
                        profile_taxid=profile_taxid,
                        match_type="alias_contains",
                        matched_value=raw_candidate,
                    )

    return None


def _resolved_profile(profile: Mapping[str, Any], profile_taxid: str, match_type: str, matched_value: str) -> Dict[str, Any]:
    resolved = dict(profile)
    resolved["taxon_id"] = str(profile_taxid)
    resolved["profile_taxon_id"] = str(profile_taxid)
    resolved["profile_match"] = {
        "type": match_type,
        "matched_value": str(matched_value),
    }
    return resolved


def _profile_match_names(profile: Mapping[str, Any]) -> list[str]:
    names = []
    if profile.get("name"):
        names.append(str(profile["name"]))
    names.extend(str(alias) for alias in profile.get("aliases", []) if alias)
    return [_normalise_profile_name(name) for name in names if _normalise_profile_name(name)]


def _normalise_profile_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def _safe_contains_match(candidate: str, profile_name: str) -> bool:
    if not candidate or not profile_name:
        return False

    ## ignore short aliases during substring matching
    if len(profile_name) < 7:
        return False

    return profile_name in candidate


def validate_profile_reference(profile: Mapping[str, Any], profiles_dir: str) -> Dict[str, Any]:
    profiles_path = Path(profiles_dir)
    reference_path = profiles_path / profile["reference"]
    model_path = profiles_path / profile["model"]

    validation = {
        "valid": True,
        "reference": str(reference_path),
        "model": str(model_path),
        "checked_loci": 0,
        "errors": [],
    }

    if not reference_path.is_file():
        validation["valid"] = False
        validation["errors"].append(f"Profile reference not found: {reference_path}")
    if not model_path.is_file():
        validation["valid"] = False
        validation["errors"].append(f"Profile model not found: {model_path}")
    if not validation["valid"]:
        return validation

    with model_path.open("r") as handle:
        model = json.load(handle)

    reference_sequences = _read_fasta(reference_path)
    for locus in model.get("loci", []):
        chrom = str(locus.get("chrom", ""))
        pos = int(locus.get("pos", 0))
        expected_ref = locus.get("ref")
        if expected_ref is None:
            continue
        expected_ref = str(expected_ref).upper()

        if chrom not in reference_sequences:
            validation["valid"] = False
            validation["errors"].append(
                f"Model locus {locus.get('locus_id', locus.get('allele_id', '?'))} uses missing contig {chrom!r}."
            )
            continue
        if pos < 1 or pos > len(reference_sequences[chrom]):
            validation["valid"] = False
            validation["errors"].append(
                f"Model locus {locus.get('locus_id', locus.get('allele_id', '?'))} position {pos} is outside contig {chrom}."
            )
            continue

        observed_ref = reference_sequences[chrom][pos - 1].upper()
        validation["checked_loci"] += 1
        if observed_ref != expected_ref:
            validation["valid"] = False
            validation["errors"].append(
                f"Model locus {locus.get('locus_id', locus.get('allele_id', '?'))} expects {expected_ref} "
                f"at {chrom}:{pos}, but reference has {observed_ref}."
            )

    return validation


def profile_paths(profile: Mapping[str, Any], profiles_dir: str) -> Dict[str, str]:
    profiles_path = Path(profiles_dir)
    return {
        "reference": str(profiles_path / profile["reference"]),
        "model": str(profiles_path / profile["model"]),
    }


def _read_fasta(fasta_path: Path) -> Dict[str, str]:
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    sequences = {}
    with opener(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[str(record.id)] = str(record.seq)
    return sequences

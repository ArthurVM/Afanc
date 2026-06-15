import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def _safe_prob(value, eps=1e-12):
    return min(max(float(value), eps), 1.0 - eps)


def load_json(json_path):
    path = Path(json_path)
    with path.open() as handle:
        return json.load(handle)


def load_model_json(model_json):
    return load_json(model_json)


def load_sample_alleles_json(sample_alleles_json):
    return load_json(sample_alleles_json)


def _validate_model(model):
    if model.get("schema_version") not in {"minimal_empirical_loso_v1", "model_schema_v1"}:
        raise ValueError("Unsupported model schema_version.")

    architecture = model.get("architecture", {})
    if architecture.get("feature_family") != "empirical":
        raise ValueError(
            "Bayesian classification requires an empirical-emission model. "
            "Use canonical_as_empirical/CEB for canonical SNP models."
        )
    if architecture.get("topology") not in {"flat", "hierarchical"}:
        raise ValueError("Model topology must be 'flat' or 'hierarchical'.")
    if architecture.get("emission_model") not in {"empirical_bayes", "semi_bayes"}:
        raise ValueError("Model emission_model must be 'empirical_bayes' or 'semi_bayes'.")

    family = _model_family(model)
    if architecture.get("topology") == "hierarchical" and family != "CEB":
        raise ValueError("Hierarchical Bayesian classification is currently supported only for CEB models.")


def _model_family(model):
    model_type = str(model.get("model_type", ""))
    architecture = model.get("architecture", {})
    emission_model = str(architecture.get("emission_model", ""))

    if model_type == "canonical_as_empirical":
        if emission_model != "semi_bayes":
            raise ValueError("CEB models must use emission_model='semi_bayes'.")
        return "CEB"
    if model_type == "empirical_loso" and emission_model == "empirical_bayes":
        return "EB"
    if model_type == "empirical_loso" and emission_model == "semi_bayes":
        return "SB"
    if emission_model == "empirical_bayes":
        return "EB"
    if emission_model == "semi_bayes":
        return "SB"
    raise ValueError(f"Unsupported Bayesian model family for model_type={model_type!r}.")


def _normalise_sample_alleles(sample_alleles_json):
    normalised = {}
    for sample_id, allele_ids in sample_alleles_json.items():
        if not isinstance(allele_ids, list):
            raise ValueError(f"Sample {sample_id!r} must map to a list of allele IDs.")
        normalised[str(sample_id)] = {str(allele_id) for allele_id in allele_ids}
    return normalised


def _parse_allele_position(allele_id):
    parts = str(allele_id).rsplit(".", 3)
    if len(parts) != 4:
        raise ValueError(
            f"Allele ID {allele_id!r} must have format 'chrom.pos.ref.alt'."
        )
    chrom, pos, _ref, _alt = parts
    return str(chrom), int(pos)


def _normalise_missing_positions(missing_positions):
    normalised = set()
    for item in missing_positions or []:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise ValueError("Missing positions must be [chrom, pos] pairs.")
        chrom, pos = item
        normalised.add((str(chrom), int(pos)))
    return normalised


def _validate_sample_payload(payload):
    if not isinstance(payload, dict):
        raise ValueError("Sample input must be a dict with 'alleles' and 'missing'.")
    if set(payload) != {"alleles", "missing"}:
        raise ValueError("Sample input must contain exactly 'alleles' and 'missing'.")
    if not isinstance(payload["alleles"], list):
        raise ValueError("Sample input 'alleles' must be a list.")
    if not isinstance(payload["missing"], list):
        raise ValueError("Sample input 'missing' must be a list.")
    _normalise_missing_positions(payload["missing"])


def _normalise_sample_inputs(sample_alleles_json, *, sample_id="sample"):
    if isinstance(sample_alleles_json, dict) and set(sample_alleles_json) == {"alleles", "missing"}:
        _validate_sample_payload(sample_alleles_json)
        return {
            str(sample_id): {
                "alleles": sample_alleles_json["alleles"],
                "missing": sample_alleles_json["missing"],
                "_has_missing_evidence": True,
            }
        }

    sample_alleles = _normalise_sample_alleles(sample_alleles_json)
    return {
        str(sample_id): {
            "alleles": sorted(allele_set),
            "missing": [],
            "_has_missing_evidence": False,
        }
        for sample_id, allele_set in sample_alleles.items()
    }


def _effective_sample_evidence(payload, model, *, infer_reference_allele_markers=True):
    allele_set = {str(allele_id) for allele_id in payload["alleles"]}
    missing_positions = _normalise_missing_positions(payload.get("missing", []))

    inferred_locus_ids = []
    can_infer_reference_markers = (
        infer_reference_allele_markers
        and model.get("model_type") == "canonical_as_empirical"
        and payload.get("_has_missing_evidence", False)
    )
    if can_infer_reference_markers:
        observed_positions = {
            _parse_allele_position(allele_id)
            for allele_id in allele_set
        }
        for locus in model.get("loci", []):
            if str(locus.get("ref")) != str(locus.get("alt")):
                continue

            position = (str(locus["chrom"]), int(locus["pos"]))
            if position in missing_positions:
                continue
            if position in observed_positions:
                continue

            locus_id = str(locus.get("locus_id", locus["allele_id"]))
            allele_set.add(locus_id)
            inferred_locus_ids.append(locus_id)

    return {
        "allele_set": allele_set,
        "missing_positions": missing_positions,
        "reference_marker_inference": {
            "enabled": bool(can_infer_reference_markers),
            "inferred_locus_count": len(inferred_locus_ids),
            "inferred_locus_ids": sorted(inferred_locus_ids),
        },
    }


def _get_emission_block(locus):
    emission = locus.get("emission", locus.get("empirical_bayes"))
    if emission is None:
        raise ValueError(f"Locus {locus.get('locus_id', '?')!r} is missing an emission block.")
    return emission


def _lineages_by_id(model):
    return {str(record["lineage_id"]): record for record in model["lineages"]}


def _loci_by_id(model):
    return {str(locus.get("locus_id", locus["allele_id"])): locus for locus in model["loci"]}


def _default_likelihood_mode(model):
    if model.get("model_type") == "canonical_as_empirical":
        return "fraction"
    return "product"


def _combine_log_terms(log_prior, target_terms, background_terms, likelihood_mode):
    if likelihood_mode == "product":
        return log_prior + sum(target_terms) + sum(background_terms)
    if likelihood_mode == "fraction":
        logp = log_prior
        if target_terms:
            logp += sum(target_terms) / len(target_terms)
        if background_terms:
            logp += sum(background_terms) / len(background_terms)
        return logp
    raise ValueError("likelihood_mode must be 'product' or 'fraction'.")


def score_samples(sample_alleles, model, *, likelihood_mode=None):
    likelihood_mode = likelihood_mode or _default_likelihood_mode(model)
    locus_ids = [str(locus.get("locus_id", locus["allele_id"])) for locus in model["loci"]]

    records = {}
    for sample_id, allele_set in sample_alleles.items():
        scores = _score_candidate_lineages(
            model["lineages"],
            {str(allele_id) for allele_id in allele_set},
            model,
            locus_ids,
            likelihood_mode=likelihood_mode,
            missing_positions=set(),
        )
        records[str(sample_id)] = {
            item["lineage_id"]: item["log_score"]
            for item in scores
        }

    score_df = pd.DataFrame.from_dict(records, orient="index")
    max_scores = score_df.max(axis=1)
    prob_df = np.exp(score_df.sub(max_scores, axis=0))
    prob_df = prob_df.div(prob_df.sum(axis=1), axis=0)
    return score_df, prob_df


def _normalise_log_scores(log_scores):
    max_score = max(log_scores.values())
    exp_scores = {
        lineage: math.exp(score - max_score)
        for lineage, score in log_scores.items()
    }
    total = sum(exp_scores.values())
    return {
        lineage: value / total
        for lineage, value in exp_scores.items()
    }


def _score_candidate_lineages(
    candidate_lineages,
    allele_set,
    model,
    locus_ids,
    *,
    likelihood_mode=None,
    missing_positions=None,
):
    likelihood_mode = likelihood_mode or _default_likelihood_mode(model)
    loci_by_id = _loci_by_id(model)
    missing_positions = set(missing_positions or ())
    candidate_ids = [str(record["lineage_id"]) for record in candidate_lineages]
    priors = {
        str(record["lineage_id"]): float(record.get("prior", 0.0) or 0.0)
        for record in candidate_lineages
    }
    prior_sum = sum(value for value in priors.values() if value > 0.0)
    if prior_sum <= 0.0:
        priors = {lineage: 1.0 / max(len(candidate_ids), 1) for lineage in candidate_ids}
    else:
        priors = {
            lineage: (priors[lineage] / prior_sum if priors[lineage] > 0.0 else 0.0)
            for lineage in candidate_ids
        }

    locus_ids = sorted({str(locus_id) for locus_id in locus_ids})
    log_scores = {}
    metrics = {}
    for lineage in candidate_ids:
        log_prior = math.log(_safe_prob(priors[lineage]))
        target_terms = []
        background_terms = []
        matched_locus_ids = []
        missing_locus_ids = []
        absent_locus_ids = []
        for locus_id in locus_ids:
            locus = loci_by_id[locus_id]
            position = (str(locus["chrom"]), int(locus["pos"]))
            is_missing = position in missing_positions
            is_target = lineage == str(locus["target_lineage"])
            if is_missing:
                if is_target:
                    missing_locus_ids.append(locus_id)
                continue

            emission = _get_emission_block(locus)
            p = emission["target_frequency"]
            if not is_target:
                p = emission["background_frequency"]
            p = _safe_prob(p)
            allele_id = str(locus["allele_id"])
            is_present = allele_id in allele_set
            if is_present and is_target:
                matched_locus_ids.append(locus_id)
            if not is_present and is_target:
                absent_locus_ids.append(locus_id)
            term = math.log(p) if is_present else math.log(1.0 - p)
            if is_target:
                target_terms.append(term)
            else:
                background_terms.append(term)
        log_scores[lineage] = _combine_log_terms(
            log_prior,
            target_terms,
            background_terms,
            likelihood_mode,
        )
        configured_count = len([
            locus_id for locus_id in locus_ids
            if str(loci_by_id[locus_id]["target_lineage"]) == lineage
        ])
        missing_count = len(missing_locus_ids)
        callable_count = configured_count - missing_count
        matched_count = len(matched_locus_ids)
        absent_count = len(absent_locus_ids)
        callable_support_fraction = matched_count / callable_count if callable_count else 0.0
        callable_fraction = callable_count / configured_count if configured_count else 0.0
        metrics[lineage] = {
            "matched_loci": matched_count,
            "callable_loci": callable_count,
            "active_loci": callable_count,
            "configured_loci": configured_count,
            "missing_loci": missing_count,
            "absent_loci": absent_count,
            "callable_support_fraction": callable_support_fraction,
            "support_fraction": callable_support_fraction,
            "callable_fraction": callable_fraction,
            "decision_score": callable_support_fraction * callable_fraction,
            "matched_locus_ids": sorted(matched_locus_ids),
            "missing_locus_ids": sorted(missing_locus_ids),
            "absent_locus_ids": sorted(absent_locus_ids),
        }

    posteriors = _normalise_log_scores(log_scores)
    records = []
    for lineage in candidate_ids:
        records.append(
            {
                "lineage_id": lineage,
                "posterior": float(posteriors[lineage]),
                "log_score": float(log_scores[lineage]),
                **metrics[lineage],
            }
        )
    records.sort(
        key=lambda item: (
            item["posterior"],
            item["matched_loci"],
            item["callable_support_fraction"],
            item["callable_fraction"],
            item["lineage_id"],
        ),
        reverse=True,
    )
    return records


def _tie_summary(sample_probs, tie_delta):
    if sample_probs.empty:
        return {
            "tie_delta": float(tie_delta),
            "is_tie": False,
            "tied_lineage_count": 0,
            "tied_lineages": [],
            "posterior_margin": None,
        }

    best_posterior = float(sample_probs.iloc[0])
    tied = sample_probs[
        sample_probs.map(lambda posterior: best_posterior - float(posterior) <= float(tie_delta))
    ]
    second_posterior = float(sample_probs.iloc[1]) if len(sample_probs) > 1 else None

    return {
        "tie_delta": float(tie_delta),
        "is_tie": len(tied) > 1,
        "tied_lineage_count": int(len(tied)),
        "tied_lineages": [str(lineage_id) for lineage_id in tied.index],
        "posterior_margin": (
            best_posterior - second_posterior
            if second_posterior is not None
            else None
        ),
    }


def _tie_summary_from_records(lineage_scores, tie_delta):
    if not lineage_scores:
        return {
            "tie_delta": float(tie_delta),
            "is_tie": False,
            "tied_lineage_count": 0,
            "tied_lineages": [],
            "posterior_margin": None,
        }

    best_posterior = float(lineage_scores[0]["posterior"])
    tied = [
        item
        for item in lineage_scores
        if best_posterior - float(item["posterior"]) <= float(tie_delta)
    ]
    second_posterior = (
        float(lineage_scores[1]["posterior"])
        if len(lineage_scores) > 1
        else None
    )
    return {
        "tie_delta": float(tie_delta),
        "is_tie": len(tied) > 1,
        "tied_lineage_count": len(tied),
        "tied_lineages": [str(item["lineage_id"]) for item in tied],
        "posterior_margin": (
            best_posterior - second_posterior
            if second_posterior is not None
            else None
        ),
    }


def _format_record_scores(lineage_scores):
    return [
        {
            "lineage_id": item["lineage_id"],
            "posterior": float(item["posterior"]),
            "log_score": float(item["log_score"]),
            "matched_loci": int(item["matched_loci"]),
            "callable_loci": int(item.get("callable_loci", item["active_loci"])),
            "active_loci": int(item["active_loci"]),
            "configured_loci": int(item.get("configured_loci", item["active_loci"])),
            "missing_loci": int(item.get("missing_loci", 0)),
            "absent_loci": int(item.get("absent_loci", 0)),
            "callable_support_fraction": float(item.get("callable_support_fraction", item["support_fraction"])),
            "support_fraction": float(item["support_fraction"]),
            "callable_fraction": float(item.get("callable_fraction", 0.0)),
            "decision_score": float(item.get("decision_score", 0.0)),
            "state": item.get("state"),
        }
        for item in lineage_scores
    ]


def _closest_lineage_fields(best_score):
    if best_score is None:
        return {
            "closest_lineage": None,
            "closest_lineage_posterior": None,
            "closest_lineage_state": None,
            "closest_lineage_decision_score": None,
            "closest_lineage_matched_loci": None,
            "closest_lineage_callable_loci": None,
            "closest_lineage_configured_loci": None,
            "closest_lineage_missing_loci": None,
            "closest_lineage_absent_loci": None,
        }

    if int(best_score.get("matched_loci", 0)) <= 0:
        return {
            "closest_lineage": None,
            "closest_lineage_posterior": None,
            "closest_lineage_state": None,
            "closest_lineage_decision_score": None,
            "closest_lineage_matched_loci": None,
            "closest_lineage_callable_loci": None,
            "closest_lineage_configured_loci": None,
            "closest_lineage_missing_loci": None,
            "closest_lineage_absent_loci": None,
        }

    return {
        "closest_lineage": str(best_score["lineage_id"]),
        "closest_lineage_posterior": float(best_score["posterior"]),
        "closest_lineage_state": best_score.get("state"),
        "closest_lineage_decision_score": float(best_score.get("decision_score", 0.0)),
        "closest_lineage_matched_loci": int(best_score["matched_loci"]),
        "closest_lineage_callable_loci": int(best_score.get("callable_loci", best_score["active_loci"])),
        "closest_lineage_configured_loci": int(best_score.get("configured_loci", best_score["active_loci"])),
        "closest_lineage_missing_loci": int(best_score.get("missing_loci", 0)),
        "closest_lineage_absent_loci": int(best_score.get("absent_loci", 0)),
    }


def _classify_hierarchical_sample(
    sample_id,
    allele_set,
    model,
    *,
    min_support=1,
    min_support_fraction=0.8,
    min_callable_fraction=0.5,
    ambiguity_score_margin=0.1,
    allow_incomplete_descent=True,
    likelihood_mode=None,
    tie_delta=0.0,
    missing_positions=None,
    extra=None,
):
    likelihood_mode = likelihood_mode or _default_likelihood_mode(model)
    missing_positions = set(missing_positions or ())
    hierarchy = model["hierarchy"]
    children_by_parent = {
        str(parent): [str(child) for child in children]
        for parent, children in hierarchy.get("children_by_parent", {}).items()
    }
    lineages_by_id = _lineages_by_id(model)
    start = str(hierarchy.get("top_classifier_node") or "__root__")
    state_priority = {
        "supported": 4,
        "possible_incomplete": 3,
        "unknown_all_missing": 2,
        "no_markers": 1,
        "rejected": 0,
    }

    def score_children(current):
        children = [child for child in children_by_parent.get(current, []) if child in lineages_by_id]
        if not children:
            return [], []

        candidate_lineages = [lineages_by_id[child] for child in children]
        node_locus_ids = []
        for child in children:
            node_locus_ids.extend(lineages_by_id[child].get("direct_locus_ids", []))
        return children, _score_candidate_lineages(
            candidate_lineages,
            allele_set,
            model,
            node_locus_ids,
            likelihood_mode=likelihood_mode,
            missing_positions=missing_positions,
        )

    def annotate_scores(node_scores):
        annotated = []
        for item in node_scores:
            record = dict(item)
            configured = int(record.get("configured_loci", 0))
            callable_loci = int(record.get("callable_loci", record.get("active_loci", 0)))
            matched = int(record.get("matched_loci", 0))
            support_fraction = float(record.get("callable_support_fraction", record.get("support_fraction", 0.0)))
            callable_fraction = float(record.get("callable_fraction", 0.0))

            if configured == 0:
                state = "no_markers"
            elif callable_loci == 0:
                state = "unknown_all_missing"
            elif matched >= int(min_support) and support_fraction >= float(min_support_fraction):
                if callable_fraction >= float(min_callable_fraction):
                    state = "supported"
                else:
                    state = "possible_incomplete"
            else:
                state = "rejected"

            record["state"] = state
            record["state_priority"] = state_priority[state]
            record["decision_score"] = support_fraction * callable_fraction
            annotated.append(record)

        annotated.sort(
            key=lambda item: (
                item["state_priority"],
                item["matched_loci"],
                item["callable_support_fraction"],
                item["callable_fraction"],
                item["posterior"],
                item["lineage_id"],
            ),
            reverse=True,
        )
        return annotated

    def compatible_alternatives(node_scores, selected_child):
        return [
            {
                "lineage_id": item["lineage_id"],
                "state": item["state"],
                "matched_loci": int(item["matched_loci"]),
                "callable_loci": int(item["callable_loci"]),
                "configured_loci": int(item["configured_loci"]),
                "missing_loci": int(item["missing_loci"]),
                "callable_support_fraction": float(item["callable_support_fraction"]),
                "callable_fraction": float(item["callable_fraction"]),
                "decision_score": float(item["decision_score"]),
            }
            for item in node_scores
            if item["lineage_id"] != selected_child
            and item["state"] in {"supported", "possible_incomplete"}
        ]

    def trace_record(current, node_scores, *, decision, selected_child=None, reason=None):
        best = node_scores[0]
        node_tie_summary = _tie_summary_from_records(node_scores, tie_delta)
        record = {
            "node": current,
            "best_child": best["lineage_id"],
            "decision": decision,
            "selected_child": selected_child,
            "reason": reason,
            "matched_loci": int(best["matched_loci"]),
            "callable_loci": int(best.get("callable_loci", best["active_loci"])),
            "active_loci": int(best["active_loci"]),
            "configured_loci": int(best.get("configured_loci", best["active_loci"])),
            "missing_loci": int(best.get("missing_loci", 0)),
            "absent_loci": int(best.get("absent_loci", 0)),
            "callable_support_fraction": float(best.get("callable_support_fraction", best["support_fraction"])),
            "support_fraction": float(best["support_fraction"]),
            "callable_fraction": float(best.get("callable_fraction", 0.0)),
            "state": best.get("state"),
            "tie_summary": node_tie_summary,
            "alternative_compatible_children": compatible_alternatives(node_scores, selected_child),
            "children": _format_record_scores(node_scores),
        }
        return record

    def decide_node(current, node_scores):
        supported = [item for item in node_scores if item["state"] == "supported"]
        possible = [item for item in node_scores if item["state"] == "possible_incomplete"]
        unknown = [
            item
            for item in node_scores
            if item["state"] == "unknown_all_missing"
            and int(item.get("configured_loci", 0)) > 0
            and children_by_parent.get(item["lineage_id"])
        ]
        transparent = [
            item
            for item in node_scores
            if item["state"] == "no_markers"
            and int(item.get("configured_loci", 0)) == 0
            and children_by_parent.get(item["lineage_id"])
        ]

        if len(supported) == 1:
            return {
                "action": "descend",
                "child": supported[0]["lineage_id"],
                "status": "supported",
                "reason": "single_supported_child",
            }

        if len(supported) > 1:
            best, second = supported[0], supported[1]
            margin = float(best["decision_score"]) - float(second["decision_score"])
            if margin >= float(ambiguity_score_margin):
                return {
                    "action": "descend",
                    "child": best["lineage_id"],
                    "status": "supported",
                    "reason": "best_supported_child_clear_margin",
                }
            return {
                "action": "stop",
                "status": "ambiguous",
                "reason": "ambiguous_supported_children",
            }

        if len(possible) == 1 and allow_incomplete_descent:
            return {
                "action": "descend",
                "child": possible[0]["lineage_id"],
                "status": "incomplete",
                "reason": "single_possible_incomplete_child",
            }

        if len(possible) == 1:
            return {
                "action": "stop",
                "status": "stopped_at_parent",
                "reason": "possible_incomplete_child_below_confidence_threshold",
            }

        if len(possible) > 1:
            return {
                "action": "stop",
                "status": "ambiguous",
                "reason": "ambiguous_incomplete_children",
            }

        if unknown:
            return {
                "action": "try_unknown_descendants",
                "children": [item["lineage_id"] for item in unknown],
                "status": "incomplete",
                "reason": "no_callable_child_evidence",
            }

        if transparent:
            return {
                "action": "try_transparent_descendants",
                "children": [item["lineage_id"] for item in transparent],
                "status": "supported",
                "reason": "markerless_child_with_descendants",
            }

        return {
            "action": "stop",
            "status": "stopped_at_parent",
            "reason": "no_supported_child",
        }

    def descend(current, called_lineage):
        children, node_scores = score_children(current)
        if not children:
            if current in lineages_by_id:
                return current, [], [], True, "supported"
            return called_lineage, [], [], False, "stopped_at_parent"

        node_scores = annotate_scores(node_scores)
        decision = decide_node(current, node_scores)

        if decision["action"] == "descend":
            next_called = decision["child"]
            final_call, child_trace, final_scores, _accepted, child_status = descend(next_called, next_called)
            status = "incomplete" if decision["status"] == "incomplete" else child_status
            return (
                final_call,
                [
                    trace_record(
                        current,
                        node_scores,
                        decision="descend",
                        selected_child=next_called,
                        reason=decision["reason"],
                    )
                ] + child_trace,
                final_scores or node_scores,
                True,
                status,
            )

        if decision["action"] == "try_unknown_descendants":
            successful = []
            for child in decision["children"]:
                final_call, child_trace, final_scores, accepted, child_status = descend(child, called_lineage)
                if accepted and final_call != called_lineage and final_scores:
                    successful.append((child, final_call, child_trace, final_scores, child_status))

            if len(successful) == 1:
                child, final_call, child_trace, final_scores, _child_status = successful[0]
                return (
                    final_call,
                    [
                        trace_record(
                            current,
                            node_scores,
                            decision="descend_incomplete",
                            selected_child=child,
                            reason="descendant_supported_through_uncallable_child",
                        )
                    ] + child_trace,
                    final_scores or node_scores,
                    True,
                    "incomplete",
                )

            if len(successful) > 1:
                stopped_call = current if current in lineages_by_id else called_lineage
                return (
                    stopped_call,
                    [
                        trace_record(
                            current,
                            node_scores,
                            decision="stop",
                            reason="ambiguous_descendants_through_uncallable_children",
                        )
                    ],
                    node_scores,
                    False,
                    "ambiguous",
                )

        if decision["action"] == "try_transparent_descendants":
            successful = []
            for child in decision["children"]:
                final_call, child_trace, final_scores, accepted, child_status = descend(child, called_lineage)
                if accepted and final_call != called_lineage and final_scores:
                    successful.append((child, final_call, child_trace, final_scores, child_status))

            if len(successful) == 1:
                child, final_call, child_trace, final_scores, child_status = successful[0]
                status = "supported" if child_status == "stopped_at_parent" else child_status
                return (
                    final_call,
                    [
                        trace_record(
                            current,
                            node_scores,
                            decision="descend_transparent",
                            selected_child=child,
                            reason="descendant_supported_through_markerless_child",
                        )
                    ] + child_trace,
                    final_scores or node_scores,
                    True,
                    status,
                )

            if len(successful) > 1:
                stopped_call = current if current in lineages_by_id else called_lineage
                return (
                    stopped_call,
                    [
                        trace_record(
                            current,
                            node_scores,
                            decision="stop",
                            reason="ambiguous_descendants_through_markerless_children",
                        )
                    ],
                    node_scores,
                    False,
                    "ambiguous",
                )

        stopped_call = current if current in lineages_by_id else "unclassified"
        return (
            stopped_call,
            [
                trace_record(
                    current,
                    node_scores,
                    decision="stop",
                    reason=decision["reason"],
                )
            ],
            node_scores,
            False,
            decision["status"],
        )

    initial_call = start if start in lineages_by_id else "unclassified"
    called_lineage, decision_trace, final_scores, _accepted, call_status = descend(start, initial_call)

    best_score = final_scores[0] if final_scores else None
    record = {
        "sample_id": str(sample_id),
        "classifier": "hierarchical_ceb",
        "model_family": _model_family(model),
        "posterior_scope": "local_sibling_nodes",
        "best_lineage": str(called_lineage),
        "best_posterior": float(best_score["posterior"]) if best_score else 0.0,
        "call_status": str(call_status),
        "tie_summary": _tie_summary_from_records(final_scores, tie_delta),
        "lineage_posteriors": _format_record_scores(final_scores),
        "topology": "hierarchical",
        "min_support": int(min_support),
        "min_support_fraction": float(min_support_fraction),
        "min_callable_fraction": float(min_callable_fraction),
        "ambiguity_score_margin": float(ambiguity_score_margin),
        "likelihood_mode": str(likelihood_mode),
        "decision_trace": decision_trace,
    }
    if best_score is not None:
        record.update(
            {
                "matched_loci": int(best_score["matched_loci"]),
                "callable_loci": int(best_score.get("callable_loci", best_score["active_loci"])),
                "active_loci": int(best_score["active_loci"]),
                "configured_loci": int(best_score.get("configured_loci", best_score["active_loci"])),
                "missing_loci": int(best_score.get("missing_loci", 0)),
                "absent_loci": int(best_score.get("absent_loci", 0)),
                "callable_support_fraction": float(best_score.get("callable_support_fraction", best_score["support_fraction"])),
                "support_fraction": float(best_score["support_fraction"]),
                "callable_fraction": float(best_score.get("callable_fraction", 0.0)),
                "matched_locus_ids": list(best_score["matched_locus_ids"]),
                "missing_locus_ids": list(best_score.get("missing_locus_ids", [])),
                "absent_locus_ids": list(best_score.get("absent_locus_ids", [])),
            }
        )
    if str(called_lineage) == "unclassified":
        record.update(_closest_lineage_fields(best_score))
    if extra:
        record.update(extra)
    return record


def _classify_flat_sample(
    sample_id,
    evidence,
    model,
    *,
    likelihood_mode=None,
    tie_delta=0.0,
    extra=None,
):
    likelihood_mode = likelihood_mode or _default_likelihood_mode(model)
    locus_ids = [str(locus.get("locus_id", locus["allele_id"])) for locus in model["loci"]]
    scores = _score_candidate_lineages(
        model["lineages"],
        evidence["allele_set"],
        model,
        locus_ids,
        likelihood_mode=likelihood_mode,
        missing_positions=evidence["missing_positions"],
    )
    best_score = scores[0] if scores else None
    record = {
        "sample_id": str(sample_id),
        "classifier": "flat_bayes",
        "model_family": _model_family(model),
        "posterior_scope": "global",
        "best_lineage": str(best_score["lineage_id"]) if best_score else "unclassified",
        "best_posterior": float(best_score["posterior"]) if best_score else 0.0,
        "tie_summary": _tie_summary_from_records(scores, tie_delta),
        "lineage_posteriors": _format_record_scores(scores),
        "topology": "flat",
        "likelihood_mode": str(likelihood_mode),
    }
    if best_score is not None:
        record.update(
            {
                "matched_loci": int(best_score["matched_loci"]),
                "callable_loci": int(best_score.get("callable_loci", best_score["active_loci"])),
                "active_loci": int(best_score["active_loci"]),
                "configured_loci": int(best_score.get("configured_loci", best_score["active_loci"])),
                "missing_loci": int(best_score.get("missing_loci", 0)),
                "absent_loci": int(best_score.get("absent_loci", 0)),
                "callable_support_fraction": float(best_score.get("callable_support_fraction", best_score["support_fraction"])),
                "support_fraction": float(best_score["support_fraction"]),
                "callable_fraction": float(best_score.get("callable_fraction", 0.0)),
                "matched_locus_ids": list(best_score["matched_locus_ids"]),
                "missing_locus_ids": list(best_score.get("missing_locus_ids", [])),
                "absent_locus_ids": list(best_score.get("absent_locus_ids", [])),
            }
        )
    if str(record["best_lineage"]) == "unclassified":
        record.update(_closest_lineage_fields(best_score))
    if extra:
        record.update(extra)
    return record


def classify_samples(
    sample_alleles_json,
    model,
    *,
    sample_id="sample",
    min_support=1,
    min_support_fraction=0.75,
    min_callable_fraction=0.5,
    ambiguity_score_margin=0.1,
    allow_incomplete_descent=True,
    infer_reference_allele_markers=True,
    likelihood_mode=None,
    tie_delta=0.0,
):
    _validate_model(model)
    likelihood_mode = likelihood_mode or _default_likelihood_mode(model)
    sample_payloads = _normalise_sample_inputs(sample_alleles_json, sample_id=sample_id)
    sample_evidence = {
        sid: _effective_sample_evidence(
            payload,
            model,
            infer_reference_allele_markers=infer_reference_allele_markers,
        )
        for sid, payload in sample_payloads.items()
    }

    if model["architecture"]["topology"] == "hierarchical":
        classifications = []
        for sid, evidence in sample_evidence.items():
            classifications.append(
                _classify_hierarchical_sample(
                    sid,
                    evidence["allele_set"],
                    model,
                    min_support=min_support,
                    min_support_fraction=min_support_fraction,
                    min_callable_fraction=min_callable_fraction,
                    ambiguity_score_margin=ambiguity_score_margin,
                    allow_incomplete_descent=allow_incomplete_descent,
                    likelihood_mode=likelihood_mode,
                    tie_delta=tie_delta,
                    missing_positions=evidence["missing_positions"],
                    extra={
                        "reference_marker_inference": evidence["reference_marker_inference"],
                    },
                )
            )
        return classifications

    classifications = []
    for sid, evidence in sample_evidence.items():
        classifications.append(
            _classify_flat_sample(
                sid,
                evidence,
                model,
                likelihood_mode=likelihood_mode,
                tie_delta=tie_delta,
                extra={
                    "reference_marker_inference": evidence["reference_marker_inference"],
                },
            )
        )
    return classifications


def write_classifications_json(classifications, output_json):
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(classifications, handle, indent=2)
    return output_path

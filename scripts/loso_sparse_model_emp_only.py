import json
import math
import re
from collections import Counter, defaultdict, namedtuple
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score, log_loss

try:
    from ardal.core.ArdalHeaderUtils import ArdalHeaderUtils
except ImportError:
    ArdalHeaderUtils = None


LosoFold = namedtuple("LosoFold", ["study", "test", "train"])
ARDAL_DOTTED_ALLELE_ID_FORMAT = "{chr}.{start}.{ref}.{alt}"
HIERARCHY_ROOT_NODE = "__root__"
DEFAULT_REPORTING_SUPPORT_EVIDENCE_BINS = (
    {"bin_id": "bin_0", "min_observed_selected_alleles": 0, "max_observed_selected_alleles": 9},
    {"bin_id": "bin_1", "min_observed_selected_alleles": 10, "max_observed_selected_alleles": 19},
    {"bin_id": "bin_2", "min_observed_selected_alleles": 20, "max_observed_selected_alleles": 29},
    {"bin_id": "bin_3", "min_observed_selected_alleles": 30, "max_observed_selected_alleles": 39},
    {"bin_id": "bin_4", "min_observed_selected_alleles": 40, "max_observed_selected_alleles": 49},
    {"bin_id": "bin_5", "min_observed_selected_alleles": 50, "max_observed_selected_alleles": None},
)


def _decode_dotted_allele_with_ardal(allele_id):
    if ArdalHeaderUtils is None:
        return None

    try:
        header_utils = ArdalHeaderUtils(
            headers={"guids": [], "alleles": [str(allele_id)]},
            meta={},
            allele_id_format=ARDAL_DOTTED_ALLELE_ID_FORMAT,
        )
        allele_positions = header_utils.get_allele_positions()
    except Exception:
        return None

    decoded = allele_positions.get(str(allele_id))
    if decoded is None:
        return None

    chrom, pos = decoded
    return str(chrom), int(pos)


def parse_ardal_allele(allele_id):
    allele_id = str(allele_id)

    ## Prefer Ardal-backed decoding for dotted IDs so contig / position parsing
    ## stays consistent with the allele headers used elsewhere in the stack.
    dotted_match = re.match(r"^(?P<chrom>.+)\.(?P<pos>\d+)\.(?P<ref>[A-ZN]+)\.(?P<alt>[A-ZN]+)$", allele_id)
    if dotted_match:
        ardal_decoded = _decode_dotted_allele_with_ardal(allele_id)
        if ardal_decoded is not None:
            chrom, pos = ardal_decoded
            ref = dotted_match.group("ref")
            alt = dotted_match.group("alt")
            return {
                "chrom": chrom,
                "pos": pos,
                "bed_start": pos - 1,
                "bed_end": pos,
                "ref": ref,
                "alt": alt,
            }

    patterns = [
        r"^(?P<chrom>.+):(?P<pos>\d+):(?P<ref>[A-ZN]+)>(?P<alt>[A-ZN]+)$",
        r"^(?P<chrom>.+)_(?P<pos>\d+)_(?P<ref>[A-ZN]+)_(?P<alt>[A-ZN]+)$",
        r"^(?P<chrom>.+)\.(?P<pos>\d+)\.(?P<ref>[A-ZN]+)\.(?P<alt>[A-ZN]+)$",
    ]

    for pattern in patterns:
        match = re.match(pattern, allele_id)
        if match:
            chrom = match.group("chrom")
            pos = int(match.group("pos"))
            ref = match.group("ref")
            alt = match.group("alt")
            return {
                "chrom": chrom,
                "pos": pos,
                "bed_start": pos - 1,
                "bed_end": pos,
                "ref": ref,
                "alt": alt,
            }

    return {
        "chrom": allele_id,
        "pos": 0,
        "bed_start": 0,
        "bed_end": 0,
        "ref": None,
        "alt": "N",
    }


def _safe_prob(p, eps=1e-12):
    return min(max(float(p), eps), 1.0 - eps)


def _to_builtin_number(value):
    if value is None:
        return None
    if isinstance(value, (np.integer, np.int64)):
        return int(value)
    if isinstance(value, (np.floating, np.float64)):
        return float(value)
    return value


def _metric_score(y_true, y_pred, prob_df, metric):
    if metric == "macro_f1":
        return f1_score(y_true, y_pred, average="macro")
    if metric == "balanced_accuracy":
        return balanced_accuracy_score(y_true, y_pred)
    if metric == "log_loss":
        labels = sorted(pd.unique(y_true))
        return -log_loss(y_true, prob_df[labels], labels=labels)
    raise ValueError(f"Unsupported metric {metric!r}.")


def build_sample_allele_sets(ard_subset, alleles, sample_ids, backend="auto"):
    guid_to_alleles = ard_subset.get.subset(
        alleles=alleles,
        chunk_size=1000,
        threads=12,
        drop_zero_cols=False,
    ).io.to_dict(backend=backend)

    sample_alleles = {}
    for sample_id in sample_ids:
        sample_alleles[sample_id] = set(guid_to_alleles.get(sample_id, []))

    return sample_alleles


def estimate_lineage_frequencies(sample_alleles, labels_by_sample, selected_alleles, smoothing_alpha=1.0):
    lineage_to_samples = defaultdict(set)
    for sample_id, lineage in labels_by_sample.items():
        lineage_to_samples[lineage].add(sample_id)

    all_samples = set(labels_by_sample)
    lineages = sorted(lineage_to_samples)

    lineage_freqs = {allele: {} for allele in selected_alleles}
    pooled_background = {}

    for allele in selected_alleles:
        carriers = {sample_id for sample_id, aset in sample_alleles.items() if allele in aset}

        ## Estimate one smoothed Bernoulli rate per lineage, then collapse all
        ## non-target lineages into a shared background emission.
        for lineage in lineages:
            lineage_samples = lineage_to_samples[lineage]
            n = len(lineage_samples)
            k = len(carriers & lineage_samples)
            freq = (k + smoothing_alpha) / (n + 2.0 * smoothing_alpha)
            lineage_freqs[allele][lineage] = float(freq)

        target_lineage = max(lineage_freqs[allele], key=lineage_freqs[allele].get)
        background_samples = all_samples - lineage_to_samples[target_lineage]
        bg_n = len(background_samples)
        bg_k = len(carriers & background_samples)
        pooled_background[allele] = float((bg_k + smoothing_alpha) / (bg_n + 2.0 * smoothing_alpha))

    return lineages, lineage_to_samples, lineage_freqs, pooled_background


def score_samples_sparse_tables(sample_alleles, selected_alleles, lineages, priors, lineage_freqs, pooled_background):
    records = {}

    for sample_id, aset in sample_alleles.items():
        log_scores = {}
        for lineage in lineages:
            logp = math.log(_safe_prob(priors[lineage]))
            ## The sparse model scores both observed and unobserved markers, so
            ## lineages with many direct markers accumulate more negative evidence.
            for allele in selected_alleles:
                target_lineage = max(lineage_freqs[allele], key=lineage_freqs[allele].get)
                p = lineage_freqs[allele][lineage] if lineage == target_lineage else pooled_background[allele]
                p = _safe_prob(p)
                logp += math.log(p) if allele in aset else math.log(1.0 - p)
            log_scores[lineage] = logp
        records[sample_id] = log_scores

    score_df = pd.DataFrame.from_dict(records, orient="index")
    max_scores = score_df.max(axis=1)
    prob_df = np.exp(score_df.sub(max_scores, axis=0))
    prob_df = prob_df.div(prob_df.sum(axis=1), axis=0)
    return score_df, prob_df


def score_samples_sparse(sample_alleles, selected_alleles, lineages, priors, lineage_freqs, pooled_background):
    _score_df, prob_df = score_samples_sparse_tables(
        sample_alleles=sample_alleles,
        selected_alleles=selected_alleles,
        lineages=lineages,
        priors=priors,
        lineage_freqs=lineage_freqs,
        pooled_background=pooled_background,
    )
    return prob_df


def _make_reporting_support_template(evidence_metric="usable_loci"):
    return {
        "status": "not_fit",
        "method": None,
        "evidence_metric": evidence_metric,
        "evidence_bins": None,
        "base_thresholds": None,
        "collision_map": None,
        "thresholds": None,
        "fallback_policy": None,
        "provenance": None,
    }


def _normalise_reporting_support_bins(evidence_bins, evidence_metric):
    if evidence_bins is None:
        evidence_bins = DEFAULT_REPORTING_SUPPORT_EVIDENCE_BINS

    min_key = f"min_{evidence_metric}"
    max_key = f"max_{evidence_metric}"
    normalised = []
    for idx, raw_bin in enumerate(evidence_bins):
        if min_key in raw_bin or max_key in raw_bin:
            min_value = raw_bin.get(min_key)
            max_value = raw_bin.get(max_key)
        else:
            min_value = raw_bin.get("min")
            max_value = raw_bin.get("max")

        bin_id = raw_bin.get("bin_id", f"bin_{idx}")
        normalised.append(
            {
                "bin_id": str(bin_id),
                min_key: int(min_value) if min_value is not None else 0,
                max_key: int(max_value) if max_value is not None else None,
            }
        )

    normalised.sort(key=lambda item: item[min_key])
    return normalised


def _assign_reporting_support_bin(value, evidence_bins, evidence_metric):
    min_key = f"min_{evidence_metric}"
    max_key = f"max_{evidence_metric}"
    value = int(value)

    for evidence_bin in evidence_bins:
        lower = evidence_bin[min_key]
        upper = evidence_bin[max_key]
        if value < lower:
            continue
        if upper is None or value <= upper:
            return evidence_bin["bin_id"]

    return evidence_bins[-1]["bin_id"]


def _top_two_labels(series):
    ranked = series.sort_values(ascending=False)
    top_label = ranked.index[0]
    top_value = float(ranked.iloc[0])

    if len(ranked) > 1:
        runner_up_label = ranked.index[1]
        runner_up_value = float(ranked.iloc[1])
    else:
        runner_up_label = None
        runner_up_value = None

    return top_label, top_value, runner_up_label, runner_up_value


def _collect_loso_prediction_records(
    fold_name,
    score_df,
    prob_df,
    true_labels,
    sample_alleles,
    selected_alleles,
):
    selected_set = set(selected_alleles)
    prediction_records = []

    for sample_id in prob_df.index:
        score_series = score_df.loc[sample_id]
        prob_series = prob_df.loc[sample_id]

        top_class, top_log_score, runner_up, runner_up_log_score = _top_two_labels(score_series)
        _top_prob_label, top_posterior, _runner_prob_label, runner_up_posterior = _top_two_labels(prob_series)

        observed_selected_alleles = len(sample_alleles.get(sample_id, set()) & selected_set)
        pair_metrics = {}
        for winner in score_series.index:
            pair_metrics[winner] = {}
            for competitor in score_series.index:
                if winner == competitor:
                    continue
                pair_metrics[winner][competitor] = {
                    "margin_log": float(score_series[winner] - score_series[competitor]),
                    "posterior": float(prob_series[winner]),
                }

        prediction_records.append(
            {
                "fold": fold_name,
                "sample_id": sample_id,
                "true_class": true_labels[sample_id],
                "top_class": top_class,
                "runner_up": runner_up,
                "top_log_score": float(top_log_score),
                "runner_up_log_score": float(runner_up_log_score) if runner_up_log_score is not None else None,
                "margin_log": float(top_log_score - runner_up_log_score) if runner_up_log_score is not None else None,
                "top_posterior": float(top_posterior),
                "runner_up_posterior": float(runner_up_posterior) if runner_up_posterior is not None else None,
                "observed_selected_alleles": int(observed_selected_alleles),
                "selected_allele_count": int(len(selected_alleles)),
                "pair_metrics": pair_metrics,
            }
        )

    return prediction_records


def _summarise_threshold_records(records, evidence_bin_id, min_samples, quantile):
    if len(records) < min_samples:
        return None

    margin_values = [float(record["margin_log"]) for record in records if record.get("margin_log") is not None]
    posterior_values = [float(record["posterior"]) for record in records if record.get("posterior") is not None]
    if not margin_values or not posterior_values:
        return None

    return {
        "evidence_bin": evidence_bin_id,
        "min_margin_log": float(np.quantile(margin_values, quantile)),
        "min_posterior": float(np.quantile(posterior_values, quantile)),
        "quantile": float(quantile),
        "n_samples": int(len(records)),
    }


def _connected_components(nodes, undirected_edges):
    adjacency = {node: set() for node in nodes}
    for left, right in undirected_edges:
        adjacency.setdefault(left, set()).add(right)
        adjacency.setdefault(right, set()).add(left)

    components = []
    seen = set()
    for node in sorted(adjacency):
        if node in seen:
            continue
        stack = [node]
        component = []
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            component.append(current)
            stack.extend(sorted(adjacency.get(current, []), reverse=True))
        components.append(sorted(component))
    return components


def _resolve_min_lineage_alleles(min_lineage_alleles, min_stable_alleles_per_lineage):
    if min_lineage_alleles is None and min_stable_alleles_per_lineage is None:
        return 0
    if min_lineage_alleles is None:
        return int(min_stable_alleles_per_lineage)
    if min_stable_alleles_per_lineage is None:
        return int(min_lineage_alleles)

    resolved = int(min_lineage_alleles)
    legacy = int(min_stable_alleles_per_lineage)
    if resolved != legacy:
        raise ValueError(
            "min_lineage_alleles and min_stable_alleles_per_lineage were both provided with different values."
        )
    return resolved


def _normalise_parent_map(parent_map):
    if parent_map is None:
        return None

    normalised = {}
    for child, parent in dict(parent_map).items():
        child_id = str(child)
        if parent is None:
            parent_id = None
        else:
            parent_id = str(parent)
            if parent_id.strip() == "":
                parent_id = None
        normalised[child_id] = parent_id

    for parent_id in list(normalised.values()):
        if parent_id is not None and parent_id not in normalised:
            normalised[parent_id] = None

    return normalised


def _build_children_by_parent(parent_map):
    children_by_parent = defaultdict(list)
    for child, parent in parent_map.items():
        if parent is None:
            continue
        children_by_parent[parent].append(child)

    return {
        parent: sorted(children)
        for parent, children in children_by_parent.items()
    }


def _find_roots(parent_map):
    return sorted([node for node, parent in parent_map.items() if parent is None])


def _validate_parent_map(parent_map, leaf_labels):
    normalised_parent_map = _normalise_parent_map(parent_map)
    if normalised_parent_map is None:
        return None

    missing_leaves = sorted(set(leaf_labels) - set(normalised_parent_map))
    if missing_leaves:
        raise ValueError(f"Leaf lineages missing from parent_map: {missing_leaves}")

    for node, parent in normalised_parent_map.items():
        if parent is not None and parent not in normalised_parent_map:
            raise ValueError(f"Parent {parent!r} for node {node!r} is not defined in parent_map.")

    visiting = set()
    visited = set()

    def visit(node):
        ## DFS cycle detection is simpler than trying to reason about malformed
        ## nested structures after the fact.
        if node in visited:
            return
        if node in visiting:
            raise ValueError(f"Cycle detected in parent_map involving node {node!r}.")
        visiting.add(node)
        parent = normalised_parent_map[node]
        if parent is not None:
            visit(parent)
        visiting.remove(node)
        visited.add(node)

    for node in sorted(normalised_parent_map):
        visit(node)

    roots = _find_roots(normalised_parent_map)
    if not roots:
        raise ValueError("parent_map did not produce any roots.")

    return normalised_parent_map


def _get_root_label(leaf, parent_map):
    current = str(leaf)
    while parent_map.get(current) is not None:
        current = parent_map[current]
    return current


def _get_direct_child_label(leaf, node, parent_map):
    current = str(leaf)
    while current in parent_map:
        parent = parent_map.get(current)
        if parent == node:
            return current
        if parent is None:
            break
        current = parent
    raise ValueError(f"Leaf {leaf!r} is not a descendant of node {node!r}.")


def _get_descendant_leaves(node, children_by_parent, leaf_set):
    descendants = []
    stack = list(reversed(children_by_parent.get(node, [])))
    while stack:
        current = stack.pop()
        if current in leaf_set:
            descendants.append(current)
            continue
        stack.extend(reversed(children_by_parent.get(current, [])))
    return sorted(set(descendants))


def _iter_internal_nodes(parent_map, leaf_set):
    children_by_parent = _build_children_by_parent(parent_map)
    internal_nodes = []
    for node in sorted(children_by_parent):
        if len(children_by_parent[node]) > 1:
            internal_nodes.append(node)
    return internal_nodes


def _build_node_training_meta(meta, lineage_col, node, parent_map, children_by_parent, leaf_set):
    if node == HIERARCHY_ROOT_NODE:
        ## The synthetic root classifier exists only when the tree has multiple
        ## disconnected roots; it relabels leaves to their top-level branch.
        descendant_leaves = sorted(leaf_set)
        node_meta = meta[meta[lineage_col].isin(descendant_leaves)].copy()
        node_meta["_node_lineage"] = node_meta[lineage_col].map(lambda leaf: _get_root_label(leaf, parent_map))
        node_children = sorted(node_meta["_node_lineage"].dropna().unique())
        return node_meta, descendant_leaves, node_children

    ## Each local node model competes only among the node's direct children, even
    ## though its training samples are drawn from all descendant leaves below it.
    descendant_leaves = _get_descendant_leaves(node, children_by_parent, leaf_set)
    node_meta = meta[meta[lineage_col].isin(descendant_leaves)].copy()
    node_meta["_node_lineage"] = node_meta[lineage_col].map(
        lambda leaf: _get_direct_child_label(leaf, node, parent_map)
    )
    node_children = sorted(node_meta["_node_lineage"].dropna().unique())
    return node_meta, descendant_leaves, node_children


def fit_reporting_support_from_loso(
    prediction_records,
    lineages,
    evidence_metric="observed_selected_alleles",
    evidence_bins=None,
    threshold_quantile=0.95,
    min_threshold_samples=2,
    collision_min_pair_count=2,
    collision_min_pair_rate=0.05,
):
    support = _make_reporting_support_template(evidence_metric=evidence_metric)
    if not prediction_records:
        return support

    normalised_bins = _normalise_reporting_support_bins(evidence_bins, evidence_metric)
    for record in prediction_records:
        record["evidence_bin"] = _assign_reporting_support_bin(
            value=record[evidence_metric],
            evidence_bins=normalised_bins,
            evidence_metric=evidence_metric,
        )

    ## Reporting support is fit only on held-out LOSO mistakes so the ambiguity
    ## thresholds reflect empirical collision patterns rather than training fit.
    incorrect_predictions = [record for record in prediction_records if record["top_class"] != record["true_class"]]
    base_by_bin = []
    for evidence_bin in normalised_bins:
        bin_records = [
            {
                "margin_log": record["margin_log"],
                "posterior": record["top_posterior"],
            }
            for record in incorrect_predictions
            if record["evidence_bin"] == evidence_bin["bin_id"]
        ]
        threshold_summary = _summarise_threshold_records(
            records=bin_records,
            evidence_bin_id=evidence_bin["bin_id"],
            min_samples=min_threshold_samples,
            quantile=threshold_quantile,
        )
        if threshold_summary is not None:
            base_by_bin.append(threshold_summary)

    base_overall = _summarise_threshold_records(
        records=[
            {
                "margin_log": record["margin_log"],
                "posterior": record["top_posterior"],
            }
            for record in incorrect_predictions
        ],
        evidence_bin_id="overall",
        min_samples=min_threshold_samples,
        quantile=threshold_quantile,
    )

    true_class_counts = Counter(record["true_class"] for record in prediction_records)
    ordered_confusions = Counter()
    unordered_confusions = Counter()
    for record in incorrect_predictions:
        winner = record["top_class"]
        truth = record["true_class"]
        ordered_confusions[(winner, truth)] += 1
        unordered_confusions[tuple(sorted((winner, truth)))] += 1

    collision_pairs = []
    recurrent_edges = []
    for (left, right), confusion_count in sorted(unordered_confusions.items()):
        denominator = true_class_counts[left] + true_class_counts[right]
        confusion_rate = float(confusion_count / denominator) if denominator else 0.0
        ordered_counts = {
            f"{left}>{right}": int(ordered_confusions.get((left, right), 0)),
            f"{right}>{left}": int(ordered_confusions.get((right, left), 0)),
        }
        is_recurrent = confusion_count >= collision_min_pair_count and confusion_rate >= collision_min_pair_rate
        collision_pairs.append(
            {
                "members": [left, right],
                "confusion_count": int(confusion_count),
                "confusion_rate": confusion_rate,
                "ordered_confusion_counts": ordered_counts,
                "is_recurrent": bool(is_recurrent),
            }
        )
        if is_recurrent:
            recurrent_edges.append((left, right))

    collision_groups = []
    group_lookup = defaultdict(list)
    for group_index, members in enumerate(_connected_components(lineages, recurrent_edges)):
        if len(members) < 2:
            continue
        group_id = f"collision_group_{group_index}"
        collision_groups.append({"group_id": group_id, "members": members})
        for member in members:
            group_lookup[member].append(group_id)

    pairwise_thresholds = []
    for winner in sorted(lineages):
        for competitor in sorted(lineages):
            if winner == competitor:
                continue

            pair_error_records = []
            for record in prediction_records:
                if record["true_class"] != competitor:
                    continue
                pair_metric = record["pair_metrics"].get(winner, {}).get(competitor)
                if pair_metric is None:
                    continue
                if pair_metric["margin_log"] <= 0.0:
                    continue
                pair_error_records.append(
                    {
                        "evidence_bin": record["evidence_bin"],
                        "margin_log": pair_metric["margin_log"],
                        "posterior": pair_metric["posterior"],
                    }
                )

            overall_threshold = _summarise_threshold_records(
                records=pair_error_records,
                evidence_bin_id="overall",
                min_samples=min_threshold_samples,
                quantile=threshold_quantile,
            )
            if overall_threshold is None:
                continue

            thresholds_by_bin = []
            for evidence_bin in normalised_bins:
                bin_threshold = _summarise_threshold_records(
                    records=[
                        record for record in pair_error_records
                        if record["evidence_bin"] == evidence_bin["bin_id"]
                    ],
                    evidence_bin_id=evidence_bin["bin_id"],
                    min_samples=min_threshold_samples,
                    quantile=threshold_quantile,
                )
                if bin_threshold is not None:
                    thresholds_by_bin.append(bin_threshold)

            pairwise_thresholds.append(
                {
                    "winner": winner,
                    "competitor": competitor,
                    "overall": overall_threshold,
                    "by_evidence_bin": thresholds_by_bin,
                    "group_ids": group_lookup.get(winner, []),
                }
            )

    support.update(
        {
            "status": "fit",
            "method": "loso_confusion_pairwise_thresholds",
            "evidence_metric": evidence_metric,
            "evidence_bins": normalised_bins,
            "base_thresholds": {
                "overall": base_overall,
                "by_evidence_bin": base_by_bin,
            },
            "collision_map": {
                "pairs": collision_pairs,
                "groups": collision_groups,
            },
            "thresholds": {
                "pairwise": pairwise_thresholds,
            },
            "fallback_policy": {
                "reporting_status": "ambiguous",
                "return_value": "ambiguity_set",
            },
            "provenance": {
                "fit_on": "loso_held_out_predictions",
                "threshold_quantile": float(threshold_quantile),
                "min_threshold_samples": int(min_threshold_samples),
                "collision_min_pair_count": int(collision_min_pair_count),
                "collision_min_pair_rate": float(collision_min_pair_rate),
                "prediction_record_count": int(len(prediction_records)),
                "incorrect_prediction_count": int(len(incorrect_predictions)),
            },
        }
    )
    return support


def build_sparse_empirical_model_from_sets(
    sample_alleles,
    labels_by_sample,
    selected_alleles,
    species_id="Plasmodium_falciparum",
    model_id="pf_geolineage_empirical_model",
    smoothing_alpha=1.0,
    uniform_priors=False,
):
    lineages, lineage_to_samples, lineage_freqs, pooled_background = estimate_lineage_frequencies(
        sample_alleles=sample_alleles,
        labels_by_sample=labels_by_sample,
        selected_alleles=selected_alleles,
        smoothing_alpha=smoothing_alpha,
    )

    lineage_counts = {lineage: len(samples) for lineage, samples in lineage_to_samples.items()}
    total_samples = sum(lineage_counts.values())
    
    if uniform_priors:
        priors = {lineage: 1/len(lineages) for lineage in lineages}
    else:
        priors = {lineage: lineage_counts[lineage] / total_samples for lineage in lineages}

    loci = []
    direct_locus_ids_by_lineage = defaultdict(list)

    for allele in selected_alleles:
        allele_meta = parse_ardal_allele(allele)
        target_lineage = max(lineage_freqs[allele], key=lineage_freqs[allele].get)

        ref_token = allele_meta["ref"] if allele_meta["ref"] is not None else "N"
        locus_id = (
            f"{allele_meta['chrom']}:{allele_meta['pos']}:"
            f"{ref_token}>{allele_meta['alt']}:{target_lineage}"
        )

        loci.append(
            {
                "locus_id": locus_id,
                "chrom": allele_meta["chrom"],
                "pos": allele_meta["pos"],
                "bed_start": allele_meta["bed_start"],
                "bed_end": allele_meta["bed_end"],
                "ref": allele_meta["ref"],
                "alt": allele_meta["alt"],
                "target_lineage": target_lineage,
                "empirical_bayes": {
                    "target_frequency": lineage_freqs[allele][target_lineage],
                    "background_frequency": pooled_background[allele],
                },
                "annotations": {
                    "source_allele_id": str(allele),
                    "lineage_frequencies": lineage_freqs[allele],
                    "smoothing_alpha": float(smoothing_alpha),
                },
            }
        )
        direct_locus_ids_by_lineage[target_lineage].append(locus_id)

    lineage_records = []
    for lineage in lineages:
        direct_locus_ids = sorted(direct_locus_ids_by_lineage.get(lineage, []))
        lineage_records.append(
            {
                "lineage_id": lineage,
                "prior": float(priors[lineage]),
                "parent_lineage": None,
                "is_root": True,
                "ancestor_lineages": [],
                "direct_children": [],
                "direct_locus_ids": direct_locus_ids,
                "inherited_locus_ids": direct_locus_ids,
                "direct_marker_count": len(direct_locus_ids),
                "inherited_marker_count": len(direct_locus_ids),
            }
        )

    created_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    model = {
        "schema_version": "1.1",
        "model_id": model_id,
        "species_id": species_id,
        "reference": {
            "name": None,
            "path": None,
            "contigs": sorted({l["chrom"] for l in loci}),
            "ref_alleles_resolved": any(l["ref"] is not None for l in loci),
        },
        "model_type": "empirical_geolineage",
        "hierarchy": {
            "is_hierarchical": False,
            "root_lineages": [x["lineage_id"] for x in lineage_records],
            "edges": [],
            "lineage_naming_rule": None,
            "note": "Geolineages are treated as a flat multiclass set.",
        },
        "lineages": lineage_records,
        "loci": loci,
        "empirical_bayes": {
            "emission_family": "bernoulli",
            "smoothing_alpha": float(smoothing_alpha),
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
        "reporting_support": _make_reporting_support_template(
            evidence_metric="observed_selected_alleles"
        ),
        "summary": {
            "lineage_count": len(lineage_records),
            "locus_count": len(loci),
            "has_resolved_reference_alleles": any(l["ref"] is not None for l in loci),
        },
        "provenance": {
            "created_at": created_at,
            "created_by": "build_geolineage_min_model_sparse",
            "input_format": "ardal_multiclass_loso_sparse",
        },
    }

    return model


def _build_candidate_pool(
    ard_train,
    train_meta,
    sample_col,
    lineage_col,
    method,
    top_k_per_lineage,
    max_model_alleles,
    threads,
    cooc_threshold,
):
    lineage_rankings = {}
    candidate_meta = {}

    for lineage in sorted(train_meta[lineage_col].unique()):
        positives = train_meta.loc[train_meta[lineage_col] == lineage, sample_col].tolist()
        ranked_items = list(ard_train.stats.allele_inform(positives, method=method).items())[:top_k_per_lineage]
        lineage_rankings[lineage] = ranked_items

        for rank, (allele, score) in enumerate(ranked_items, start=1):
            meta = candidate_meta.setdefault(
                allele,
                {
                    "best_lineage": lineage,
                    "best_score": float(score),
                    "best_rank": rank,
                    "lineage_scores": {},
                },
            )
            meta["lineage_scores"][lineage] = float(score)
            if (score > meta["best_score"]) or (math.isclose(score, meta["best_score"]) and rank < meta["best_rank"]):
                meta["best_lineage"] = lineage
                meta["best_score"] = float(score)
                meta["best_rank"] = rank

    ## Interleave per-lineage rankings so one lineage cannot monopolise the
    ## candidate pool before LOSO has a chance to evaluate balance explicitly.
    ordered_candidates = []
    seen = set()
    max_rank_depth = max((len(ranked) for ranked in lineage_rankings.values()), default=0)
    for rank_idx in range(max_rank_depth):
        for lineage in sorted(lineage_rankings):
            ranked_items = lineage_rankings[lineage]
            if rank_idx >= len(ranked_items):
                continue
            allele = ranked_items[rank_idx][0]
            if allele in seen:
                continue
            seen.add(allele)
            ordered_candidates.append(allele)

    print(f"candidate_union={len(ordered_candidates)}")

    if cooc_threshold is not None and ordered_candidates:
        cooc_map = ard_train.stats.allele_cooc(
            ordered_candidates,
            threshold=cooc_threshold,
            threads=threads,
        )
        cooc_map = {allele: set(partners) for allele, partners in cooc_map.items()}

        pruned_candidates = []
        for allele in ordered_candidates:
            allele_lineage = candidate_meta[allele]["best_lineage"]
            redundant = False
            for kept in pruned_candidates:
                if candidate_meta[kept]["best_lineage"] != allele_lineage:
                    continue
                kept_partners = cooc_map.get(kept, set())
                allele_partners = cooc_map.get(allele, set())
                if allele in kept_partners or kept in allele_partners:
                    redundant = True
                    break
            if not redundant:
                pruned_candidates.append(allele)

        ordered_candidates = pruned_candidates
        print(f"cooc_pruned_candidates={len(ordered_candidates)}")

    ordered_candidates = ordered_candidates[:max_model_alleles]
    print(f"ordered_candidates={len(ordered_candidates)}")

    return ordered_candidates, candidate_meta


def _select_fold_winner(prefix_results, score_tolerance):
    best_score = max(result["score"] for result in prefix_results)
    competitive_results = [
        result for result in prefix_results
        if result["score"] >= best_score - score_tolerance
    ]

    return max(
        competitive_results,
        key=lambda result: (
            result["represented_lineages"],
            result["min_lineage_markers"],
            result["num_alleles"],
            result["score"],
        ),
    )


def _consensus_target_lineages(allele_target_support):
    consensus = {}
    for allele, lineage_counts in allele_target_support.items():
        if lineage_counts:
            consensus[allele] = lineage_counts.most_common(1)[0][0]
    return consensus


def _stable_allele_order(alleles, allele_support, allele_target_support):
    consensus = _consensus_target_lineages(allele_target_support)
    return sorted(
        alleles,
        key=lambda allele: (
            -allele_support[allele],
            -max(allele_target_support[allele].values(), default=0),
            consensus.get(allele, ""),
            str(allele),
        ),
    )


def _supplement_lineage_markers(
    stable_alleles,
    lineages,
    allele_support,
    allele_target_support,
    min_lineage_alleles,
):
    if min_lineage_alleles <= 0:
        return stable_alleles

    consensus = _consensus_target_lineages(allele_target_support)
    stable = list(stable_alleles)
    stable_set = set(stable)
    lineage_counts = Counter(consensus.get(allele) for allele in stable if consensus.get(allele) is not None)

    ## This is a floor, not a cap: supported markers are added for underpowered
    ## classes, but richer lineages are allowed to keep more than the minimum.
    for lineage in sorted(lineages):
        needed = max(0, min_lineage_alleles - lineage_counts.get(lineage, 0))
        if needed == 0:
            continue

        lineage_candidates = [
            allele for allele, counts in allele_target_support.items()
            if counts.get(lineage, 0) > 0 and allele not in stable_set
        ]
        lineage_candidates.sort(
            key=lambda allele: (
                -allele_target_support[allele][lineage],
                -allele_support[allele],
                str(allele),
            ),
        )

        for allele in lineage_candidates[:needed]:
            stable.append(allele)
            stable_set.add(allele)
            lineage_counts[lineage] += 1

    return stable


def _model_id_for_node(model_id, node_name):
    if node_name is None:
        return model_id
    safe_node_name = re.sub(r"[^0-9A-Za-z_.-]+", "_", str(node_name))
    return f"{model_id}__{safe_node_name}"


def _build_flat_sparse_node_model(
    ard_obj,
    meta_df,
    sample_col="Sample",
    lineage_col="GeoLineage",
    study_col="Country",
    method="kullbackleibler",
    top_k_per_lineage=200,
    max_model_alleles=500,
    min_lineage_samples=10,
    min_fold_support=0.5,
    smoothing_alpha=1.0,
    uniform_priors=False,
    threads=16,
    metric="macro_f1",
    species_id="Plasmodium_falciparum",
    model_id="pf_geolineage_empirical_model",
    cooc_threshold=None,
    score_tolerance=0.005,
    min_lineage_alleles=0,
    reporting_support_bins=None,
    reporting_support_quantile=0.95,
    reporting_support_min_samples=2,
    collision_min_pair_count=2,
    collision_min_pair_rate=0.05,
    node_name=None,
):
    meta = meta_df[[sample_col, lineage_col, study_col]].dropna().copy()
    lineage_counts = meta[lineage_col].value_counts()
    keep_lineages = set(lineage_counts[lineage_counts >= min_lineage_samples].index)
    meta = meta[meta[lineage_col].isin(keep_lineages)].copy()
    if meta.empty or meta[lineage_col].nunique() < 2:
        return None, []

    ## LOSO folds are only valid when every held-out label is represented in the
    ## corresponding training split; otherwise that fold cannot score all classes.
    folds = []
    for study in sorted(meta[study_col].unique()):
        test_ids = meta.loc[meta[study_col] == study, sample_col].tolist()
        train_ids = meta.loc[meta[study_col] != study, sample_col].tolist()

        train_lineages = set(meta.loc[meta[sample_col].isin(train_ids), lineage_col])
        test_lineages = set(meta.loc[meta[sample_col].isin(test_ids), lineage_col])

        if test_lineages.issubset(train_lineages):
            folds.append(LosoFold(study, test_ids, train_ids))

    if not folds:
        return None, []

    allele_support = Counter()
    allele_target_support = defaultdict(Counter)
    fold_results = []
    loso_prediction_records = []
    local_model_id = _model_id_for_node(model_id, node_name)

    for fold in folds:
        print(f"Starting {fold.study} training with {len(fold.train)} samples...")

        train_meta = meta[meta[sample_col].isin(fold.train)].copy()
        test_meta = meta[meta[sample_col].isin(fold.test)].copy()

        ard_train = ard_obj.get.subset(
            guids=train_meta[sample_col].tolist(),
            chunk_size=1000,
            threads=12,
            drop_zero_cols=False,
        )
        ard_test = ard_obj.get.subset(
            guids=test_meta[sample_col].tolist(),
            chunk_size=1000,
            threads=12,
            drop_zero_cols=False,
        )

        ordered_candidates, _candidate_meta = _build_candidate_pool(
            ard_train=ard_train,
            train_meta=train_meta,
            sample_col=sample_col,
            lineage_col=lineage_col,
            method=method,
            top_k_per_lineage=top_k_per_lineage,
            max_model_alleles=max_model_alleles,
            threads=threads,
            cooc_threshold=cooc_threshold,
        )

        train_sample_alleles = build_sample_allele_sets(
            ard_subset=ard_train,
            alleles=ordered_candidates,
            sample_ids=train_meta[sample_col].tolist(),
        )
        print(f"train_sample_alleles={len(train_sample_alleles)}")

        test_sample_alleles = build_sample_allele_sets(
            ard_subset=ard_test,
            alleles=ordered_candidates,
            sample_ids=test_meta[sample_col].tolist(),
        )
        print(f"test_sample_alleles={len(test_sample_alleles)}")

        y_train = dict(zip(train_meta[sample_col], train_meta[lineage_col]))
        y_test = pd.Series(
            [dict(zip(test_meta[sample_col], test_meta[lineage_col]))[sid] for sid in test_meta[sample_col]],
            index=test_meta[sample_col],
        )

        selected = []
        prefix_results = []

        ## Evaluate prefix panels in order so LOSO chooses the smallest strong
        ## local model rather than taking the whole candidate pool by default.
        for allele in ordered_candidates:
            selected.append(allele)

            model = build_sparse_empirical_model_from_sets(
                sample_alleles=train_sample_alleles,
                labels_by_sample=y_train,
                selected_alleles=selected,
                species_id=species_id,
                model_id=local_model_id,
                smoothing_alpha=smoothing_alpha,
                uniform_priors=uniform_priors,
            )

            priors = {x["lineage_id"]: x["prior"] for x in model["lineages"]}
            lineages = [x["lineage_id"] for x in model["lineages"]]

            lineage_freqs = {
                locus["annotations"]["source_allele_id"]: locus["annotations"]["lineage_frequencies"]
                for locus in model["loci"]
            }
            pooled_background = {
                locus["annotations"]["source_allele_id"]: locus["empirical_bayes"]["background_frequency"]
                for locus in model["loci"]
            }

            prob_df = score_samples_sparse(
                sample_alleles=test_sample_alleles,
                selected_alleles=selected,
                lineages=lineages,
                priors=priors,
                lineage_freqs=lineage_freqs,
                pooled_background=pooled_background,
            )

            y_pred = prob_df.idxmax(axis=1)
            score = _metric_score(y_test.loc[prob_df.index], y_pred, prob_df, metric)
            lineage_marker_counts = {
                lineage_record["lineage_id"]: lineage_record["direct_marker_count"]
                for lineage_record in model["lineages"]
            }
            allele_target_lineages = {
                locus["annotations"]["source_allele_id"]: locus["target_lineage"]
                for locus in model["loci"]
            }

            prefix_results.append(
                {
                    "study": fold.study,
                    "num_alleles": len(selected),
                    "metric": metric,
                    "score": float(score),
                    "selected_alleles": selected.copy(),
                    "lineage_marker_counts": lineage_marker_counts,
                    "represented_lineages": sum(
                        1 for count in lineage_marker_counts.values() if count > 0
                    ),
                    "min_lineage_markers": min(lineage_marker_counts.values()) if lineage_marker_counts else 0,
                    "allele_target_lineages": allele_target_lineages,
                }
            )

        if not prefix_results:
            continue

        best_result = _select_fold_winner(prefix_results, score_tolerance=score_tolerance)
        best_result["best_fold_score"] = max(result["score"] for result in prefix_results)
        best_result["score_tolerance"] = float(score_tolerance)

        allele_support.update(best_result["selected_alleles"])
        for allele, target_lineage in best_result["allele_target_lineages"].items():
            if allele in best_result["selected_alleles"]:
                allele_target_support[allele][target_lineage] += 1

        best_model = build_sparse_empirical_model_from_sets(
            sample_alleles=train_sample_alleles,
            labels_by_sample=y_train,
            selected_alleles=best_result["selected_alleles"],
            species_id=species_id,
            model_id=local_model_id,
            smoothing_alpha=smoothing_alpha,
            uniform_priors=uniform_priors,
        )
        best_priors = {x["lineage_id"]: x["prior"] for x in best_model["lineages"]}
        best_lineages = [x["lineage_id"] for x in best_model["lineages"]]
        best_lineage_freqs = {
            locus["annotations"]["source_allele_id"]: locus["annotations"]["lineage_frequencies"]
            for locus in best_model["loci"]
        }
        best_pooled_background = {
            locus["annotations"]["source_allele_id"]: locus["empirical_bayes"]["background_frequency"]
            for locus in best_model["loci"]
        }
        best_score_df, best_prob_df = score_samples_sparse_tables(
            sample_alleles=test_sample_alleles,
            selected_alleles=best_result["selected_alleles"],
            lineages=best_lineages,
            priors=best_priors,
            lineage_freqs=best_lineage_freqs,
            pooled_background=best_pooled_background,
        )
        y_test_dict = dict(zip(test_meta[sample_col], test_meta[lineage_col]))
        loso_prediction_records.extend(
            _collect_loso_prediction_records(
                fold_name=fold.study,
                score_df=best_score_df,
                prob_df=best_prob_df,
                true_labels=y_test_dict,
                sample_alleles=test_sample_alleles,
                selected_alleles=best_result["selected_alleles"],
            )
        )

        fold_results.append(best_result)

    if not fold_results:
        return None, []

    ## Stable alleles are the fold-winning markers that recur across enough LOSO
    ## splits, then optionally topped up to satisfy the per-lineage minimum.
    min_support_count = max(1, math.ceil(min_fold_support * len(fold_results)))
    stable_alleles = [
        allele for allele, count in allele_support.items()
        if count >= min_support_count
    ]
    stable_alleles = _stable_allele_order(stable_alleles, allele_support, allele_target_support)
    stable_alleles = _supplement_lineage_markers(
        stable_alleles=stable_alleles,
        lineages=sorted(meta[lineage_col].unique()),
        allele_support=allele_support,
        allele_target_support=allele_target_support,
        min_lineage_alleles=min_lineage_alleles,
    )

    print(f"stable_alleles={len(stable_alleles)}")
    if not stable_alleles:
        best_idx = int(np.argmax([x["score"] for x in fold_results]))
        stable_alleles = fold_results[best_idx]["selected_alleles"]

    full_ids = meta[sample_col].tolist()
    ard_full = ard_obj.get.subset(guids=full_ids, chunk_size=1000, threads=12, drop_zero_cols=False)
    full_sample_alleles = build_sample_allele_sets(
        ard_subset=ard_full,
        alleles=stable_alleles,
        sample_ids=full_ids,
    )
    print(f"full_sample_alleles={len(full_sample_alleles)}")
    full_labels = dict(zip(meta[sample_col], meta[lineage_col]))

    final_model = build_sparse_empirical_model_from_sets(
        sample_alleles=full_sample_alleles,
        labels_by_sample=full_labels,
        selected_alleles=stable_alleles,
        species_id=species_id,
        model_id=local_model_id,
        smoothing_alpha=smoothing_alpha,
        uniform_priors=uniform_priors,
    )

    for locus in final_model["loci"]:
        allele = locus["annotations"]["source_allele_id"]
        locus["annotations"]["loso_support"] = {
            "fold_count": int(allele_support.get(allele, 0)),
            "target_lineage_fold_counts": dict(allele_target_support.get(allele, {})),
        }

    final_model["provenance"]["loso"] = {
        "study_column": study_col,
        "lineage_column": lineage_col,
        "sample_column": sample_col,
        "num_folds": len(fold_results),
        "selection_metric": metric,
        "min_fold_support": float(min_fold_support),
        "cooc_threshold": None if cooc_threshold is None else float(cooc_threshold),
        "score_tolerance": float(score_tolerance),
        "min_lineage_alleles": int(min_lineage_alleles),
        "reporting_support_quantile": float(reporting_support_quantile),
        "reporting_support_min_samples": int(reporting_support_min_samples),
        "collision_min_pair_count": int(collision_min_pair_count),
        "collision_min_pair_rate": float(collision_min_pair_rate),
        "candidate_strategy": "round_robin_top_k_per_lineage",
    }
    final_model["provenance"]["fold_results"] = fold_results
    final_model["provenance"]["loso_prediction_records"] = loso_prediction_records
    final_model["summary"]["selected_locus_count_after_loso"] = len(stable_alleles)
    final_model["reporting_support"] = fit_reporting_support_from_loso(
        prediction_records=loso_prediction_records,
        lineages=[record["lineage_id"] for record in final_model["lineages"]],
        evidence_metric="observed_selected_alleles",
        evidence_bins=reporting_support_bins,
        threshold_quantile=reporting_support_quantile,
        min_threshold_samples=reporting_support_min_samples,
        collision_min_pair_count=collision_min_pair_count,
        collision_min_pair_rate=collision_min_pair_rate,
    )

    return final_model, fold_results


def build_geolineage_min_model_sparse(
    ard_obj,
    meta_df,
    parent_map=None,
    sample_col="Sample",
    lineage_col="Lineage",
    study_col="Country",
    method="kullbackleibler",
    top_k_per_lineage=200,
    max_model_alleles=500,
    min_lineage_samples=10,
    min_fold_support=0.5,
    smoothing_alpha=1.0,
    uniform_priors=False,
    threads=16,
    metric="macro_f1",
    species_id="Plasmodium_falciparum",
    model_id="pf_geolineage_empirical_model",
    cooc_threshold=None,
    score_tolerance=0.005,
    min_lineage_alleles=None,
    min_stable_alleles_per_lineage=None,
    reporting_support_bins=None,
    reporting_support_quantile=0.95,
    reporting_support_min_samples=2,
    collision_min_pair_count=2,
    collision_min_pair_rate=0.05,
):
    resolved_min_lineage_alleles = _resolve_min_lineage_alleles(
        min_lineage_alleles=min_lineage_alleles,
        min_stable_alleles_per_lineage=min_stable_alleles_per_lineage,
    )

    if parent_map is None:
        ## Preserve the original flat behavior unless an explicit hierarchy is
        ## provided.
        final_model, fold_results = _build_flat_sparse_node_model(
            ard_obj=ard_obj,
            meta_df=meta_df,
            sample_col=sample_col,
            lineage_col=lineage_col,
            study_col=study_col,
            method=method,
            top_k_per_lineage=top_k_per_lineage,
            max_model_alleles=max_model_alleles,
            min_lineage_samples=min_lineage_samples,
            min_fold_support=min_fold_support,
            smoothing_alpha=smoothing_alpha,
            uniform_priors=uniform_priors,
            threads=threads,
            metric=metric,
            species_id=species_id,
            model_id=model_id,
            cooc_threshold=cooc_threshold,
            score_tolerance=score_tolerance,
            min_lineage_alleles=resolved_min_lineage_alleles,
            reporting_support_bins=reporting_support_bins,
            reporting_support_quantile=reporting_support_quantile,
            reporting_support_min_samples=reporting_support_min_samples,
            collision_min_pair_count=collision_min_pair_count,
            collision_min_pair_rate=collision_min_pair_rate,
        )
        if final_model is None:
            raise ValueError("No valid LOSO folds were produced.")
        return final_model, fold_results

    meta = meta_df[[sample_col, lineage_col, study_col]].dropna().copy()
    lineage_counts = meta[lineage_col].value_counts()
    keep_lineages = set(lineage_counts[lineage_counts >= min_lineage_samples].index)
    meta = meta[meta[lineage_col].isin(keep_lineages)].copy()
    if meta.empty:
        raise ValueError("No samples remain after applying min_lineage_samples.")

    leaf_set = set(meta[lineage_col].unique())
    normalised_parent_map = _validate_parent_map(parent_map, leaf_set)
    children_by_parent = _build_children_by_parent(normalised_parent_map)
    roots = _find_roots(normalised_parent_map)

    ## Train one local classifier per internal node; each node only competes
    ## among its direct children.
    nodes_to_train = []
    if len(roots) > 1:
        nodes_to_train.append(HIERARCHY_ROOT_NODE)
    for node in _iter_internal_nodes(normalised_parent_map, leaf_set):
        if node not in nodes_to_train:
            nodes_to_train.append(node)

    node_models = {}
    fold_results_by_node = {}

    for node in nodes_to_train:
        node_meta, descendant_leaves, node_children = _build_node_training_meta(
            meta=meta,
            lineage_col=lineage_col,
            node=node,
            parent_map=normalised_parent_map,
            children_by_parent=children_by_parent,
            leaf_set=leaf_set,
        )
        if len(node_children) < 2:
            continue

        local_model, local_fold_results = _build_flat_sparse_node_model(
            ard_obj=ard_obj,
            meta_df=node_meta,
            sample_col=sample_col,
            lineage_col="_node_lineage",
            study_col=study_col,
            method=method,
            top_k_per_lineage=top_k_per_lineage,
            max_model_alleles=max_model_alleles,
            min_lineage_samples=min_lineage_samples,
            min_fold_support=min_fold_support,
            smoothing_alpha=smoothing_alpha,
            uniform_priors=uniform_priors,
            threads=threads,
            metric=metric,
            species_id=species_id,
            model_id=model_id,
            cooc_threshold=cooc_threshold,
            score_tolerance=score_tolerance,
            min_lineage_alleles=resolved_min_lineage_alleles,
            reporting_support_bins=reporting_support_bins,
            reporting_support_quantile=reporting_support_quantile,
            reporting_support_min_samples=reporting_support_min_samples,
            collision_min_pair_count=collision_min_pair_count,
            collision_min_pair_rate=collision_min_pair_rate,
            node_name=node,
        )
        if local_model is None:
            continue

        local_model["provenance"]["node"] = {
            "node_name": node,
            "node_children": list(node_children),
            "node_descendant_leaves": list(descendant_leaves),
            "min_lineage_alleles": int(resolved_min_lineage_alleles),
        }
        node_models[node] = local_model
        fold_results_by_node[node] = local_fold_results

    if not node_models:
        raise ValueError("No hierarchical node produced a valid model.")

    serialised_parent_map = dict(normalised_parent_map)
    serialised_children = dict(children_by_parent)
    serialised_roots = list(roots)
    top_classifier_node = None
    if len(roots) > 1:
        ## Multiple roots need a synthetic top-level classifier so inference has
        ## an explicit first decision point.
        serialised_parent_map[HIERARCHY_ROOT_NODE] = None
        serialised_children[HIERARCHY_ROOT_NODE] = list(roots)
        for root in roots:
            serialised_parent_map[root] = HIERARCHY_ROOT_NODE
        serialised_roots = [HIERARCHY_ROOT_NODE]
        top_classifier_node = HIERARCHY_ROOT_NODE
    elif len(roots) == 1 and roots[0] in node_models:
        top_classifier_node = roots[0]

    created_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    hierarchical_model = {
        "schema_version": "1.1",
        "model_id": model_id,
        "species_id": species_id,
        "model_type": "hierarchical_empirical_geolineage",
        "hierarchy": {
            "is_hierarchical": True,
            "parent_map": serialised_parent_map,
            "root_lineages": serialised_roots,
            "children_by_parent": serialised_children,
            "leaf_lineages": sorted(leaf_set),
            "top_classifier_node": top_classifier_node,
            "lineage_naming_rule": None,
            "note": "Lineage hierarchy is explicit. Each node model classifies among direct children only.",
        },
        "node_models": node_models,
        "summary": {
            "node_model_count": len(node_models),
            "leaf_lineage_count": len(leaf_set),
            "has_reporting_support": any(
                model.get("reporting_support", {}).get("status") == "fit"
                for model in node_models.values()
            ),
        },
        "provenance": {
            "created_at": created_at,
            "created_by": "build_geolineage_min_model_sparse",
            "input_format": "ardal_multiclass_loso_sparse",
            "sample_column": sample_col,
            "lineage_column": lineage_col,
            "study_column": study_col,
            "min_lineage_alleles": int(resolved_min_lineage_alleles),
            "min_lineage_samples": int(min_lineage_samples),
            "selection_metric": metric,
            "candidate_strategy": "round_robin_top_k_per_lineage",
        },
    }

    return hierarchical_model, fold_results_by_node


def write_model_json(model, output_json):
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as handle:
        json.dump(model, handle, indent=2)
    return output_json

import json
import math
import re
from collections import Counter, defaultdict, namedtuple
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score, log_loss
from scipy import optimize

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
DEFAULT_FULL_BAYES_PRIORS = {
    "target_mu": {"alpha": 7.0, "beta": 5.0},
    "background_mu": {"alpha": 2.0, "beta": 8.0},
    "kappa": {"alpha": 2.0, "beta": 1.0},
}


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

    ## prefer Ardal backed decoding for dotted IDs so contig / position parsing
    ## stays consistent with the allele headers used elsewhere in the stack
    ## TODO: Ardal has shifted its intermediate files schema in to v3 so might need attention at some point
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


def score_samples_full_bayes_tables(sample_alleles, model):
    priors = {x["lineage_id"]: x["prior"] for x in model["lineages"]}
    lineages = [x["lineage_id"] for x in model["lineages"]]
    loci = list(model.get("loci", []))

    records = {}
    for sample_id, aset in sample_alleles.items():
        log_scores = {}
        for lineage in lineages:
            logp = math.log(_safe_prob(priors[lineage]))
            for locus in loci:
                side = "target" if lineage == locus["target_lineage"] else "background"
                p = _safe_prob((locus.get("full_bayes") or {})[f"{side}_mu"])
                allele = locus["annotations"]["source_allele_id"]
                logp += math.log(p) if allele in aset else math.log(1.0 - p)
            log_scores[lineage] = logp
        records[sample_id] = log_scores

    score_df = pd.DataFrame.from_dict(records, orient="index")
    max_scores = score_df.max(axis=1)
    prob_df = np.exp(score_df.sub(max_scores, axis=0))
    prob_df = prob_df.div(prob_df.sum(axis=1), axis=0)
    return score_df, prob_df


def _logit(p):
    p = _safe_prob(p)
    return math.log(p) - math.log(1.0 - p)


def _expit(x):
    if x >= 0.0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


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
        child_ardal_kwargs={"roaring" : True}
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

        ## estimate one smoothed Bernoulli rate per lineage
        ## collapse non target lineages into a shared background emission
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
            ## scores both observed and unobserved
            ## works better
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


def _resolve_full_bayes_priors(full_bayes_priors):
    resolved = json.loads(json.dumps(DEFAULT_FULL_BAYES_PRIORS))
    if full_bayes_priors is None:
        return resolved

    for key, default_value in resolved.items():
        override = full_bayes_priors.get(key)
        if override is None:
            continue
        resolved[key]["alpha"] = float(override["alpha"])
        resolved[key]["beta"] = float(override["beta"])
    return resolved


def _log_beta_binomial_pmf(k, n, alpha, beta):
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


def _smoothed_frequency_from_group_counts(group_counts, smoothing_alpha):
    n_total = sum(int(group["n"]) for group in group_counts)
    k_total = sum(int(group["k"]) for group in group_counts)
    return float((k_total + smoothing_alpha) / (n_total + 2.0 * smoothing_alpha))


def _negative_log_grouped_beta_binomial_posterior(params, group_counts, mu_prior, kappa_prior):
    eta, log_kappa = params
    mu = _expit(float(eta))
    kappa = math.exp(float(log_kappa))
    alpha = mu * kappa
    beta = (1.0 - mu) * kappa

    log_posterior = 0.0
    for group in group_counts:
        n = int(group["n"])
        k = int(group["k"])
        log_posterior += _log_beta_binomial_pmf(k, n, alpha, beta)

    log_posterior += (
        (float(mu_prior["alpha"]) - 1.0) * math.log(mu)
        + (float(mu_prior["beta"]) - 1.0) * math.log(1.0 - mu)
    )
    log_posterior += (
        (float(kappa_prior["alpha"]) - 1.0) * math.log(kappa)
        - float(kappa_prior["beta"]) * kappa
    )

    return -float(log_posterior)


def _fit_beta_binomial_map(group_counts, mu_prior, kappa_prior, default_kappa):
    empirical_mu = _smoothed_frequency_from_group_counts(group_counts, smoothing_alpha=1.0)
    initial = np.array([_logit(empirical_mu), math.log(float(default_kappa))], dtype=float)

    optimisation = optimize.minimize(
        _negative_log_grouped_beta_binomial_posterior,
        x0=initial,
        args=(group_counts, mu_prior, kappa_prior),
        method="L-BFGS-B",
        bounds=[(-12.0, 12.0), (math.log(1e-3), math.log(1e4))],
    )
    if not optimisation.success:
        raise RuntimeError(str(optimisation.message))

    mu = _expit(float(optimisation.x[0]))
    kappa = math.exp(float(optimisation.x[1]))
    return {
        "mu": float(mu),
        "kappa": float(kappa),
        "alpha": float(mu * kappa),
        "beta": float((1.0 - mu) * kappa),
        "optimiser_status": str(optimisation.message),
        "optimiser_success": True,
    }


def _fit_beta_binomial_map_with_fallback(
    group_counts,
    mu_prior,
    kappa_prior,
    smoothing_alpha,
    default_kappa,
    min_groups_for_full_bayes,
):
    if len(group_counts) < int(min_groups_for_full_bayes):
        mu = _smoothed_frequency_from_group_counts(group_counts, smoothing_alpha=smoothing_alpha)
        kappa = float(default_kappa)
        return {
            "mu": mu,
            "kappa": kappa,
            "alpha": float(mu * kappa),
            "beta": float((1.0 - mu) * kappa),
            "fit_status": "fallback_too_few_groups",
            "group_count": int(len(group_counts)),
        }

    try:
        fitted = _fit_beta_binomial_map(
            group_counts=group_counts,
            mu_prior=mu_prior,
            kappa_prior=kappa_prior,
            default_kappa=default_kappa,
        )
    except Exception as exc:
        mu = _smoothed_frequency_from_group_counts(group_counts, smoothing_alpha=smoothing_alpha)
        kappa = float(default_kappa)
        return {
            "mu": mu,
            "kappa": kappa,
            "alpha": float(mu * kappa),
            "beta": float((1.0 - mu) * kappa),
            "fit_status": "fallback_optimisation_failed",
            "group_count": int(len(group_counts)),
            "error": str(exc),
        }

    fitted["fit_status"] = "fit"
    fitted["group_count"] = int(len(group_counts))
    return fitted


def _build_grouped_side_counts(sample_alleles, meta, sample_col, lineage_col, group_col, allele, target_lineage, side):
    if side == "target":
        side_meta = meta[meta[lineage_col] == target_lineage].copy()
    elif side == "background":
        side_meta = meta[meta[lineage_col] != target_lineage].copy()
    else:
        raise ValueError(f"Unsupported side {side!r}.")

    grouped_counts = []
    for group_id, group_df in side_meta.groupby(group_col, sort=True):
        sample_ids = list(group_df[sample_col])
        n = len(sample_ids)
        if n == 0:
            continue
        k = sum(1 for sample_id in sample_ids if allele in sample_alleles.get(sample_id, set()))
        grouped_counts.append(
            {
                "group_id": str(group_id),
                "n": int(n),
                "k": int(k),
            }
        )
    return grouped_counts


def _fit_full_bayes_on_model(
    model,
    sample_alleles,
    labels_by_sample,
    meta,
    sample_col,
    lineage_col,
    group_col,
    smoothing_alpha,
    full_bayes_priors,
    min_groups_for_full_bayes,
    default_kappa,
):
    if group_col not in meta.columns:
        raise ValueError(f"full_bayes_group_col {group_col!r} was not found in the training metadata.")

    group_meta = meta[[sample_col, lineage_col, group_col]].dropna().copy()
    group_meta = group_meta[group_meta[sample_col].isin(labels_by_sample)].copy()

    for locus in model["loci"]:
        allele = locus["annotations"]["source_allele_id"]
        target_lineage = locus["target_lineage"]

        target_group_counts = _build_grouped_side_counts(
            sample_alleles=sample_alleles,
            meta=group_meta,
            sample_col=sample_col,
            lineage_col=lineage_col,
            group_col=group_col,
            allele=allele,
            target_lineage=target_lineage,
            side="target",
        )
        background_group_counts = _build_grouped_side_counts(
            sample_alleles=sample_alleles,
            meta=group_meta,
            sample_col=sample_col,
            lineage_col=lineage_col,
            group_col=group_col,
            allele=allele,
            target_lineage=target_lineage,
            side="background",
        )

        target_fit = _fit_beta_binomial_map_with_fallback(
            group_counts=target_group_counts,
            mu_prior=full_bayes_priors["target_mu"],
            kappa_prior=full_bayes_priors["kappa"],
            smoothing_alpha=smoothing_alpha,
            default_kappa=default_kappa,
            min_groups_for_full_bayes=min_groups_for_full_bayes,
        )
        background_fit = _fit_beta_binomial_map_with_fallback(
            group_counts=background_group_counts,
            mu_prior=full_bayes_priors["background_mu"],
            kappa_prior=full_bayes_priors["kappa"],
            smoothing_alpha=smoothing_alpha,
            default_kappa=default_kappa,
            min_groups_for_full_bayes=min_groups_for_full_bayes,
        )

        locus["full_bayes"] = {
            "target_mu": float(target_fit["mu"]),
            "target_kappa": float(target_fit["kappa"]),
            "target_alpha": float(target_fit["alpha"]),
            "target_beta": float(target_fit["beta"]),
            "background_mu": float(background_fit["mu"]),
            "background_kappa": float(background_fit["kappa"]),
            "background_alpha": float(background_fit["alpha"]),
            "background_beta": float(background_fit["beta"]),
            "target_fit_status": target_fit["fit_status"],
            "background_fit_status": background_fit["fit_status"],
            "target_group_count": int(target_fit["group_count"]),
            "background_group_count": int(background_fit["group_count"]),
        }

    model["full_bayes"] = {
        "status": "fit",
        "model_family": "beta_binomial_hierarchical",
        "fit_method": "map_grouped_beta_binomial",
        "group_column": str(group_col),
        "min_groups_for_full_bayes": int(min_groups_for_full_bayes),
        "priors": {
            "target_mu_prior": dict(full_bayes_priors["target_mu"]),
            "background_mu_prior": dict(full_bayes_priors["background_mu"]),
            "kappa_prior": dict(full_bayes_priors["kappa"]),
        },
        "defaults": {
            "target_kappa": float(default_kappa),
            "background_kappa": float(default_kappa),
        },
        "posterior_summary": "map_point_estimates",
        "draw_store": None,
    }
    return model


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
        ## dfs cycle detection to prevent malformed
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
    if current == node:
        return current
    while current in parent_map:
        parent = parent_map.get(current)
        if parent == node:
            return current
        if parent is None:
            break
        current = parent
    raise ValueError(f"Leaf {leaf!r} is not a descendant of node {node!r}.")


def _node_depth(node, parent_map):
    depth = 0
    current = node
    while current in parent_map and parent_map[current] is not None:
        current = parent_map[current]
        depth += 1
    return depth


def _get_descendant_leaves(node, children_by_parent, leaf_set):
    descendants = []
    if node in leaf_set:
        descendants.append(node)
    stack = list(reversed(children_by_parent.get(node, [])))
    while stack:
        current = stack.pop()
        if current in leaf_set:
            descendants.append(current)
            continue
        stack.extend(reversed(children_by_parent.get(current, [])))
    return sorted(set(descendants))


def _get_effective_node_children(node, children_by_parent, leaf_set):
    effective_children = set(children_by_parent.get(node, []))
    if node in leaf_set:
        effective_children.add(node)
    return sorted(effective_children)


def _iter_internal_nodes(parent_map, leaf_set):
    children_by_parent = _build_children_by_parent(parent_map)
    internal_nodes = []
    for node in sorted(children_by_parent):
        if len(_get_effective_node_children(node, children_by_parent, leaf_set)) > 1:
            internal_nodes.append(node)
    return internal_nodes


def _build_node_training_meta(meta, lineage_col, node, parent_map, children_by_parent, leaf_set):
    if node == HIERARCHY_ROOT_NODE:
        ## synthetic root classifier exists only when the tree has multiple disconnected roots
        ## relabels leaves to their top-level branch.
        descendant_leaves = sorted(leaf_set)
        node_meta = meta[meta[lineage_col].isin(descendant_leaves)].copy()
        node_meta["_node_lineage"] = node_meta[lineage_col].map(lambda leaf: _get_root_label(leaf, parent_map))
        node_children = sorted(node_meta["_node_lineage"].dropna().unique())
        return node_meta, descendant_leaves, node_children

    ## each local node model competes only among the nodes direct children
    ## even though training samples from all descendants
    descendant_leaves = _get_descendant_leaves(node, children_by_parent, leaf_set)
    node_meta = meta[meta[lineage_col].isin(descendant_leaves)].copy()
    node_meta["_node_lineage"] = node_meta[lineage_col].map(
        lambda leaf: _get_direct_child_label(leaf, node, parent_map)
    )
    node_children = _get_effective_node_children(node, children_by_parent, leaf_set)
    return node_meta, descendant_leaves, node_children


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
        "schema_version": "1.2",
        "model_id": model_id,
        "species_id": species_id,
        "reference": {
            "name": None,
            "path": None,
            "contigs": sorted({l["chrom"] for l in loci}),
            "ref_alleles_resolved": any(l["ref"] is not None for l in loci),
        },
        "architecture": {
            "feature_family": "empirical",
            "topology": "flat",
            "classification_workflow": "empirical_flat",
            "emission_model": "empirical_bayes",
        },
        "model_type": "empirical_geolineage",
        "hierarchy": {
            "is_hierarchical": False,
            "parent_map": None,
            "children_by_parent": {},
            "root_lineages": [x["lineage_id"] for x in lineage_records],
            "leaf_lineages": [x["lineage_id"] for x in lineage_records],
            "top_classifier_node": None,
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
            "fit_method": None,
            "group_column": None,
            "min_groups_for_full_bayes": None,
            "priors": {
                "mu_prior": None,
                "kappa_prior": None,
            },
            "defaults": None,
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

    ## interleave per lin rankings so one lin cannot monopolise the candidate pool before LOSO
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

    ## IMPORTANT: this is a floor rather than a hard threshold
    ## supported markers are added for underpowered classes
    ## richer lineages are allowed to keep more than the minimum
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
    emission_model="empirical_bayes",
    full_bayes_group_col=None,
    full_bayes_priors=None,
    min_groups_for_full_bayes=3,
    default_kappa=20.0,
    reporting_support_bins=None,
    reporting_support_quantile=0.95,
    reporting_support_min_samples=2,
    collision_min_pair_count=2,
    collision_min_pair_rate=0.05,
    node_name=None,
):
    resolved_group_col = full_bayes_group_col or study_col
    meta_columns = [sample_col, lineage_col, study_col]
    if resolved_group_col not in meta_columns:
        meta_columns.append(resolved_group_col)
    meta = meta_df[meta_columns].dropna().copy()
    lineage_counts = meta[lineage_col].value_counts()
    keep_lineages = set(lineage_counts[lineage_counts >= min_lineage_samples].index)
    meta = meta[meta[lineage_col].isin(keep_lineages)].copy()
    if meta.empty or meta[lineage_col].nunique() < 2:
        return None, []

    ## loso folds are only valid when every held out label is represented in the training split
    ## otherwise that fold cannot score all classes
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
            child_ardal_kwargs={"roaring" : True},
        )
        ard_test = ard_obj.get.subset(
            guids=test_meta[sample_col].tolist(),
            chunk_size=1000,
            threads=12,
            drop_zero_cols=False,
            child_ardal_kwargs={"roaring" : True},
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

        ## evaluate prefix panels in order so loso chooses the smallest strong local model
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

            if emission_model == "full_bayes_map":
                model = _fit_full_bayes_on_model(
                    model=model,
                    sample_alleles=train_sample_alleles,
                    labels_by_sample=y_train,
                    meta=train_meta,
                    sample_col=sample_col,
                    lineage_col=lineage_col,
                    group_col=resolved_group_col,
                    smoothing_alpha=smoothing_alpha,
                    full_bayes_priors=_resolve_full_bayes_priors(full_bayes_priors),
                    min_groups_for_full_bayes=min_groups_for_full_bayes,
                    default_kappa=default_kappa,
                )
                _score_df, prob_df = score_samples_full_bayes_tables(
                    sample_alleles=test_sample_alleles,
                    model=model,
                )
            elif emission_model == "empirical_bayes":
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
            else:
                raise ValueError(f"Unsupported emission_model {emission_model!r}.")

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
        if emission_model == "full_bayes_map":
            best_model = _fit_full_bayes_on_model(
                model=best_model,
                sample_alleles=train_sample_alleles,
                labels_by_sample=y_train,
                meta=train_meta,
                sample_col=sample_col,
                lineage_col=lineage_col,
                group_col=resolved_group_col,
                smoothing_alpha=smoothing_alpha,
                full_bayes_priors=_resolve_full_bayes_priors(full_bayes_priors),
                min_groups_for_full_bayes=min_groups_for_full_bayes,
                default_kappa=default_kappa,
            )

        fold_results.append(best_result)

    if not fold_results:
        return None, []

    ## stable alleles are the fold winning markers that recur across enough loso splits
    ## can be topped up to satisfy the per lin min
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
    ard_full = ard_obj.get.subset(guids=full_ids, chunk_size=1000, threads=12, drop_zero_cols=False, child_ardal_kwargs={"roaring" : True})
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
    if emission_model == "full_bayes_map":
        final_model = _fit_full_bayes_on_model(
            model=final_model,
            sample_alleles=full_sample_alleles,
            labels_by_sample=full_labels,
            meta=meta,
            sample_col=sample_col,
            lineage_col=lineage_col,
            group_col=resolved_group_col,
            smoothing_alpha=smoothing_alpha,
            full_bayes_priors=_resolve_full_bayes_priors(full_bayes_priors),
            min_groups_for_full_bayes=min_groups_for_full_bayes,
            default_kappa=default_kappa,
        )
    elif emission_model != "empirical_bayes":
        raise ValueError(f"Unsupported emission_model {emission_model!r}.")

    final_model["architecture"] = {
        **final_model.get("architecture", {}),
        "feature_family": "empirical",
        "topology": "flat",
        "classification_workflow": "empirical_flat",
        "emission_model": emission_model,
    }

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
        "candidate_strategy": "round_robin_top_k_per_lineage",
        "emission_model": emission_model,
        "full_bayes_group_col": str(resolved_group_col),
        "min_groups_for_full_bayes": int(min_groups_for_full_bayes),
        "default_kappa": float(default_kappa),
        "full_bayes_priors": _resolve_full_bayes_priors(full_bayes_priors),
    }
    final_model["provenance"]["fold_results"] = fold_results
    final_model["summary"]["selected_locus_count_after_loso"] = len(stable_alleles)
    final_model["hierarchy"] = {
        "is_hierarchical": False,
        "parent_map": None,
        "root_lineages": sorted(meta[lineage_col].unique()),
        "children_by_parent": {},
        "leaf_lineages": sorted(meta[lineage_col].unique()),
        "top_classifier_node": None,
        "lineage_naming_rule": None,
        "note": "Geolineages are treated as a flat multiclass set.",
    }
    final_model["node_models"] = {}

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
    emission_model="empirical_bayes",
    full_bayes_group_col=None,
    full_bayes_priors=None,
    min_groups_for_full_bayes=3,
    default_kappa=20.0,
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
        ## flat unless there is an explicit hierarchy
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
            emission_model=emission_model,
            full_bayes_group_col=full_bayes_group_col,
            full_bayes_priors=full_bayes_priors,
            min_groups_for_full_bayes=min_groups_for_full_bayes,
            default_kappa=default_kappa,
            reporting_support_bins=reporting_support_bins,
            reporting_support_quantile=reporting_support_quantile,
            reporting_support_min_samples=reporting_support_min_samples,
            collision_min_pair_count=collision_min_pair_count,
            collision_min_pair_rate=collision_min_pair_rate,
        )
        if final_model is None:
            raise ValueError("No valid LOSO folds were produced.")
        return final_model, fold_results

    resolved_group_col = full_bayes_group_col or study_col
    meta_columns = [sample_col, lineage_col, study_col]
    if resolved_group_col not in meta_columns:
        meta_columns.append(resolved_group_col)
    meta = meta_df[meta_columns].dropna().copy()
    lineage_counts = meta[lineage_col].value_counts()
    keep_lineages = set(lineage_counts[lineage_counts >= min_lineage_samples].index)
    meta = meta[meta[lineage_col].isin(keep_lineages)].copy()
    if meta.empty:
        raise ValueError("No samples remain after applying min_lineage_samples.")

    leaf_set = set(meta[lineage_col].unique())
    normalised_parent_map = _validate_parent_map(parent_map, leaf_set)
    children_by_parent = _build_children_by_parent(normalised_parent_map)
    roots = _find_roots(normalised_parent_map)

    ## train one local classifier per internal node
    ## each node only competes among its direct children
    nodes_to_train = []
    if len(roots) > 1:
        nodes_to_train.append(HIERARCHY_ROOT_NODE)
    for node in _iter_internal_nodes(normalised_parent_map, leaf_set):
        if node not in nodes_to_train:
            nodes_to_train.append(node)

    node_models = {}
    fold_results_by_node = {}
    fold_results = []

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
            emission_model=emission_model,
            full_bayes_group_col=full_bayes_group_col,
            full_bayes_priors=full_bayes_priors,
            min_groups_for_full_bayes=min_groups_for_full_bayes,
            default_kappa=default_kappa,
            reporting_support_bins=reporting_support_bins,
            reporting_support_quantile=reporting_support_quantile,
            reporting_support_min_samples=reporting_support_min_samples,
            collision_min_pair_count=collision_min_pair_count,
            collision_min_pair_rate=collision_min_pair_rate,
            node_name=node,
        )
        if local_model is None:
            continue

        for result in local_fold_results:
            result["node"] = node

        local_model["provenance"]["node"] = {
            "node_name": node,
            "node_children": list(node_children),
            "node_descendant_leaves": list(descendant_leaves),
            "min_lineage_alleles": int(resolved_min_lineage_alleles),
        }
        node_models[node] = local_model
        fold_results_by_node[node] = local_fold_results
        fold_results.extend(local_fold_results)

    if not node_models:
        raise ValueError("No hierarchical node produced a valid model.")

    serialised_parent_map = dict(normalised_parent_map)
    serialised_children = dict(children_by_parent)
    serialised_roots = list(roots)
    top_classifier_node = None

    if len(roots) > 1:
        ## multiple roots need a synthetic top level classifier so inference has an explicit first decision point
        serialised_parent_map[HIERARCHY_ROOT_NODE] = None
        serialised_children[HIERARCHY_ROOT_NODE] = list(roots)
        
        for root in roots:
            serialised_parent_map[root] = HIERARCHY_ROOT_NODE
        serialised_roots = [HIERARCHY_ROOT_NODE]
        top_classifier_node = HIERARCHY_ROOT_NODE
        
    elif len(roots) == 1 and roots[0] in node_models:
        top_classifier_node = roots[0]

    if top_classifier_node is None and node_models:
        top_classifier_node = min(
            node_models,
            key=lambda node: (
                _node_depth(node, normalised_parent_map),
                node,
            ),
        )

    created_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    if top_classifier_node is not None:
        top_model = copy.deepcopy(node_models[top_classifier_node])
        top_model["schema_version"] = "1.2"
        top_model["architecture"] = {
            "feature_family": "empirical",
            "topology": "hierarchical",
            "classification_workflow": "empirical_hierarchical",
            "emission_model": emission_model,
        }
        top_model["model_type"] = "hierarchical_empirical_geolineage"
        top_model["hierarchy"] = {
            "is_hierarchical": True,
            "parent_map": serialised_parent_map,
            "root_lineages": serialised_roots,
            "children_by_parent": serialised_children,
            "leaf_lineages": sorted(leaf_set),
            "top_classifier_node": top_classifier_node,
            "lineage_naming_rule": None,
            "note": "Lineage hierarchy is explicit. Each node model classifies among direct children only.",
        }
        top_model["node_models"] = node_models
        top_model["summary"] = {
            **top_model.get("summary", {}),
            "node_model_count": len(node_models),
            "leaf_lineage_count": len(leaf_set),
            "has_reporting_support": any(
                model.get("reporting_support", {}).get("status") == "fit"
                for model in node_models.values()
            ),
        }
        top_model["provenance"] = {
            **top_model.get("provenance", {}),
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
            "emission_model": emission_model,
            "full_bayes_group_col": str(resolved_group_col),
            "min_groups_for_full_bayes": int(min_groups_for_full_bayes),
            "default_kappa": float(default_kappa),
            "full_bayes_priors": _resolve_full_bayes_priors(full_bayes_priors),
        }
        hierarchical_model = top_model
    else:
        hierarchical_model = {
            "schema_version": "1.2",
            "model_id": model_id,
            "species_id": species_id,
            "architecture": {
                "feature_family": "empirical",
                "topology": "hierarchical",
                "classification_workflow": "empirical_hierarchical",
                "emission_model": emission_model,
            },
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
                "emission_model": emission_model,
                "full_bayes_group_col": str(resolved_group_col),
                "min_groups_for_full_bayes": int(min_groups_for_full_bayes),
                "default_kappa": float(default_kappa),
                "full_bayes_priors": _resolve_full_bayes_priors(full_bayes_priors),
            },
        }

    return hierarchical_model, fold_results


def write_model_json(model, output_json):
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as handle:
        json.dump(model, handle, indent=2)
    return output_json

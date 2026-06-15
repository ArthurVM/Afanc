import importlib.util
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts" / "loso_sparse_model.py"
if not MODULE_PATH.exists():
    pytest.skip("loso_sparse_model.py is not distributed with this tree", allow_module_level=True)
SPEC = importlib.util.spec_from_file_location("loso_sparse_model", MODULE_PATH)
loso_sparse_model = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(loso_sparse_model)


class FakeAlleleSubset:
    def __init__(self, sample_alleles, alleles):
        self._sample_alleles = sample_alleles
        self._alleles = set(alleles)
        self.io = self

    def to_dict(self, backend="auto"):
        return {
            sample_id: sorted(alleles & self._alleles)
            for sample_id, alleles in self._sample_alleles.items()
        }


class FakeSubsetGetter:
    def __init__(self, sample_alleles):
        self._sample_alleles = sample_alleles

    def subset(self, guids=None, alleles=None, **kwargs):
        if guids is not None:
            return FakeArdSubset({guid: self._sample_alleles[guid] for guid in guids})
        if alleles is not None:
            return FakeAlleleSubset(self._sample_alleles, alleles)
        raise ValueError("Expected either guids or alleles.")


class FakeStats:
    def __init__(self, sample_alleles):
        self._sample_alleles = sample_alleles

    def allele_inform(self, positives, method="kullbackleibler"):
        positives = set(positives)
        negatives = set(self._sample_alleles) - positives
        all_alleles = sorted(set().union(*self._sample_alleles.values()))
        scores = {}

        for allele in all_alleles:
            pos_rate = sum(allele in self._sample_alleles[s] for s in positives) / max(1, len(positives))
            neg_rate = sum(allele in self._sample_alleles[s] for s in negatives) / max(1, len(negatives))
            scores[allele] = pos_rate - neg_rate

        ordered = sorted(scores.items(), key=lambda item: (item[1], item[0]), reverse=True)
        return dict(ordered)

    def allele_cooc(self, alleles, threshold=None, threads=1):
        return {allele: [] for allele in alleles}


class FakeArdSubset:
    def __init__(self, sample_alleles):
        self._sample_alleles = {sample_id: set(alleles) for sample_id, alleles in sample_alleles.items()}
        self.get = FakeSubsetGetter(self._sample_alleles)
        self.stats = FakeStats(self._sample_alleles)


class FakeArdObj:
    def __init__(self, sample_alleles):
        self._sample_alleles = {sample_id: set(alleles) for sample_id, alleles in sample_alleles.items()}
        self.get = FakeSubsetGetter(self._sample_alleles)


def _make_fake_training_data():
    sample_alleles = {
        "s1": {"a1", "a1b", "a_root"},
        "s2": {"a2", "a2b", "a_root"},
        "s3": {"b1", "b1b"},
        "s4": {"a1", "a1b", "a_root"},
        "s5": {"a2", "a2b", "a_root"},
        "s6": {"b1", "b1b"},
    }
    meta = pd.DataFrame(
        [
            {"Sample": "s1", "GeoLineage": "A1", "Country": "C1"},
            {"Sample": "s2", "GeoLineage": "A2", "Country": "C1"},
            {"Sample": "s3", "GeoLineage": "B1", "Country": "C1"},
            {"Sample": "s4", "GeoLineage": "A1", "Country": "C2"},
            {"Sample": "s5", "GeoLineage": "A2", "Country": "C2"},
            {"Sample": "s6", "GeoLineage": "B1", "Country": "C2"},
        ]
    )
    return FakeArdObj(sample_alleles), meta


def test_normalise_parent_map_adds_missing_root_nodes():
    parent_map = {"A1": "A", "A2": "A", "B1": "B"}
    normalised = loso_sparse_model._validate_parent_map(parent_map, {"A1", "A2", "B1"})
    assert normalised["A"] is None
    assert normalised["B"] is None


def test_flat_min_lineage_alleles_is_a_minimum():
    ard_obj, meta = _make_fake_training_data()
    model, fold_results = loso_sparse_model.build_geolineage_min_model_sparse(
        ard_obj=ard_obj,
        meta_df=meta,
        lineage_col="GeoLineage",
        study_col="Country",
        top_k_per_lineage=4,
        max_model_alleles=10,
        min_lineage_samples=1,
        min_lineage_alleles=2,
        uniform_priors=True,
    )

    lineage_counts = {lineage["lineage_id"]: lineage["direct_marker_count"] for lineage in model["lineages"]}
    assert set(lineage_counts) == {"A1", "A2", "B1"}
    assert all(count >= 2 for count in lineage_counts.values())
    assert len(fold_results) == 2


def test_hierarchical_parent_map_trains_node_local_models():
    ard_obj, meta = _make_fake_training_data()
    model, fold_results_by_node = loso_sparse_model.build_geolineage_min_model_sparse(
        ard_obj=ard_obj,
        meta_df=meta,
        lineage_col="GeoLineage",
        study_col="Country",
        top_k_per_lineage=4,
        max_model_alleles=10,
        min_lineage_samples=1,
        min_lineage_alleles=2,
        uniform_priors=True,
        parent_map={"A1": "A", "A2": "A", "B1": "B"},
    )

    assert model["model_type"] == "hierarchical_empirical_geolineage"
    assert set(model["node_models"]) == {loso_sparse_model.HIERARCHY_ROOT_NODE, "A"}
    assert set(fold_results_by_node) == {loso_sparse_model.HIERARCHY_ROOT_NODE, "A"}

    root_model = model["node_models"][loso_sparse_model.HIERARCHY_ROOT_NODE]
    assert {lineage["lineage_id"] for lineage in root_model["lineages"]} == {"A", "B"}

    a_model = model["node_models"]["A"]
    assert {lineage["lineage_id"] for lineage in a_model["lineages"]} == {"A1", "A2"}
    assert a_model["provenance"]["node"]["min_lineage_alleles"] == 2


def test_fit_beta_binomial_map_with_fallback_uses_default_kappa_for_too_few_groups():
    fit = loso_sparse_model._fit_beta_binomial_map_with_fallback(
        group_counts=[{"group_id": "C1", "n": 5, "k": 4}],
        mu_prior={"alpha": 7.0, "beta": 5.0},
        kappa_prior={"alpha": 2.0, "beta": 1.0},
        smoothing_alpha=1.0,
        default_kappa=17.0,
        min_groups_for_full_bayes=2,
    )

    assert fit["fit_status"] == "fallback_too_few_groups"
    assert fit["kappa"] == 17.0
    assert fit["group_count"] == 1


def test_flat_full_bayes_map_populates_full_bayes_parameters():
    ard_obj, meta = _make_fake_training_data()
    model, fold_results = loso_sparse_model.build_geolineage_min_model_sparse(
        ard_obj=ard_obj,
        meta_df=meta,
        lineage_col="GeoLineage",
        study_col="Country",
        top_k_per_lineage=4,
        max_model_alleles=10,
        min_lineage_samples=1,
        min_lineage_alleles=2,
        uniform_priors=True,
        emission_model="full_bayes_map",
        min_groups_for_full_bayes=2,
        default_kappa=11.0,
    )

    assert model["full_bayes"]["status"] == "fit"
    assert model["full_bayes"]["fit_method"] == "map_grouped_beta_binomial"
    assert model["provenance"]["loso"]["emission_model"] == "full_bayes_map"
    assert len(fold_results) == 2

    for locus in model["loci"]:
        assert "full_bayes" in locus
        assert "target_mu" in locus["full_bayes"]
        assert "target_kappa" in locus["full_bayes"]
        assert "background_mu" in locus["full_bayes"]
        assert "background_kappa" in locus["full_bayes"]
        assert locus["full_bayes"]["target_group_count"] >= 2

import json

from Afanc.screen.profiles import load_profiles_manifest, resolve_profile_for_event, validate_profile_reference
from Afanc.screen.report.final_report import _index_detection_event_by_taxid, makeFinalReport
from Afanc.screen.variant_profiler.bayesian_profile import run_lineage_classification


def _write_model(model_path, ref="A"):
    model = {
        "schema_version": "model_schema_v1",
        "model_id": "toy",
        "model_type": "empirical_loso",
        "architecture": {
            "feature_family": "empirical",
            "topology": "flat",
            "emission_model": "empirical_bayes",
        },
        "lineages": [
            {"lineage_id": "L1", "prior": 0.5},
            {"lineage_id": "L2", "prior": 0.5},
        ],
        "loci": [
            {
                "allele_id": f"chr1.2.{ref}.G",
                "chrom": "chr1",
                "pos": 2,
                "bed_start": 1,
                "bed_end": 2,
                "ref": ref,
                "alt": "G",
                "locus_id": "chr1.2.A.G",
                "target_lineage": "L1",
                "emission": {"target_frequency": 0.99, "background_frequency": 0.01},
            }
        ],
    }
    model_path.write_text(json.dumps(model))


def test_validate_profile_reference_checks_model_ref_alleles(tmp_path):
    profiles_dir = tmp_path / "profiles"
    profile_dir = profiles_dir / "Toy"
    profile_dir.mkdir(parents=True)
    reference = profile_dir / "ref.fa"
    model = profile_dir / "model.json"
    reference.write_text(">chr1\nAAC\n")
    _write_model(model, ref="A")

    validation = validate_profile_reference(
        {"reference": "Toy/ref.fa", "model": "Toy/model.json"},
        str(profiles_dir),
    )

    assert validation["valid"] is True
    assert validation["checked_loci"] == 1


def test_validate_profile_reference_rejects_ref_mismatch(tmp_path):
    profiles_dir = tmp_path / "profiles"
    profile_dir = profiles_dir / "Toy"
    profile_dir.mkdir(parents=True)
    reference = profile_dir / "ref.fa"
    model = profile_dir / "model.json"
    reference.write_text(">chr1\nAAC\n")
    _write_model(model, ref="T")

    validation = validate_profile_reference(
        {"reference": "Toy/ref.fa", "model": "Toy/model.json"},
        str(profiles_dir),
    )

    assert validation["valid"] is False
    assert "expects T" in validation["errors"][0]


def test_resolve_profile_prefers_detected_taxon_over_closest_variant():
    manifest = {
        "profiles": {
            "1773": {"name": "Mycobacterium tuberculosis", "enabled": True},
            "123": {"name": "variant", "enabled": True},
        }
    }
    event = {
        "taxon_id": "1773",
        "closest_variant": {"taxon_id": "123"},
    }

    profile = resolve_profile_for_event(event, manifest)

    assert profile["taxon_id"] == "1773"


def test_resolve_profile_matches_renamed_mpox_clade_by_alias_contains():
    manifest = {
        "profiles": {
            "10244": {
                "name": "Monkeypox virus",
                "aliases": ["Mpox_IIb", "Monkeypox virus", "Mpox", "MPXV", "clade IIb"],
                "enabled": True,
            }
        }
    }
    event = {
        "taxon_id": 3431483,
        "name": "Orthopoxvirus monkeypox",
        "closest_variant": {
            "taxon_id": 3706793,
            "name": "Monkeypox virus cladeIIb",
        },
    }

    profile = resolve_profile_for_event(event, manifest)

    assert profile["taxon_id"] == "10244"
    assert profile["profile_match"] == {
        "type": "alias_contains",
        "matched_value": "Monkeypox virus cladeIIb",
    }


def test_run_lineage_classification_writes_summary(tmp_path):
    model = tmp_path / "model.json"
    snps = tmp_path / "sample.snps.json"
    _write_model(model, ref="A")
    snps.write_text(json.dumps({"alleles": ["chr1.2.A.G"], "missing": []}))

    class Args:
        reportsDir = tmp_path
        output_prefix = "sample"
        lineage_min_support = 1
        lineage_min_support_fraction = 0.75
        lineage_min_callable_fraction = 0.5
        lineage_ambiguity_margin = 0.1
        lineage_disable_incomplete_descent = False
        lineage_disable_reference_marker_inference = False
        lineage_tie_delta = 0.0

    result = run_lineage_classification(
        Args,
        {"ACC": {"snp_json": str(snps)}},
        {
            "ACC": {
                "lineage_profile": {
                    "taxon_id": "1",
                    "name": "Toy",
                    "model": str(model),
                    "reference": "ref.fa",
                }
            }
        },
    )

    assert result["ACC"]["best_lineage"] == "L1"
    assert (tmp_path / "ACC.lineage_classification.json").exists()


def test_run_lineage_classification_reports_no_profile(tmp_path):
    snps = tmp_path / "sample.snps.json"
    snps.write_text(json.dumps({"alleles": [], "missing": []}))

    class Args:
        reportsDir = tmp_path
        output_prefix = "sample"

    result = run_lineage_classification(
        Args,
        {"ACC": {"snp_json": str(snps), "snp_count": 0, "missing_count": 1}},
        {"ACC": {}},
    )

    assert result["ACC"]["status"] == "not_run"
    assert result["ACC"]["reason"] == "no_profile_model"


def test_load_profiles_manifest_finds_sibling_profiles_dir(tmp_path):
    database = tmp_path / "DB"
    profiles = tmp_path / "profiles"
    database.mkdir()
    profiles.mkdir()
    (profiles / "profiles.json").write_text(json.dumps({"profiles": {}}))

    profiles_dir, manifest = load_profiles_manifest(str(database))

    assert profiles_dir == str(profiles)
    assert manifest == {"profiles": {}}


def test_detection_events_are_indexed_by_terminal_clustering_call_with_profiles():
    event = {
        "name": "Mycobacterium tuberculosis",
        "taxon_id": 1773,
        "accession": "ASM19595v2",
        "assembly": "H37Rv.fasta",
        "closest_variant": {
            "name": "Mycobacterium tuberculosis variant bovis BCG",
            "taxon_id": 33892,
            "accession": "ASM23472v1",
            "assembly": "33892_GCA_000234725.1_ASM23472v1_genomic.fna",
        },
    }
    snp_profile = {
        "H37Rv": {"snp_count": 2},
        "33892_GCA_000234725.1_ASM23472v1": {"snp_count": 3},
    }
    lineage_profile = {
        "H37Rv": {"best_lineage": "lineage4"},
    }

    by_taxid = {}
    _index_detection_event_by_taxid(
        by_taxid,
        event,
        snp_profile=snp_profile,
        lineage_profile=lineage_profile,
    )

    assert sorted(by_taxid) == ["33892"]
    assert by_taxid["33892"]["name"] == "Mycobacterium tuberculosis variant bovis BCG"
    assert by_taxid["33892"]["clustering"]["parent_context"]["taxon_id"] == 1773
    assert by_taxid["33892"]["snp_profile"]["H37Rv"]["snp_count"] == 2
    assert by_taxid["33892"]["lineage_profile"]["H37Rv"]["best_lineage"] == "lineage4"
    assert by_taxid["33892"]["snp_profile"]["33892_GCA_000234725.1_ASM23472v1"]["snp_count"] == 3


def test_mpox_parent_and_clade_are_not_duplicated_in_final_report_index():
    event = {
        "name": "Orthopoxvirus monkeypox",
        "taxon_id": 3431483,
        "accession": "NA",
        "assembly": "NC_063383.fasta",
        "closest_variant": {
            "name": "Monkeypox virus cladeIIb",
            "taxon_id": 3706793,
            "accession": "",
            "assembly": "NC_063383.fasta",
        },
    }
    snp_profile = {"NC_063383": {"snp_count": 66}}
    lineage_profile = {"NC_063383": {"best_lineage": "B.1.5"}}

    by_taxid = {}
    _index_detection_event_by_taxid(
        by_taxid,
        event,
        snp_profile=snp_profile,
        lineage_profile=lineage_profile,
    )

    assert sorted(by_taxid) == ["3706793"]
    assert by_taxid["3706793"]["name"] == "Monkeypox virus cladeIIb"
    assert by_taxid["3706793"]["clustering"]["parent_context"]["taxon_id"] == 3431483
    assert by_taxid["3706793"]["snp_profile"]["NC_063383"]["snp_count"] == 66
    assert by_taxid["3706793"]["lineage_profile"]["NC_063383"]["best_lineage"] == "B.1.5"


def test_final_report_is_taxid_indexed_and_writes_subreports(tmp_path, monkeypatch):
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir()
    (reports_dir / "sample.k2.json").write_text(
        json.dumps(
            {
                "Detection_events": [
                    {
                        "name": "Mycobacterium tuberculosis",
                        "taxon_id": 1773,
                        "accession": "ASM19595v2",
                        "assembly": "H37Rv.fasta",
                    }
                ]
            }
        )
    )
    (reports_dir / "H37Rv.mapstats.json").write_text(
        json.dumps(
            {
                "warnings": [],
                "map_data": {
                    "mean_DOC": 12,
                    "median_DOC": 11,
                    "proportion_cov": 0.9,
                    "gini": 0.1,
                },
            }
        )
    )

    class Args:
        database = "db"
        fastq = ["reads_1.fq", "reads_2.fq"]
        num_threshold = 10
        pct_threshold = 0.1
        upper_bound = 1.0
        lower_bound = 0.0
        variant_caller = "freebayes"
        lineage_profile_compound = False
        snp_min_qual = 20
        snp_min_dp = 5
        snp_min_missing_depth = 5
        snp_allow_filtered = False
        snp_accept_missing_qual = False
        no_lineage_classify = False
        profiles_dir = None
        lineage_min_support = 1
        lineage_min_support_fraction = 0.75
        lineage_min_callable_fraction = 0.5
        lineage_ambiguity_margin = 0.1
        lineage_tie_delta = 0
        lineage_disable_reference_marker_inference = False
        lineage_disable_incomplete_descent = False
        output_prefix = "sample"
        reportsDir = str(reports_dir)

    monkeypatch.chdir(tmp_path)
    final_report = makeFinalReport(
        Args,
        reports={"H37Rv": str(reports_dir / "H37Rv.mapstats.json")},
        snp_profile={"H37Rv": {"snp_count": 1}},
        lineage_profile={"H37Rv": {"best_lineage": "lineage4"}},
    )

    data = json.loads((tmp_path / "sample.json").read_text())
    detection_events = data["results"]["Detection_events"]

    assert final_report == str(tmp_path / "sample.json")
    assert sorted(detection_events) == ["1773"]
    assert detection_events["1773"]["name"] == "Mycobacterium tuberculosis"
    assert detection_events["1773"]["clustering"]["mean_DOC"] == 12
    assert detection_events["1773"]["snp_profile"]["H37Rv"]["snp_count"] == 1
    assert detection_events["1773"]["lineage_profile"]["H37Rv"]["best_lineage"] == "lineage4"
    assert "Clustering_results" not in detection_events
    assert "SNP_profile" not in detection_events
    assert (reports_dir / "sample.clustering_results.json").exists()
    assert (reports_dir / "sample.snp_profile.json").exists()
    assert (reports_dir / "sample.lineage_profile.json").exists()

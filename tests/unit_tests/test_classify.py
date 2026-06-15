import json
from io import StringIO
from types import SimpleNamespace

import pytest

from Afanc.classify.allele_json import load_allele_json
from Afanc.classify.reference_validation import validate_allele_reference, validate_vcf_reference
from Afanc.classify.runFuncs import runClassify
from Afanc.parser import base_parser
from Afanc.screen.profiles import load_profiles_manifest, resolve_profile_for_species


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
                "locus_id": f"chr1.2.{ref}.G",
                "target_lineage": "L1",
                "emission": {"target_frequency": 0.99, "background_frequency": 0.01},
            }
        ],
    }
    model_path.write_text(json.dumps(model))


def _make_profiles(tmp_path):
    profiles_dir = tmp_path / "profiles"
    profile_dir = profiles_dir / "Toy"
    profile_dir.mkdir(parents=True)
    (profile_dir / "ref.fa").write_text(">chr1\nAAC\n")
    _write_model(profile_dir / "model.json")
    (profiles_dir / "profiles.json").write_text(
        json.dumps(
            {
                "profiles": {
                    "1": {
                        "name": "Toy species",
                        "aliases": ["Toy alias"],
                        "enabled": True,
                        "reference": "Toy/ref.fa",
                        "model": "Toy/model.json",
                    }
                }
            }
        )
    )
    return profiles_dir


def _make_args(tmp_path, profiles_dir, **overrides):
    args = SimpleNamespace(
        species="Toy species",
        profiles_dir=str(profiles_dir),
        database=str(tmp_path),
        output_prefix=str(tmp_path / "sample"),
        sample_id="sample",
        vcf=None,
        allele_json=None,
        allele_id_format=None,
        depth_bed=None,
        snp_min_qual=30.0,
        snp_min_dp=None,
        snp_min_missing_depth=10,
        snp_allow_filtered=False,
        snp_accept_missing_qual=False,
        lineage_min_support=1,
        lineage_min_support_fraction=0.75,
        lineage_min_callable_fraction=0.5,
        lineage_ambiguity_margin=0.1,
        lineage_tie_delta=0.0,
        lineage_disable_reference_marker_inference=False,
        lineage_disable_incomplete_descent=False,
        stdout=StringIO(),
        stderr=StringIO(),
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def test_resolve_profile_for_species_matches_alias(tmp_path):
    profiles_dir = _make_profiles(tmp_path)
    _profiles_dir, manifest = load_profiles_manifest(str(tmp_path), profiles_dir=str(profiles_dir))

    profile = resolve_profile_for_species("Toy alias", manifest)

    assert profile["taxon_id"] == "1"
    assert profile["profile_match"]["type"] == "name"


def test_vcf_reference_validation_rejects_length_mismatch(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAC\n")
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=4>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
            ]
        )
        + "\n"
    )

    validation = validate_vcf_reference(vcf, reference)

    assert validation["valid"] is False
    assert "length is 4" in validation["errors"][0]


def test_vcf_reference_validation_rejects_ref_mismatch(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAC\n")
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=3>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                "chr1\t2\t.\tT\tG\t99\tPASS\tDP=12\tGT:DP\t1/1:12",
            ]
        )
        + "\n"
    )

    validation = validate_vcf_reference(vcf, reference)

    assert validation["valid"] is False
    assert "VCF REF allele" in validation["errors"][0]


def test_allele_reference_validation_requires_ref_field(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAC\n")

    validation = validate_allele_reference(["chr1.2.G"], "{chrom}.{pos}.{alt}", reference)

    assert validation["valid"] is False
    assert any("ref" in error for error in validation["errors"])


def test_allele_reference_validation_rejects_ref_mismatch(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nAAC\n")

    validation = validate_allele_reference(["chr1.2.T.G"], "{chrom}.{pos}.{ref}.{alt}", reference)

    assert validation["valid"] is False
    assert "expects ref T" in validation["errors"][0]


def test_load_allele_json_accepts_ardal_allele_field(tmp_path):
    allele_json = tmp_path / "alleles.json"
    allele_json.write_text(json.dumps({"allele": ["chr1.2.A.G"], "missing": ["chr1.3"]}))

    payload = load_allele_json(allele_json)

    assert payload == {"alleles": ["chr1.2.A.G"], "missing": [["chr1", "3"]]}


def test_run_classify_from_allele_json_writes_summary(tmp_path, capsys):
    profiles_dir = _make_profiles(tmp_path)
    allele_json = tmp_path / "alleles.json"
    allele_json.write_text(json.dumps({"allele": ["chr1.2.A.G"], "missing": []}))

    args = _make_args(
        tmp_path,
        profiles_dir,
        allele_json=str(allele_json),
        allele_id_format="{chrom}.{pos}.{ref}.{alt}",
    )

    runClassify(args)

    summary = json.loads((tmp_path / "sample.classify.json").read_text())
    assert summary["classification"]["status"] == "classified"
    assert summary["classification"]["best_lineage"] == "L1"
    assert summary["input"]["type"] == "allele_json"
    assert "Running Bayesian lineage classification" in capsys.readouterr().out
    assert "Running Bayesian lineage classification" in args.stdout.getvalue()


def test_run_classify_from_vcf_writes_summary(tmp_path):
    profiles_dir = _make_profiles(tmp_path)
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=3>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                "chr1\t2\t.\tA\tG\t99\tPASS\tDP=12\tGT:DP\t1/1:12",
            ]
        )
        + "\n"
    )
    args = _make_args(tmp_path, profiles_dir, vcf=str(vcf))

    runClassify(args)

    summary = json.loads((tmp_path / "sample.classify.json").read_text())
    assert summary["classification"]["best_lineage"] == "L1"
    assert summary["input"]["type"] == "vcf"


def test_run_classify_wrong_species_reports_supported_names(tmp_path):
    profiles_dir = _make_profiles(tmp_path)
    allele_json = tmp_path / "alleles.json"
    allele_json.write_text(json.dumps({"allele": ["chr1.2.A.G"], "missing": []}))
    args = _make_args(
        tmp_path,
        profiles_dir,
        species="Wrong species",
        allele_json=str(allele_json),
        allele_id_format="{chrom}.{pos}.{ref}.{alt}",
    )

    with pytest.raises(ValueError) as exc:
        runClassify(args)

    assert "Supported reference species/profile names" in str(exc.value)
    assert "Toy species" in str(exc.value)
    assert "Toy alias" in str(exc.value)
    assert "Supported reference species/profile names" in args.stderr.getvalue()


def test_classify_parser_accepts_allele_json(tmp_path):
    profiles_dir = _make_profiles(tmp_path)
    allele_json = tmp_path / "alleles.json"
    allele_json.write_text(json.dumps({"allele": [], "missing": []}))

    args = base_parser.parse_args(
        [
            "classify",
            "--species",
            "Toy species",
            "--profiles-dir",
            str(profiles_dir),
            "--allele-json",
            str(allele_json),
            "--allele-id-format",
            "{chrom}.{pos}.{ref}.{alt}",
        ]
    )

    assert args.command == "classify"
    assert args.allele_json == allele_json.resolve()
    assert args.allele_id_format == "{chrom}.{pos}.{ref}.{alt}"

import importlib.util
from argparse import Namespace
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_canonical_model.py"
SPEC = importlib.util.spec_from_file_location("build_canonical_model", SCRIPT_PATH)
build_canonical_model = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(build_canonical_model)


def _make_args(input_bed):
    return Namespace(
        input_bed=input_bed,
        output_json=None,
        species_id="Mycobacterium_tuberculosis",
        model_id=None,
        reference_name=None,
        reference_fasta=None,
        target_frequency=0.99,
        background_frequency=0.01,
        allele_id_format="{chrom}.{start}.{ref}.{alt}",
    )


def test_build_model_accepts_five_column_barcode(tmp_path):
    barcode = tmp_path / "canonical.tsv"
    barcode.write_text("Chromosome\t99\t100\tL1\tG\n")

    model, output_json = build_canonical_model.build_model(_make_args(barcode))

    assert output_json == tmp_path / "canonical.canonical_model.json"
    assert model["schema_version"] == "1.2"
    assert model["architecture"] == {
        "feature_family": "canonical",
        "topology": "flat",
        "classification_workflow": "canonical_flat",
        "emission_model": "none",
    }
    assert model["provenance"]["input_format"] == "canonical_barcode_bed_5plus"
    assert model["summary"]["locus_count"] == 1
    assert model["reporting_support"] == {
        "status": "not_fit",
        "method": None,
        "evidence_metric": "usable_loci",
        "base_thresholds": None,
        "collision_map": None,
        "thresholds": None,
        "fallback_policy": None,
        "provenance": None,
    }

    locus = model["loci"][0]
    assert locus["chrom"] == "Chromosome"
    assert locus["bed_start"] == 99
    assert locus["bed_end"] == 100
    assert locus["target_lineage"] == "L1"
    assert locus["alt"] == "G"
    assert locus["annotations"]["clade_label"] is None
    assert locus["annotations"]["spoligotype"] is None
    assert locus["annotations"]["region_of_difference"] is None
    assert locus["annotations"]["source_line_number"] == 1


def test_build_model_preserves_optional_legacy_annotations(tmp_path):
    barcode = tmp_path / "canonical_legacy.tsv"
    barcode.write_text("Chromosome\t199\t200\tL2\tA\tclade2\tspol2\trd2\textra\n")

    model, _ = build_canonical_model.build_model(_make_args(barcode))

    locus = model["loci"][0]
    assert locus["annotations"]["clade_label"] == "clade2"
    assert locus["annotations"]["spoligotype"] == "spol2"
    assert locus["annotations"]["region_of_difference"] == "rd2"


def test_build_model_decodes_encoded_snp_ids_when_bed_columns_are_not_coordinates(tmp_path):
    barcode = tmp_path / "geolineage_sparse.tsv"
    barcode.write_text("Pf3D7_02_v3.305438.T.G\t-1\t0\tSA\tN\n")

    model, _ = build_canonical_model.build_model(_make_args(barcode))

    locus = model["loci"][0]
    assert locus["chrom"] == "Pf3D7_02_v3"
    assert locus["pos"] == 305438
    assert locus["bed_start"] == 305437
    assert locus["bed_end"] == 305438
    assert locus["ref"] == "T"
    assert locus["alt"] == "G"
    assert locus["target_lineage"] == "SA"
    assert locus["locus_id"] == "Pf3D7_02_v3:305438:T>G:SA"

import pytest
from pathlib import Path
from os import path
from Afanc.parser import base_parser
from Afanc._version import __version__

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"
EMPTY_TEST_FILE = TEST_DATA_DIR / "empty_test_file.txt"


## TODO : mocking argparse doesnt work as expected
# @pytest.fixture
# def mock_file_checks(mocker):
#     mocker.patch("Afanc.parser.isFile", return_value="mocked_file_path")
#     mocker.patch("Afanc.parser.isDir", return_value="mocked_dir_path")
#     mocker.patch("Afanc.parser.checkNcbiTaxDB", return_value="mocked_taxdb_path")


def test_parser_version():
    """Test the version flag."""
    with pytest.raises(SystemExit) as e:
        base_parser.parse_args(["-v"])
    assert e.value.code == 0 ## successful exit
    

def test_parser_get_dataset_defaults():
    """Test get_dataset subcommand with default arguments."""
    args = base_parser.parse_args(["get_dataset", str(EMPTY_TEST_FILE)])
    assert args.command == "get_dataset"
    assert args.ID_file == EMPTY_TEST_FILE.resolve()
    assert args.accessions is False
    assert args.num_assemblies == 1
    assert args.output_prefix == "assemblies"
    assert hasattr(args, "func")


def test_parser_get_dataset_custom():
    """Test get_dataset subcommand with custom arguments."""
    args = base_parser.parse_args([
        "get_dataset", str(EMPTY_TEST_FILE),
        "-a",
        "-n", "10",
        "-o", "my_custom_assemblies"
    ])
    assert args.command == "get_dataset"
    assert args.ID_file == EMPTY_TEST_FILE.resolve()
    assert args.accessions is True
    assert args.num_assemblies == 10
    assert args.output_prefix == "my_custom_assemblies"


def test_parser_autodatabase_defaults():
    """Test autodatabase subcommand with default arguments."""
    args = base_parser.parse_args(["autodatabase", str(TEST_DATA_DIR)])
    assert args.command == "autodatabase"
    assert args.fastaDir == TEST_DATA_DIR.resolve()
    assert args.output_prefix == "Afanc_autodb"
    assert args.ncbi_date == "2026-05-01"
    assert args.mode_range == 0.1
    assert args.variant_index_method == "mash"
    assert args.ncbi_tax_db is False
    assert args.threads == 4
    assert args.clean is False
    assert args.superclean is False
    assert hasattr(args, "func")


def _snp_filter_values(args):
    return {
        "snp_min_qual": args.snp_min_qual,
        "snp_min_dp": args.snp_min_dp,
        "snp_min_missing_depth": args.snp_min_missing_depth,
        "snp_allow_filtered": args.snp_allow_filtered,
        "snp_accept_missing_qual": args.snp_accept_missing_qual,
    }


def test_screen_and_classify_share_snp_filter_defaults(tmp_path):
    database = tmp_path / "db"
    database.mkdir()
    r1 = tmp_path / "reads_1.fq"
    r2 = tmp_path / "reads_2.fq"
    vcf = tmp_path / "sample.vcf"
    r1.write_text("")
    r2.write_text("")
    vcf.write_text("")

    screen_args = base_parser.parse_args(["screen", str(database), str(r1), str(r2)])
    classify_args = base_parser.parse_args(["classify", "--species", "Toy species", "--vcf", str(vcf)])

    assert _snp_filter_values(screen_args) == _snp_filter_values(classify_args)


def test_screen_and_classify_parse_snp_filter_overrides_identically(tmp_path):
    database = tmp_path / "db"
    database.mkdir()
    r1 = tmp_path / "reads_1.fq"
    r2 = tmp_path / "reads_2.fq"
    vcf = tmp_path / "sample.vcf"
    r1.write_text("")
    r2.write_text("")
    vcf.write_text("")
    snp_options = [
        "--snp-min-qual", "12",
        "--snp-min-dp", "7",
        "--snp-min-missing-depth", "3",
        "--snp-allow-filtered",
        "--snp-accept-missing-qual",
    ]

    screen_args = base_parser.parse_args(["screen", str(database), str(r1), str(r2), *snp_options])
    classify_args = base_parser.parse_args(["classify", "--species", "Toy species", "--vcf", str(vcf), *snp_options])

    assert _snp_filter_values(screen_args) == _snp_filter_values(classify_args)

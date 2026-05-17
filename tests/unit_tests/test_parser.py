import pytest
from pathlib import Path
from os import path
from Afanc.parser import base_parser
from Afanc._version import __version__


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
    args = base_parser.parse_args(["get_dataset", "./test_data/empty_test_file.txt"])
    assert args.command == "get_dataset"
    assert args.ID_file == Path(path.abspath("./test_data/empty_test_file.txt"))
    assert args.accessions is False
    assert args.num_assemblies == 1
    assert args.output_prefix == "assemblies"
    assert hasattr(args, "func")


def test_parser_get_dataset_custom():
    """Test get_dataset subcommand with custom arguments."""
    args = base_parser.parse_args([
        "get_dataset", "./test_data/empty_test_file.txt",
        "-a",
        "-n", "10",
        "-o", "my_custom_assemblies"
    ])
    assert args.command == "get_dataset"
    assert args.ID_file == Path(path.abspath("./test_data/empty_test_file.txt"))
    assert args.accessions is True
    assert args.num_assemblies == 10
    assert args.output_prefix == "my_custom_assemblies"


def test_parser_autodatabase_defaults():
    """Test autodatabase subcommand with default arguments."""
    args = base_parser.parse_args(["autodatabase", "test_data"])
    assert args.command == "autodatabase"
    assert args.fastaDir == Path(path.abspath("test_data"))
    assert args.output_prefix == "Afanc_autodb"
    assert args.ncbi_date == "2022-05-01"
    assert args.mode_range == 0.1
    assert args.variant_index_method == "mash"
    assert args.ncbi_tax_db is False
    assert args.threads == 4
    assert args.clean is False
    assert args.superclean is False
    assert hasattr(args, "func")

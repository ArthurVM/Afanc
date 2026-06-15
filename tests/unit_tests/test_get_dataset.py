import json
from io import StringIO
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from Afanc.get_dataset.download_assemblies import download_genome, parse_names_file, runGet_dataset


def _args(**overrides):
    args = SimpleNamespace(
        accessions=False,
        num_assemblies=2,
        output_prefix="assemblies",
        ID_file=None,
        stdout=StringIO(),
        stderr=StringIO(),
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def test_parse_names_file_strips_whitespace_and_ignores_blank_lines(tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("\n  Mycobacterium tuberculosis  \n\nMycobacterium avium\n")

    assert parse_names_file(names) == ["Mycobacterium tuberculosis", "Mycobacterium avium"]


def test_download_genome_rejects_deprecated_accession_mode(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError, match="deprecated"):
        download_genome("GCA_000001", _args(accessions=True))


def test_download_genome_selects_highest_n50_accessions_and_moves_fastas(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    summary_payload = {
        "total_count": 3,
        "reports": [
            _report("GCA_LOW", 10),
            _report("GCA_HIGH", 30),
            _report("GCA_MID", 20),
        ],
    }
    command_calls = []

    class FakeCommand:
        def __init__(self, command, subprocess_id):
            self.command = command
            self.subprocess_id = subprocess_id
            command_calls.append((command, subprocess_id))

        def run_comm_quiet(self, *_args, **_kwargs):
            if self.subprocess_id == "DOWNLOAD_HITS":
                return json.dumps(summary_payload).encode(), b""
            if self.command[0] == "unzip":
                extract_dir = tmp_path / "Mycobacterium_tuberculosis" / self.command[3]
                data_dir = extract_dir / "ncbi_dataset" / "data" / "GCA_HIGH"
                data_dir.mkdir(parents=True)
                (data_dir / "GCA_HIGH.fna").write_text(">chr1\nACGT\n")
                return b"", b""
            return b"", b""

    with patch("Afanc.get_dataset.download_assemblies.command", FakeCommand):
        selected = download_genome("Mycobacterium tuberculosis", _args(num_assemblies=2))

    assert selected == ["GCA_HIGH", "GCA_MID"]
    assert (tmp_path / "Mycobacterium_tuberculosis" / "GCA_HIGH.fna").exists()
    download_call = [call for call in command_calls if call[1] == "DOWNLOAD"][0][0]
    assert download_call[:5] == ["datasets", "download", "genome", "accession", "GCA_HIGH"]
    assert "GCA_MID" in download_call
    assert "GCA_LOW" not in download_call


def test_run_get_dataset_skips_invalid_ids_and_processes_valid_names(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    names = tmp_path / "names.txt"
    names.write_text("Good species\nBad/species\n\n")
    processed = []

    def fake_download_genome(assembly, _args):
        processed.append(assembly)

    args = _args(ID_file=names, output_prefix="out")
    with patch("Afanc.get_dataset.download_assemblies.download_genome", fake_download_genome):
        runGet_dataset(args)

    assert processed == ["Good species"]
    assert "Found 2 IDs" in args.stdout.getvalue()
    assert "Invalid characters found in Bad/species" in args.stdout.getvalue()
    assert (tmp_path / "out").is_dir()


def _report(accession, n50):
    return {
        "accession": accession,
        "organism": {"organism_name": accession},
        "average_nucleotide_identity": {
            "best_ani_match": {
                "ani": 99.0,
                "organism_name": accession,
            }
        },
        "assembly_info": {"assembly_level": "Complete Genome"},
        "assembly_stats": {"scaffold_n50": n50},
    }

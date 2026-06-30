from pathlib import Path
from unittest.mock import patch

import pytest

from Afanc.screen.mapping.bwa import (
    build_mapping_pipeline_command,
    build_read_group,
    map_reads_to_bam,
    prepare_reference,
)
from Afanc.screen.maths.mappingMetrics import breadthofCoverage
from Afanc.screen.mapHits import (
    build_combined_reference,
    make_accessions_dict,
    make_low_genome_coverage_warning,
)


def test_build_read_group():
    assert build_read_group("SAMPLE") == "@RG\\tID:SAMPLE\\tSM:SAMPLE"


@patch("Afanc.screen.mapping.bwa.command")
@patch("Afanc.screen.mapping.bwa.shutil.which", return_value="/usr/bin/tool")
def test_prepare_reference_runs_missing_indexes(mock_which, mock_command_class, tmp_path):
    ref_fasta = tmp_path / "ref.fa"
    ref_fasta.write_text(">chr1\nACGT\n")

    mock_command_instance = mock_command_class.return_value

    result = prepare_reference(ref_fasta)

    assert result["reference_fasta"] == ref_fasta
    assert result["fai_path"] == Path(str(ref_fasta) + ".fai")
    assert mock_command_class.call_count == 2
    assert mock_command_class.call_args_list[0].args == (["samtools", "faidx", str(ref_fasta)], "MAP")
    assert mock_command_class.call_args_list[1].args == (["bwa", "index", str(ref_fasta)], "MAP")
    assert mock_command_instance.run_comm_quiet.call_count == 2


@patch("Afanc.screen.mapping.bwa.command")
@patch("Afanc.screen.mapping.bwa.shutil.which", return_value="/usr/bin/tool")
def test_prepare_reference_skips_existing_indexes(mock_which, mock_command_class, tmp_path):
    ref_fasta = tmp_path / "ref.fa"
    ref_fasta.write_text(">chr1\nACGT\n")
    Path(str(ref_fasta) + ".fai").write_text("chr1\t4\t6\t4\t5\n")
    for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
        Path(str(ref_fasta) + suffix).write_text("")

    prepare_reference(ref_fasta)

    mock_command_class.assert_not_called()


def test_build_mapping_pipeline_command(tmp_path):
    ref_fasta = tmp_path / "ref.fa"
    r1_fastq = tmp_path / "reads_1.fq.gz"
    r2_fastq = tmp_path / "reads_2.fq.gz"
    output_bam = tmp_path / "sample.bam"

    command = build_mapping_pipeline_command(
        ref_fasta=ref_fasta,
        r1_fastq=r1_fastq,
        r2_fastq=r2_fastq,
        output_bam=output_bam,
        sample_name="SAMPLE",
        cpus=8,
        sort_threads=4,
        sort_mem="2G",
        tmpdir=tmp_path,
    )

    assert "bwa mem -Y -M" in command
    assert "-R '@RG\\tID:SAMPLE\\tSM:SAMPLE'" in command
    assert f"{ref_fasta}.fai" in command
    assert "samclip --max 10" in command
    assert "samtools fixmate -m" in command
    assert "samtools markdup -r -s" in command
    assert f"> {output_bam}" in command


@patch("Afanc.screen.mapping.bwa.command")
@patch("Afanc.screen.mapping.bwa.prepare_reference")
@patch("Afanc.screen.mapping.bwa._require_executable")
def test_map_reads_to_bam_runs_pipeline_and_index(
    mock_require_executable,
    mock_prepare_reference,
    mock_command_class,
    tmp_path,
):
    ref_fasta = tmp_path / "ref.fa"
    r1_fastq = tmp_path / "reads_1.fq.gz"
    r2_fastq = tmp_path / "reads_2.fq.gz"
    output_bam = tmp_path / "sample.bam"

    ref_fasta.write_text(">chr1\nACGT\n")
    r1_fastq.write_text("")
    r2_fastq.write_text("")

    mock_prepare_reference.return_value = {
        "reference_fasta": ref_fasta,
        "fai_path": Path(str(ref_fasta) + ".fai"),
        "bwa_index_paths": {},
    }
    mock_command_instance = mock_command_class.return_value

    result = map_reads_to_bam(
        ref_fasta=ref_fasta,
        r1_fastq=r1_fastq,
        r2_fastq=r2_fastq,
        output_bam=output_bam,
        sample_name="SAMPLE",
        cpus=8,
        sort_threads=4,
        sort_mem="2G",
        tmpdir=tmp_path,
    )

    mock_prepare_reference.assert_called_once()
    mock_require_executable.assert_called_once_with("samclip")
    assert mock_command_class.call_count == 2
    assert mock_command_class.call_args_list[0].args[1] == "MAP"
    assert mock_command_class.call_args_list[0].kwargs == {
        "shell": True,
        "pipefail": True,
    }
    assert mock_command_class.call_args_list[1].args == (["samtools", "index", str(output_bam)], "MAP")
    assert mock_command_instance.run_comm_quiet.call_count == 2

    assert result["bam_path"] == output_bam
    assert result["bai_path"] == Path(str(output_bam) + ".bai")


def test_map_reads_to_bam_requires_existing_inputs(tmp_path):
    ref_fasta = tmp_path / "missing_ref.fa"
    r1_fastq = tmp_path / "reads_1.fq.gz"
    r2_fastq = tmp_path / "reads_2.fq.gz"

    with pytest.raises(Exception):
        map_reads_to_bam(
            ref_fasta=ref_fasta,
            r1_fastq=r1_fastq,
            r2_fastq=r2_fastq,
            output_bam=tmp_path / "sample.bam",
            sample_name="SAMPLE",
        )


def test_build_combined_reference_marks_contigs_by_accession(tmp_path):
    assembly_a = tmp_path / "ACC_A.fna"
    assembly_b = tmp_path / "ACC_B.fna"
    combined = tmp_path / "combined.fa"

    assembly_a.write_text(">chr1 first contig\nACGT\n")
    assembly_b.write_text(">contig2 second contig\nTGCA\n")

    result = build_combined_reference(
        {
            "ACC_A": str(assembly_a),
            "ACC_B": str(assembly_b),
        },
        output_fasta=str(combined),
    )

    combined_text = combined.read_text()

    assert result == str(combined)
    assert ">chr1___ACC_A" in combined_text
    assert ">contig2___ACC_B" in combined_text


def test_make_accessions_dict_includes_fasta_extension(tmp_path):
    class Args:
        mappingWDir = tmp_path

    reference = tmp_path / "H37Rv.fasta"
    reference.write_text(">chr1\nACGT\n")

    assemblies = make_accessions_dict(Args)

    assert assemblies["H37Rv"] == "H37Rv.fasta"


def test_low_genome_coverage_warning_is_not_created_at_threshold():
    warning = make_low_genome_coverage_warning(0.05, "reference.fna")

    assert warning is None


def test_breadth_of_coverage_ignores_trailing_blank_depth_row():
    coverage_rows = [
        ["chr1", "1", "5"],
        ["chr1", "2", "3"],
        [""],
    ]

    assert breadthofCoverage(coverage_rows, genomesize=4) == 0.5


@pytest.mark.parametrize(
    ("coverage_fraction", "expected_severity"),
    [
        (0.049, "low"),
        (0.009, "very_low"),
    ],
)
def test_low_genome_coverage_warning_marks_possible_spurious_result(
    coverage_fraction,
    expected_severity,
):
    warning = make_low_genome_coverage_warning(coverage_fraction, "reference.fna")

    assert warning["code"] == "low_genome_coverage"
    assert warning["flag"] is True
    assert warning["severity"] == expected_severity
    assert warning["genome_coverage_fraction"] == coverage_fraction
    assert warning["warning_threshold"] == 0.05
    assert warning["possible_spurious_result"] is True
    assert "could be a spurious result" in warning["message"]

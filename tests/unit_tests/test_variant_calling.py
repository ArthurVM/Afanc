from pathlib import Path
from unittest.mock import patch

import pytest

from Afanc.screen.variant_calling.callers import (
    build_bcftools_call_command,
    build_bcftools_mpileup_command,
    build_freebayes_call_command,
    build_freebayes_filter_command,
    build_freebayes_regions_command,
    run_variant_caller,
)
from Afanc.screen.variant_calling.profile import run_snp_profiling, vcf_to_snp_json


def test_build_freebayes_regions_command(tmp_path):
    ref_fasta = tmp_path / "ref.fa"
    out_regions = tmp_path / "sample.regions.txt"

    command = build_freebayes_regions_command(ref_fasta, out_regions, 50000)

    assert "fasta_generate_regions.py" in command
    assert f"{ref_fasta}.fai" in command
    assert "50000" in command
    assert str(out_regions) in command


def test_build_freebayes_call_command(tmp_path):
    bam_path = tmp_path / "sample.bam"
    ref_fasta = tmp_path / "ref.fa"
    regions_file = tmp_path / "sample.regions.txt"
    raw_vcf = tmp_path / "sample.raw.vcf"

    command = build_freebayes_call_command(
        bam_path=bam_path,
        reference_fasta=ref_fasta,
        regions_file=regions_file,
        output_vcf=raw_vcf,
        cpus=8,
    )

    assert "freebayes-parallel" in command
    assert f"{regions_file} 8" in command
    assert "--strict-vcf" in command
    assert str(ref_fasta) in command
    assert str(bam_path) in command
    assert str(raw_vcf) in command


def test_build_freebayes_filter_command(tmp_path):
    raw_vcf = tmp_path / "sample.raw.vcf"
    ref_fasta = tmp_path / "ref.fa"
    filt_vcf = tmp_path / "sample.filt.vcf"

    command = build_freebayes_filter_command(raw_vcf, ref_fasta, filt_vcf)

    assert "bcftools view" in command
    assert 'FMT/GT="1/1"' in command
    assert "FMT/AO" in command
    assert "bcftools norm" in command
    assert str(filt_vcf) in command


def test_build_bcftools_commands(tmp_path):
    bam_path = tmp_path / "sample.bam"
    ref_fasta = tmp_path / "ref.fa"
    raw_vcf = tmp_path / "sample.raw.vcf"

    mpileup_command = build_bcftools_mpileup_command(bam_path, ref_fasta)
    call_command = build_bcftools_call_command(raw_vcf)

    assert "bcftools mpileup" in mpileup_command
    assert "FORMAT/AD,FORMAT/DP" in mpileup_command
    assert str(bam_path) in mpileup_command
    assert "bcftools call" in call_command
    assert str(raw_vcf) in call_command


@patch("Afanc.screen.variant_calling.callers.run_freebayes_variant_caller")
def test_run_variant_caller_defaults_to_freebayes(mock_run_freebayes_variant_caller, tmp_path):
    bam_path = tmp_path / "sample.bam"
    ref_fasta = tmp_path / "ref.fa"
    output_prefix = tmp_path / "sample"

    run_variant_caller(
        bam_path=bam_path,
        reference_fasta=ref_fasta,
        output_prefix=output_prefix,
    )

    mock_run_freebayes_variant_caller.assert_called_once()


def test_run_variant_caller_rejects_unknown_backend(tmp_path):
    bam_path = tmp_path / "sample.bam"
    ref_fasta = tmp_path / "ref.fa"
    output_prefix = tmp_path / "sample"

    with pytest.raises(ValueError):
        run_variant_caller(
            bam_path=bam_path,
            reference_fasta=ref_fasta,
            output_prefix=output_prefix,
            caller_name="unknown",
        )


def test_vcf_to_snp_json_writes_classifier_ready_json(tmp_path):
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                "chr1\t42\t.\tA\tG\t99\tPASS\tDP=12\tGT:DP:AD\t1/1:12:1,11",
                "chr1\t44\t.\tC\tT\t20\tPASS\tDP=12\tGT:DP:AD\t1/1:12:1,11",
                "chr1\t50\t.\tAT\tA\t99\tPASS\t.\tGT:DP:AD\t1/1:12:1,11",
                "",
            ]
        )
    )
    depth_path = tmp_path / "sample.depth.bed"
    depth_path.write_text("chr1\t41\t10\nchr1\t42\t12\nchr1\t43\t0\n")
    output_json = tmp_path / "sample.snps.json"

    snp_json = vcf_to_snp_json(
        vcf_path=vcf_path,
        depth_path=depth_path,
        output_json=output_json,
        min_qual=30,
        min_dp=10,
        min_missing_depth=10,
    )

    assert output_json.exists()
    assert set(snp_json) == {"alleles", "missing"}
    assert snp_json["alleles"] == ["chr1.42.A.G"]
    assert snp_json["missing"] == [["chr1", "43"]]


@patch("Afanc.screen.variant_calling.profile.prepare_variant_reference")
@patch("Afanc.screen.variant_calling.profile.vcf_to_snp_json")
@patch("Afanc.screen.variant_calling.profile.write_depth_bed")
@patch("Afanc.screen.variant_calling.profile.run_variant_caller")
def test_run_snp_profiling_passes_cpus_only_to_freebayes(
    mock_run_variant_caller,
    mock_write_depth_bed,
    mock_vcf_to_snp_json,
    mock_prepare_variant_reference,
    tmp_path,
):
    class Args:
        profilerWDir = tmp_path
        reportsDir = tmp_path
        output_prefix = "sample"
        variant_caller = "freebayes"
        threads = 8
        stdout = None
        stderr = None
        snp_min_qual = 30.0
        snp_allow_filtered = False
        snp_accept_missing_qual = False
        snp_min_dp = None
        snp_min_missing_depth = 10

    mock_prepare_variant_reference.return_value = tmp_path / "ref.fa"
    mock_run_variant_caller.return_value = {
        "raw_vcf": tmp_path / "ACC.raw.vcf",
        "filtered_vcf": tmp_path / "ACC.filt.vcf",
    }
    mock_vcf_to_snp_json.return_value = {"alleles": ["chr1.42.A.G"], "missing": [["chr1", "43"]]}

    run_snp_profiling(
        Args,
        {"ACC": {"bam": tmp_path / "ACC.bam", "assembly": tmp_path / "ACC.fa"}},
    )

    mock_run_variant_caller.assert_called_once()
    assert mock_run_variant_caller.call_args.kwargs["caller_name"] == "freebayes"
    assert mock_run_variant_caller.call_args.kwargs["cpus"] == 8


@patch("Afanc.screen.variant_calling.profile.prepare_variant_reference")
@patch("Afanc.screen.variant_calling.profile.vcf_to_snp_json")
@patch("Afanc.screen.variant_calling.profile.write_depth_bed")
@patch("Afanc.screen.variant_calling.profile.run_variant_caller")
def test_run_snp_profiling_does_not_pass_cpus_to_bcftools(
    mock_run_variant_caller,
    mock_write_depth_bed,
    mock_vcf_to_snp_json,
    mock_prepare_variant_reference,
    tmp_path,
):
    class Args:
        profilerWDir = tmp_path
        reportsDir = tmp_path
        output_prefix = "sample"
        variant_caller = "bcftools"
        threads = 8
        stdout = None
        stderr = None
        snp_min_qual = 30.0
        snp_allow_filtered = False
        snp_accept_missing_qual = False
        snp_min_dp = None
        snp_min_missing_depth = 10

    mock_prepare_variant_reference.return_value = tmp_path / "ref.fa"
    mock_run_variant_caller.return_value = {
        "raw_vcf": tmp_path / "ACC.raw.vcf",
        "filtered_vcf": None,
    }
    mock_vcf_to_snp_json.return_value = {"alleles": ["chr1.42.A.G"], "missing": [["chr1", "43"]]}

    run_snp_profiling(
        Args,
        {"ACC": {"bam": tmp_path / "ACC.bam", "assembly": tmp_path / "ACC.fa"}},
    )

    mock_run_variant_caller.assert_called_once()
    assert mock_run_variant_caller.call_args.kwargs["caller_name"] == "bcftools"
    assert "cpus" not in mock_run_variant_caller.call_args.kwargs

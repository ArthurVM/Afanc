from __future__ import annotations

import shlex
import shutil
from pathlib import Path
from typing import Any, Dict, Optional, Union

from Afanc.utilities.exceptions import FileNotFoundErrorAfanc
from Afanc.utilities.runCommands import command


PathLike = Union[str, Path]


def build_freebayes_regions_command(
    reference_fasta: PathLike,
    output_regions: PathLike,
    chunk_size: int,
    fasta_generate_regions_executable: str = "fasta_generate_regions.py",
) -> str:
    reference_fasta = Path(reference_fasta)
    output_regions = Path(output_regions)

    return (
        f"{shlex.quote(fasta_generate_regions_executable)} "
        f"{shlex.quote(str(reference_fasta) + '.fai')} "
        f"{int(chunk_size)} > {shlex.quote(str(output_regions))}"
    )


def build_freebayes_call_command(
    bam_path: PathLike,
    reference_fasta: PathLike,
    regions_file: PathLike,
    output_vcf: PathLike,
    cpus: int = 1,
    freebayes_parallel_executable: str = "freebayes-parallel",
    ploidy: int = 2,
    min_alternate_fraction: float = 0.05,
    min_alternate_count: int = 2,
    min_coverage: int = 10,
    min_repeat_entropy: float = 1.0,
    min_base_quality: int = 13,
    min_mapping_quality: int = 60,
) -> str:
    bam_path = Path(bam_path)
    reference_fasta = Path(reference_fasta)
    regions_file = Path(regions_file)
    output_vcf = Path(output_vcf)

    ## retain permissive calls as classification evidence
    return (
        f"{shlex.quote(freebayes_parallel_executable)} "
        f"{shlex.quote(str(regions_file))} {int(cpus)} "
        f"-p {int(ploidy)} "
        f"-P 0 "
        f"-C {int(min_alternate_count)} "
        f"-F {min_alternate_fraction} "
        f"--min-coverage {int(min_coverage)} "
        f"--min-repeat-entropy {min_repeat_entropy} "
        f"-q {int(min_base_quality)} "
        f"-m {int(min_mapping_quality)} "
        f"--strict-vcf "
        f"-f {shlex.quote(str(reference_fasta))} "
        f"{shlex.quote(str(bam_path))} > {shlex.quote(str(output_vcf))}"
    )


def build_freebayes_filter_command(
    raw_vcf: PathLike,
    reference_fasta: PathLike,
    output_vcf: PathLike,
    bcftools_executable: str = "bcftools",
    min_quality: int = 100,
    min_depth: int = 10,
    min_alt_fraction: float = 0.9,
    require_hom_alt: bool = True,
) -> str:
    raw_vcf = Path(raw_vcf)
    reference_fasta = Path(reference_fasta)
    output_vcf = Path(output_vcf)

    gt_filter = 'FMT/GT="1/1" && ' if require_hom_alt else ""
    include_expr = (
        f"{gt_filter}"
        f"QUAL>={int(min_quality)} && "
        f"FMT/DP>={int(min_depth)} && "
        f"(FMT/AO)/(FMT/DP)>={min_alt_fraction}"
    )

    ## retain the raw VCF for evidence extraction
    return (
        f"{shlex.quote(bcftools_executable)} view --include {shlex.quote(include_expr)} "
        f"{shlex.quote(str(raw_vcf))} "
        f"| {shlex.quote(bcftools_executable)} norm -f {shlex.quote(str(reference_fasta))} - "
        f"> {shlex.quote(str(output_vcf))}"
    )


def build_bcftools_mpileup_command(
    bam_path: PathLike,
    reference_fasta: PathLike,
    targets_bed: Optional[PathLike] = None,
    bcftools_executable: str = "bcftools",
) -> str:
    bam_path = Path(bam_path)
    reference_fasta = Path(reference_fasta)

    command_box = [
        shlex.quote(bcftools_executable),
        "mpileup",
        "-f",
        shlex.quote(str(reference_fasta)),
        "-a",
        "FORMAT/AD,FORMAT/DP",
        "-Ou",
    ]
    if targets_bed is not None:
        command_box.extend(["-T", shlex.quote(str(targets_bed))])
    command_box.append(shlex.quote(str(bam_path)))

    return " ".join(command_box)


def build_bcftools_call_command(
    output_vcf: PathLike,
    bcftools_executable: str = "bcftools",
) -> str:
    output_vcf = Path(output_vcf)
    return (
        f"{shlex.quote(bcftools_executable)} call -m -Ov "
        f"-o {shlex.quote(str(output_vcf))}"
    )


def run_freebayes_variant_caller(
    bam_path: PathLike,
    reference_fasta: PathLike,
    output_prefix: PathLike,
    cpus: int = 1,
    chunk_size: int = 100000,
    freebayes_parallel_executable: str = "freebayes-parallel",
    fasta_generate_regions_executable: str = "fasta_generate_regions.py",
    bcftools_executable: str = "bcftools",
    write_filtered_vcf: bool = True,
    **freebayes_kwargs: Any,
) -> Dict[str, Path]:
    bam_path = Path(bam_path)
    reference_fasta = Path(reference_fasta)
    output_prefix = Path(output_prefix)

    _check_variant_calling_inputs(bam_path, reference_fasta)
    _require_executable(freebayes_parallel_executable)
    _require_executable(fasta_generate_regions_executable)
    _require_executable(bcftools_executable)

    regions_file = output_prefix.with_suffix(".regions.txt")
    raw_vcf = output_prefix.with_suffix(".raw.vcf")
    filtered_vcf = output_prefix.with_suffix(".filt.vcf")

    command(
        build_freebayes_regions_command(
            reference_fasta=reference_fasta,
            output_regions=regions_file,
            chunk_size=chunk_size,
            fasta_generate_regions_executable=fasta_generate_regions_executable,
        ),
        "VARCALL",
        shell=True,
    ).run_comm_quiet(0)

    command(
        build_freebayes_call_command(
            bam_path=bam_path,
            reference_fasta=reference_fasta,
            regions_file=regions_file,
            output_vcf=raw_vcf,
            cpus=cpus,
            freebayes_parallel_executable=freebayes_parallel_executable,
            **freebayes_kwargs,
        ),
        "VARCALL",
        shell=True,
    ).run_comm_quiet(0)

    if write_filtered_vcf:
        command(
            build_freebayes_filter_command(
                raw_vcf=raw_vcf,
                reference_fasta=reference_fasta,
                output_vcf=filtered_vcf,
                bcftools_executable=bcftools_executable,
            ),
            "VARCALL",
            shell=True,
            pipefail=True,
        ).run_comm_quiet(0)
    else:
        filtered_vcf = None

    return {
        "caller": Path(freebayes_parallel_executable),
        "regions_file": regions_file,
        "raw_vcf": raw_vcf,
        "filtered_vcf": filtered_vcf,
    }


def run_bcftools_variant_caller(
    bam_path: PathLike,
    reference_fasta: PathLike,
    output_prefix: PathLike,
    targets_bed: Optional[PathLike] = None,
    bcftools_executable: str = "bcftools",
) -> Dict[str, Optional[Path]]:
    bam_path = Path(bam_path)
    reference_fasta = Path(reference_fasta)
    output_prefix = Path(output_prefix)

    _check_variant_calling_inputs(bam_path, reference_fasta)
    _require_executable(bcftools_executable)

    raw_vcf = output_prefix.with_suffix(".raw.vcf")
    mpileup_command = build_bcftools_mpileup_command(
        bam_path=bam_path,
        reference_fasta=reference_fasta,
        targets_bed=targets_bed,
        bcftools_executable=bcftools_executable,
    )
    call_command = build_bcftools_call_command(
        output_vcf=raw_vcf,
        bcftools_executable=bcftools_executable,
    )

    ## alternate caller for compatibility and benchmarking
    command(
        f"{mpileup_command} | {call_command}",
        "VARCALL",
        shell=True,
        pipefail=True,
    ).run_comm_quiet(0)

    return {
        "caller": Path(bcftools_executable),
        "regions_file": None,
        "raw_vcf": raw_vcf,
        "filtered_vcf": None,
    }


def run_variant_caller(
    bam_path: PathLike,
    reference_fasta: PathLike,
    output_prefix: PathLike,
    caller_name: str = "freebayes",
    **caller_kwargs: Any,
) -> Dict[str, Optional[Path]]:
    if caller_name == "freebayes":
        return run_freebayes_variant_caller(
            bam_path=bam_path,
            reference_fasta=reference_fasta,
            output_prefix=output_prefix,
            **caller_kwargs,
        )

    if caller_name == "bcftools":
        return run_bcftools_variant_caller(
            bam_path=bam_path,
            reference_fasta=reference_fasta,
            output_prefix=output_prefix,
            **caller_kwargs,
        )

    raise ValueError(f"Unsupported caller_name={caller_name!r}. Expected 'freebayes' or 'bcftools'.")


def _check_variant_calling_inputs(bam_path: Path, reference_fasta: Path) -> None:
    for required_path in [bam_path, reference_fasta]:
        if not required_path.is_file():
            raise FileNotFoundErrorAfanc(required_path)

    fai_path = Path(str(reference_fasta) + ".fai")
    if not fai_path.is_file():
        raise FileNotFoundErrorAfanc(
            fai_path,
            message=f"Reference FASTA index not found: {fai_path}",
        )


def _require_executable(executable_name: str) -> None:
    if shutil.which(executable_name) is None:
        raise FileNotFoundErrorAfanc(
            executable_name,
            message=f"Required executable not found on PATH: {executable_name}",
        )

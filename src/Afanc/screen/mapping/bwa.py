from __future__ import annotations

import shlex
import shutil
from pathlib import Path
from typing import Dict, Optional, Union

from Afanc.utilities.exceptions import FileNotFoundErrorAfanc
from Afanc.utilities.runCommands import command


PathLike = Union[str, Path]


def build_read_group(sample_name: str) -> str:
    """Build a minimal read group string for bwa mem."""

    return f"@RG\\tID:{sample_name}\\tSM:{sample_name}"


def prepare_reference(
    ref_fasta: PathLike,
    bwa_executable: str = "bwa",
    samtools_executable: str = "samtools",
) -> Dict[str, object]:
    """Prepare the reference FASTA for bwa/samtools mapping.

    Returns a small metadata dictionary describing the prepared assets.
    """

    ref_fasta = Path(ref_fasta)
    if not ref_fasta.is_file():
        raise FileNotFoundErrorAfanc(ref_fasta)

    _require_executable(bwa_executable)
    _require_executable(samtools_executable)

    fai_path = ref_fasta.with_suffix(ref_fasta.suffix + ".fai")
    bwa_index_paths = _bwa_index_paths(ref_fasta)

    if not fai_path.exists():
        command([samtools_executable, "faidx", str(ref_fasta)], "MAP").run_comm_quiet(0)

    if not all(index_path.exists() for index_path in bwa_index_paths.values()):
        command([bwa_executable, "index", str(ref_fasta)], "MAP").run_comm_quiet(0)

    return {
        "reference_fasta": ref_fasta,
        "fai_path": fai_path,
        "bwa_index_paths": bwa_index_paths,
    }


def build_mapping_pipeline_command(
    ref_fasta: PathLike,
    r1_fastq: PathLike,
    r2_fastq: PathLike,
    output_bam: PathLike,
    sample_name: str,
    cpus: int = 1,
    sort_threads: Optional[int] = None,
    sort_mem: str = "1G",
    tmpdir: PathLike = "/tmp",
    bwa_executable: str = "bwa",
    samtools_executable: str = "samtools",
    samclip_executable: str = "samclip",
) -> str:
    """Construct the bwa/samclip/samtools mapping pipeline."""

    ref_fasta = Path(ref_fasta)
    r1_fastq = Path(r1_fastq)
    r2_fastq = Path(r2_fastq)
    output_bam = Path(output_bam)
    tmpdir = Path(tmpdir)
    sort_threads = sort_threads if sort_threads is not None else cpus
    read_group = build_read_group(sample_name)

    pipeline = (
        f"{shlex.quote(bwa_executable)} mem -Y -M "
        f"-R {shlex.quote(read_group)} "
        f"-t {int(cpus)} "
        f"{shlex.quote(str(ref_fasta))} "
        f"{shlex.quote(str(r1_fastq))} "
        f"{shlex.quote(str(r2_fastq))} "
        f"| {shlex.quote(samclip_executable)} --max 10 --ref {shlex.quote(str(ref_fasta) + '.fai')} "
        f"| {shlex.quote(samtools_executable)} sort -n -l 0 "
        f"-T {shlex.quote(str(tmpdir))} "
        f"--threads {int(sort_threads)} "
        f"-m {shlex.quote(str(sort_mem))} "
        f"| {shlex.quote(samtools_executable)} fixmate -m "
        f"--threads {int(sort_threads)} - - "
        f"| {shlex.quote(samtools_executable)} sort -l 0 "
        f"-T {shlex.quote(str(tmpdir))} "
        f"--threads {int(sort_threads)} "
        f"-m {shlex.quote(str(sort_mem))} "
        f"| {shlex.quote(samtools_executable)} markdup -r -s "
        f"-T {shlex.quote(str(tmpdir))} "
        f"--threads {int(sort_threads)} - - "
        f"> {shlex.quote(str(output_bam))}"
    )

    return pipeline


def map_reads_to_bam(
    ref_fasta: PathLike,
    r1_fastq: PathLike,
    r2_fastq: PathLike,
    output_bam: PathLike,
    sample_name: str,
    cpus: int = 1,
    sort_threads: Optional[int] = None,
    sort_mem: str = "1G",
    tmpdir: PathLike = "/tmp",
    bwa_executable: str = "bwa",
    samtools_executable: str = "samtools",
    samclip_executable: str = "samclip",
) -> Dict[str, Path]:
    """Map paired reads to a reference FASTA and produce an indexed BAM."""

    ref_fasta = Path(ref_fasta)
    r1_fastq = Path(r1_fastq)
    r2_fastq = Path(r2_fastq)
    output_bam = Path(output_bam)

    for required_path in [ref_fasta, r1_fastq, r2_fastq]:
        if not required_path.is_file():
            raise FileNotFoundErrorAfanc(required_path)

    prepared_reference = prepare_reference(
        ref_fasta=ref_fasta,
        bwa_executable=bwa_executable,
        samtools_executable=samtools_executable,
    )

    _require_executable(samclip_executable)

    output_bam.parent.mkdir(parents=True, exist_ok=True)

    pipeline = build_mapping_pipeline_command(
        ref_fasta=ref_fasta,
        r1_fastq=r1_fastq,
        r2_fastq=r2_fastq,
        output_bam=output_bam,
        sample_name=sample_name,
        cpus=cpus,
        sort_threads=sort_threads,
        sort_mem=sort_mem,
        tmpdir=tmpdir,
        bwa_executable=bwa_executable,
        samtools_executable=samtools_executable,
        samclip_executable=samclip_executable,
    )

    command(
        pipeline,
        "MAP",
        shell=True,
        pipefail=True,
    ).run_comm_quiet(0)

    command([samtools_executable, "index", str(output_bam)], "MAP").run_comm_quiet(0)

    return {
        "reference_fasta": ref_fasta,
        "fai_path": prepared_reference["fai_path"],
        "bam_path": output_bam,
        "bai_path": Path(str(output_bam) + ".bai"),
    }


def _bwa_index_paths(ref_fasta: Path) -> Dict[str, Path]:
    return {
        "amb": Path(str(ref_fasta) + ".amb"),
        "ann": Path(str(ref_fasta) + ".ann"),
        "bwt": Path(str(ref_fasta) + ".bwt"),
        "pac": Path(str(ref_fasta) + ".pac"),
        "sa": Path(str(ref_fasta) + ".sa"),
    }


def _require_executable(executable_name: str) -> None:
    if shutil.which(executable_name) is None:
        raise FileNotFoundErrorAfanc(
            executable_name,
            message=f"Required executable not found on PATH: {executable_name}",
        )

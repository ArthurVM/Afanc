#!/usr/bin/env python3
"""Regenerate the deterministic synthetic Afanc reviewer dataset."""

from pathlib import Path
import random


SEED = 20260629
GENOME_LENGTH = 24_000
READ_LENGTH = 150
FRAGMENT_LENGTH = 350
READ_PAIRS = 250

ROOT = Path(__file__).resolve().parent
ASSEMBLY = (
    ROOT
    / "assemblies"
    / "Toybacter_alpha"
    / "GCF_900000001.1_toy_genomic.fna"
)
READ_1 = ROOT / "reads" / "toy_R1.fastq"
READ_2 = ROOT / "reads" / "toy_R2.fastq"


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def wrap(sequence, width=80):
    return "\n".join(
        sequence[offset : offset + width]
        for offset in range(0, len(sequence), width)
    )


def main():
    rng = random.Random(SEED)
    genome = "".join(rng.choices("ACGT", weights=(27, 23, 23, 27), k=GENOME_LENGTH))

    ASSEMBLY.parent.mkdir(parents=True, exist_ok=True)
    READ_1.parent.mkdir(parents=True, exist_ok=True)
    ASSEMBLY.write_text(f">toy_chromosome\n{wrap(genome)}\n")

    max_start = GENOME_LENGTH - FRAGMENT_LENGTH
    starts = [
        round(index * max_start / (READ_PAIRS - 1))
        for index in range(READ_PAIRS)
    ]
    quality = "I" * READ_LENGTH

    with READ_1.open("w") as read_1, READ_2.open("w") as read_2:
        for pair_number, start in enumerate(starts, start=1):
            fragment = genome[start : start + FRAGMENT_LENGTH]
            r1_sequence = fragment[:READ_LENGTH]
            r2_sequence = reverse_complement(fragment[-READ_LENGTH:])
            read_name = f"toy_{pair_number:04d}_start_{start}"
            read_1.write(f"@{read_name}/1\n{r1_sequence}\n+\n{quality}\n")
            read_2.write(f"@{read_name}/2\n{r2_sequence}\n+\n{quality}\n")

    print(f"Wrote {ASSEMBLY}")
    print(f"Wrote {READ_1}")
    print(f"Wrote {READ_2}")


if __name__ == "__main__":
    main()

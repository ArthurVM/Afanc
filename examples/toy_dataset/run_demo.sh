#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_dir="${1:-$PWD/afanc-toy-run}"

required_executables=(
    afanc
    mash
    kraken2
    kraken2-build
    kraken2-inspect
    ktImportText
    bwa
    samtools
    samclip
    freebayes-parallel
    fasta_generate_regions.py
    bcftools
)

for executable in "${required_executables[@]}"; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "Required executable not found on PATH: $executable" >&2
        exit 1
    fi
done

if [[ -e "$output_dir" ]]; then
    echo "Output path already exists; choose a new path: $output_dir" >&2
    exit 1
fi

mkdir -p "$output_dir"
output_dir="$(cd "$output_dir" && pwd)"

echo ""
echo "=========== Building toy Afanc database in $output_dir/toy_db ==========="
echo ""
(
    cd "$output_dir"
    afanc autodatabase \
        "$script_dir/assemblies" \
        --ncbi_tax_db "$script_dir/taxonomy" \
        --output_prefix toy_db \
        --threads 2
)

echo ""
echo ""
echo "=========== Screening toy paired-end reads ==========="
echo ""
(
    cd "$output_dir"
    afanc screen \
        toy_db \
        "$script_dir/reads/toy_R1.fastq" \
        "$script_dir/reads/toy_R2.fastq" \
        --output_prefix toy_screen \
        --threads 2 \
        --pct_threshold 5 \
        --num_threshold 20
)

result="$output_dir/toy_screen/toy_screen.json"
echo ""
echo "=========== Toy demonstration complete: $result ==========="
echo ""


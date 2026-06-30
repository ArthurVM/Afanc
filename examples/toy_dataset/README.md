# Afanc toy dataset

This directory contains a deterministic, entirely synthetic dataset for a quick offline demonstration of:

1. `afanc autodatabase`
2. `afanc screen`

The data describe one fictional species, `Toybacter alpha`, with a 24 kbp assembly and 250 paired-end read pairs.

## One-command demonstration

Activate an environment containing Afanc and its command-line dependencies, then run from anywhere:

```bash
examples/toy_dataset/run_demo.sh
```

By default, outputs are written to `./afanc-toy-run`. An alternative new output path can be supplied as the first argument:

```bash
examples/toy_dataset/run_demo.sh /tmp/afanc-toy-run
```

The script builds `toy_db`, runs the full screen (including mapping and SNP profiling), and writes the main result to:

```text
<output-directory>/toy_screen/toy_screen.json
```

The expected detection event is taxID `1000001`, `Toybacter alpha`.

## Equivalent commands

```bash
mkdir afanc-toy-run
cd afanc-toy-run

afanc autodatabase \
  /path/to/Afanc/examples/toy_dataset/assemblies \
  --ncbi_tax_db /path/to/Afanc/examples/toy_dataset/taxonomy \
  --output_prefix toy_db \
  --threads 2

afanc screen \
  toy_db \
  /path/to/Afanc/examples/toy_dataset/reads/toy_R1.fastq \
  /path/to/Afanc/examples/toy_dataset/reads/toy_R2.fastq \
  --output_prefix toy_screen \
  --threads 2 \
  --pct_threshold 5 \
  --num_threshold 20
```

For an even shorter taxonomy-only smoke test, add `--no_map`. That still runs Kraken2, deconvolution, JSON reporting, and Krona reporting, but skips mapping and variant calling.

## Rebuilding the inputs

The committed FASTA and FASTQ files can be reproduced exactly with:

```bash
python examples/toy_dataset/generate_toy_data.py
```

The bundled minimal NCBI-format taxonomy avoids any taxonomy download. It contains only the lineages required by this fictional organism.
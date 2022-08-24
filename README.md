# Afanc: High Resolution Metagenomic Disambiguator
A tool which can perform variant level metageonomics disambiguation of NGS reads.

## Installation
Installation is quick and simple, simply run

  git clone git@github.com:ArthurVM/Afanc.git
  cd Afanc
  pip3 install ./

## Sub-Modules
Afanc is split into 3 sub-modules:

  get_dataset         Download a dataset of genome assemblies from GenBank.
  autodatabase        Generate a database from a FASTA directory structure with autodatabase.
  screen              High-resolution metagenomic screening of short read data using a database
                      constructed by autodatabase.

These modules enable the user to download a dataset of genome assemblies using species names or accession IDs (get_dataset), perform quality-control and construct
a Kraken2 database (autodatabase), and screen reads against this database to produce a high-resolution metagenomic report (screen).

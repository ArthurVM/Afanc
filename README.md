# Afanc: High Resolution Metagenomic Disambiguator
A toolkit which can perform variant level metageonomics disambiguation of NGS reads.

## Installation
Installation is quick and easy, simply run
```
  git clone git@github.com:ArthurVM/Afanc.git
  cd Afanc
  pip3 install ./
```

## Description

### Sub-Modules
Afanc is split into 3 sub-modules:
```
  get_dataset         Download a dataset of genome assemblies from GenBank.
  autodatabase        Generate a database from a FASTA directory structure with autodatabase.
  screen              High-resolution metagenomic screening of short read data using a database
                      constructed by autodatabase.
```

These modules enable the user to download a dataset of genome assemblies using species names or accession IDs with `get_dataset`, perform quality-control and construct a Kraken2 database with `autodatabase`, and screen reads against this database to produce a high-resolution metagenomic report with `screen`.

### get_dataset
This is a general ease-of-use module, allowing the user to provide a list of species IDs (e.g. Escherichia coli) or accessions (e.g. NC_000913.3) in a line separated text file to download and deposit within a directory structure which is acceptable to autodatabase. The user can also specify how many assemblies of each species to download (not available if accessions are provided). 

Given the text file
```
Mycobacterium tuberculosis
Mycobacterium tuberculosis variant bovis BCG
Mycobacterium avium
Mycobacterium simiae
```
the directory structure will be
```
.
|
├── Mycobacterium_avium
│   ├── assembly_1.fa
│   ├── assembly_2.fa
│   └── assembly_3.fa
├── Mycobacterium_simiae
│   ├── assembly_1.fa
│   ├── assembly_2.fa
│   └── assembly_3.fa
└── Mycobacterium_tuberculosis
    ├── assembly_1.fa
    ├── assembly_2.fa
    ├── assembly_3.fa
    └── Mycobacterium_tuberculosis_variant_bovis_BCG
        ├── assembly_1.fa
        ├── assembly_2.fa
        └── assembly_3.fa
```
If accessions are provided, they will just be deposited into a single directory.

### autodatabase
Autodatabase automates the process of constructing a Kraken2 database. This is a pythonic reimagination of the nextflow pipeline https://github.com/annacprice/autodatabase

This module takes a directory structure as described in above, in the get_dataset section. It must contain directories for each species level taxon, where subdirectories within each species directory pertain to subspecies/variants/strains, or any other taxonomic rank lower than species (hereafter referred to simply as variants).

There are six stages to the workflow of this module:

    1) Download the NCBI taxonomy
    2) Add the taxonomic ID to the sequence IDs and the filenames
    3) Create a mash matrix for each taxon
    4) Use the mash matrix to select high quality assemblies for database construction
    6) Build the Kraken2 database
    7) Create a Krona chart showing the composition of the Kraken2 database

By default, it will use the ncbi database from 2020-05-01. If a species or variant is not found within the ncbi database, Afanc will attempt to add it to the database and assign it an ncbi taxonomy ID.

### screen
This module takes a database produced by the autodatabase module, and paired end read data in .fastq format, and performs metagenomic analysis upon it. It produces a report in .json format.

There are five stages to the workflow of this module:
  
    1) Check autodatabase structure is correct
    2) Run Kraken2 on the input dataset
    3) Parse and filter the K2 report to determine the species contained within the dataset, and the most likely variants ("hits")
    4) Map reads to hit assemblies
    5) Construct report
    
## Running Afanc

Running Afanc should, in general, be done in the order of modules presented above. The `get_dataset` module is not necessary if you already have genome assemblies in the directory structure outlined previously.

### Step 1: Create Assembly Directory
```
  afanc get_dataset species_list.txt -n 5 -o my_assemblies_dir
```
This will create a directory structure containing up to 5 (if enough are available on GenBank) assemblies of each species/variant downloaded from GenBank. This can then be fed into the autodatabase module

### Step 2: Create a Database
```
  afanc autodatabase my_assemblies_dir -o my_assemblies_DB
```
This will create a directory structure, which constitutes the database for screening reads against.

### Step 3: Screen Reads
```
  afanc screen my_assemblies_DB my_reads_1.fq.gz my_reads_2.fq.gz -o my_analysis
```
Results will be deposited in a directory structure within `my_analysis`.

## Dependencies
Afanc has a number of dependancies which must be satisfied for full functionality. All software must be in PATH.
```
perl
Python 3.7
Entrez Direct E-utilities
Mash v2.3
Kraken2 v2.1.2
ncbi-blast+
Krona
bowtie2
```
### Entrez Direct
Install instructions for Entrez Direct E-utilities can be found at https://www.ncbi.nlm.nih.gov/books/NBK179288/

### Mash
```
  wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar \
  tar -xf mash-Linux64-v2.3.tar \
  mv mash-Linux64-v2.3/mash /usr/local/bin \
```

### ncbi-blast+
```
  apt-get update
  apt-get install ncbi-blast+
```

### Kraken2
```
  https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
  wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz \
  tar -xzf v2.1.2.tar.gz \
  cd kraken2-2.1.2 \
  ./install_kraken2.sh /usr/local/bin
```

### Krona
```
  git clone https://github.com/marbl/Krona \
  mkdir -p Krona/KronaTools/taxonomy \
  cd /Krona/KronaTools \
  ./install.pl \
  ./updateTaxonomy.sh
```

### Bowtie2
```
  curl -fsSL https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${bowtie2_version}/bowtie2-2.3.4.1-source.zip -o bowtie2-2.3.4.1-source.zip
  unzip bowtie2-2.3.4.1-source.zip 
  make -C bowtie2-2.3.4.1 prefix=/usr/local install
  rm -r bowtie2-2.3.4.1
  rm bowtie2-2.3.4.1-source.zip
```

import re
import pandas as pd
from collections import defaultdict
from os import mkdir, chdir, path, walk


def getTaxidNames(taxname, names):

    # remove underscore from taxname
    taxname = taxname.replace("_", " ")

    # read names.dmp into dataframe
    df = pd.read_csv(names, header=None, sep="|")
    # strip leading and trailing whitespace from column 2
    df[1] = df[1].str.strip()

    # find row with species name
    find_row = df[df[1] == taxname]

    # return the tax ID, throw error if taxon is not found
    if find_row.empty:
        print(f"{taxname} taxon not found!")
        taxid = None
        # raise ValueError("Taxon not found!")
    else:
        taxid = str(find_row.iloc[0, 0])

    return taxid


def editFasta(infasta, outdir, taxname, names):

    # get the taxonomic ID
    taxid = getTaxidNames(taxname, names)

    ## errorStrategy : Ignore equivilent
    ## TODO: add taxa to names and nodes
    if taxid == None:
        return 1

    # edit the headers in each fasta file then output with the taxID
    with open(infasta) as file:
         fastData = file.read()
         fastEdit = re.sub(r"(>.+)",r"\1|kraken:taxid|{0}".format(taxid), fastData)
         fastFile = outdir + "/" + taxid + "_" + path.basename(infasta)
         fastFile = open(fastFile, "w")
         fastFile.write(fastEdit)

    return 0


def taxadd_Main(args, fasta_db_path, names):

    fasta_dict = defaultdict(list)

    for dir, subdirs, fastas in walk(fasta_db_path):
        taxname = dir.split("/")[-1]
        for fasta in fastas:
            infasta = path.join(dir, fasta)
            fasta_dict[taxname].append(infasta)
            editFasta(infasta, args.fasta_WDir, taxname, names)

    return fasta_dict

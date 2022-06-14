import re
import pandas as pd
from collections import defaultdict
from os import mkdir, chdir, path, walk, rename


def editFasta(infasta, outdir, taxid):
    """ Constructs fastas to build the kraken2 database from using the ncbi taxonomy ID
    """

    # edit the headers in each fasta file then output with the taxID
    with open(infasta) as file:
         fasta_data = file.read()
         fasta_edit = re.sub(r"(>.+)",r"\1|kraken:taxid|{0}".format(taxid), fasta_data)
         fasta_file = f"{outdir}/{taxid}_{path.basename(infasta)}"
         fasta_file_fin = open(fasta_file, "w")
         fasta_file_fin.write(fasta_edit)

    return 0


def addTaxon(taxname, mother_clade_taxid, names_df, nodes_path):
    """ takes a taxon missing from the ncbi database, and the names dataframe, and adds the taxon
    to the database.
    """
    max_taxid = max(names_df[0])
    taxid = max_taxid+1
    rank = "no rank"
    # names_dbline = f"\n{taxid}\t|\t{taxname}\t|\t\t|\tscientific name\t|"
    names_dbline = pd.DataFrame([ [taxid, taxname, "", "scientific name", "NaN"] ])
    nodes_dbline = f"\n{taxid}\t|\t{mother_clade_taxid}\t|\t{rank}\t|\t\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|"

    # with open(names_path, 'a') as names_fin:
    #     print(names_dbline, file=names_fin)

    names_df = names_df.append(names_dbline)

    with open(nodes_path, 'a') as nodes_fin:
        print(nodes_dbline, file=nodes_fin)

    print(f"Added {taxname} to ncbi taxonomy database.")

    return str(taxid), names_df


def search_taxon(taxname, names_df):
    """ Takes a taxonomy id and searches for it in the names dataframe
    """
    ## find row with species name
    taxon_row = names_df[names_df[1] == taxname]

    ## return None if taxon cannot be found
    if taxon_row.empty:
        # print(f"{taxname} taxon not found!", end=" ")
        return None
    ## return the ncbi taxon ID number of it can be found
    else:
        return str(taxon_row.iloc[0, 0])


def getTaxidNames(taxname, mother_clade, names_df, nodes_path):
    """ Get ncbi taxonomy ID for this taxon
    """

    # remove underscore from taxname
    taxname = taxname.replace("_", " ")

    taxid = search_taxon(taxname, names_df)

    ## block for dealing with taxon missing from the ncbi taxonomy database
    if taxid == None:

        ## check to see if a mother clade is given
        if mother_clade != None:
            ## find taxid for the mother clade
            mother_clade_unformatted = mother_clade.replace("_", " ")
            # print(f"Attempting to find {mother_clade_unformatted} in ncbi taxonomy database...", end=" ")

            mother_taxid = search_taxon(mother_clade_unformatted, names_df)

            ## if no taxon exists for the mother clade, find taxid for the genus
            if mother_taxid == None:
                genus = mother_clade.split("_")[0]
                # print(f"Attempting to find {genus} in ncbi taxonomy database...")
                genus_taxid = search_taxon(genus, names_df)

                ## if the genus does not exist within the ncbi taxonomy database, then call a fail
                if genus_taxid == None:
                    print(f"{genus} taxon not found! Failed to add to ncbi taxonomy database...")

                ## else if the genus does exist within the database, add both the mother and daughter taxa
                else:
                    # print("Found.", end=" ")
                    mother_taxid, names_df = addTaxon(mother_clade_unformatted, genus_taxid, names_df, nodes_path)
                    taxid, names_df = addTaxon(taxname, mother_taxid, names_df, nodes_path)

            else:
                # print("Found.", end=" ")
                taxid, names_df = addTaxon(taxname, mother_taxid, names_df, nodes_path)

        ## if no mother clade is given, try to find the genus in the taxonomy database
        else:
            genus = taxname.split(" ")[0]
            # print(f"Attempting to find genus {genus} in ncbi taxonomy database...", end=" ")
            genus_taxid = search_taxon(genus, names_df)

            ## if the genus does not exist within the ncbi taxonomy database, then call a fail
            if genus_taxid == None:
                print(f"{genus} taxon not found! Failed to add to ncbi taxonomy database...")

            ## else if the genus does exist within the database, add both the mother and daughter taxa
            else:
                # print("Found.", end=" ")
                taxid, names_df = addTaxon(taxname, genus_taxid, names_df, nodes_path)

    return taxid, names_df


def writeNames(names_df, outpath):
    """ Writes a new names.dmp file
    """

    with open(outpath, "w") as fout:
        for index, row in names_df.iterrows():
            tmp_row = list([str(r) for r in row])[:-1]
            tmp_row.append("")
            newrow = "\t|\t".join(tmp_row)
            print(newrow, file=fout)


def taxadd_Main(fasta_WDir, fasta_db_path, names_path, nodes_path):

    fasta_dict = defaultdict(list)

    # read names.dmp into dataframe
    names_df = pd.read_csv(names_path, header=None, sep="|")
    # strip leading and trailing whitespace from all columns
    names_df = names_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    for dir, subdirs, fastas in walk(fasta_db_path):
        taxon_rank = dir.split(fasta_db_path)[-1].split("/")
        taxname = taxon_rank[-1]

        if len(taxon_rank) > 1:
            mother_clade = taxon_rank[-2]
        else:
            mother_clade = None

        for fasta in fastas:
            infasta = path.join(dir, fasta)
            fasta_dict[taxname].append(infasta)

            # get the taxonomic ID
            taxid, names_df = getTaxidNames(taxname, mother_clade, names_df, nodes_path)

            ## errorStrategy : Ignore equivilent
            ## TODO: add taxa to names and nodes
            if taxid == None:
                continue

            editFasta(infasta, fasta_WDir, taxid)

    ## rename the old names.dmp file
    dirpath = "/".join(names_path.split("/")[:-1])
    rename(names_path, f"{dirpath}/names.dmp.old")
    ## write the new nodes.dmp file
    writeNames(names_df, f"{dirpath}/names.dmp")

    return fasta_dict

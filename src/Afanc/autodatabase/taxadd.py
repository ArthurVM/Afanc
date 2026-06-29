import re
import gzip
import pandas as pd
from Bio import SeqIO, SeqRecord
from collections import defaultdict
from os import mkdir, chdir, path, walk, rename
from urllib.parse import unquote


FASTA_EXTS = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
NCBI_ALIAS_NAME_CLASSES = {"scientific name", "synonym", "equivalent name", "authority"}


def normaliseTaxonName(taxname):
    """ Normalise directory-derived taxon names for NCBI lookup."""
    taxname = unquote(str(taxname))
    taxname = taxname.replace("_", " ")
    taxname = re.sub(r"\s+", " ", taxname)

    return taxname.strip()


def editFasta(infasta, outdir, taxid, taxname):
    """ Constructs fastas to build the kraken2 database from using the ncbi taxonomy ID
    """

    if infasta.endswith(".gz"):
        fasta_fin = gzip.open(infasta, "rt")
        fasta_file = f"{outdir}/{taxid}_{path.basename(infasta).split('.gz')[0]}"
    else:
        fasta_fin = open(infasta, "r")
        fasta_file = f"{outdir}/{taxid}_{path.basename(infasta)}"

    fasta_fout = open(fasta_file, "w")

    for rec in SeqIO.parse(fasta_fin, "fasta"):
        header = f">{rec.id}|kraken:taxid|{taxid} {taxname}\n"
        fasta_fout.write(header)
        fasta_fout.write(str(rec.seq)+"\n")

    fasta_fin.close()
    fasta_fout.close()

    return 0


def addTaxon(taxname, mother_clade_taxid, names_df, nodes_df, name_index=None, taxadd_state=None):
    """ takes a taxon missing from the ncbi database, and the names dataframe, and adds the taxon
    to the database.
    """

    if taxadd_state is not None:
        taxid = taxadd_state["next_taxid"]
        taxadd_state["next_taxid"] += 1
        taxadd_state["taxonomy_changed"] = True
    else:
        max_taxid = max(names_df[0])
        taxid = int(max_taxid + 1)

    ## infer rank from binomial or trinomial name
    if taxname.count(" ")>1:
        rank = "subspecies"
    elif taxname.count(" ")==1:
        rank = "species"
    else:
        rank = "no rank"

    names_dbline = pd.DataFrame([ [taxid, taxname, " ", "scientific name"] ])
    nodes_dbline = pd.DataFrame([ [taxid, mother_clade_taxid, rank, "", "1", "1", "1", "1", "1" ,"1", "1", "1", " "] ])

    names_df = pd.concat([names_df, names_dbline], ignore_index=True)
    nodes_df = pd.concat([nodes_df, nodes_dbline], ignore_index=True)

    if name_index is not None:
        name_index[normaliseTaxonName(taxname)] = str(taxid)

    print(f"Added {taxname} to ncbi taxonomy database.")

    return str(taxid), names_df, nodes_df


def buildNameIndex(names_df):
    """ Build a fast taxon-name:taxID lookup from names.dmp."""
    name_index = {}

    for taxid, taxname in zip(names_df[0], names_df[1]):
        normalised_taxname = normaliseTaxonName(taxname)

        if normalised_taxname == "":
            continue

        ## keep the first taxid for duplicate names
        if normalised_taxname not in name_index:
            name_index[normalised_taxname] = str(taxid)

    return name_index


def taxonRank(taxname):
    """Infer the expected NCBI rank from a normalized input taxon name."""
    tokens = normaliseTaxonName(taxname).split()

    if len(tokens) == 2:
        return "species"

    return None


def taxonPrefix(taxname):
    """Return the conservative taxon-name prefix used for authority matching."""
    tokens = normaliseTaxonName(taxname).split()

    if len(tokens) == 2:
        return " ".join(tokens[:2])

    return ""


def getTaxonRank(taxid, nodes_df):
    """Return the NCBI rank for a taxID from nodes.dmp data."""
    taxid = str(taxid)
    taxon_row = nodes_df[nodes_df[0].astype(str) == taxid]

    if taxon_row.empty:
        return None

    return str(taxon_row.iloc[0, 2])


def getScientificName(taxid, names_df):
    """Return the scientific name for a taxID from names.dmp data."""
    taxid = str(taxid)
    taxon_rows = names_df[
        (names_df[0].astype(str) == taxid)
        & (names_df[3].astype(str) == "scientific name")
    ]

    if taxon_rows.empty:
        return None

    return normaliseTaxonName(taxon_rows.iloc[0, 1])


def searchTaxonAuthorityPrefix(taxname, names_df, nodes_df):
    """Resolve authority-bearing NCBI names before adding a novel taxon.

    NCBI often stores old genus combinations as authority names, e.g.
    "Mycobacterium chelonae Bergey ...", while the scientific name is now
    "Mycobacteroides chelonae". This fallback only runs for species-level
    binomials to avoid collapsing below-species labels such as lineages onto
    unrelated authority names.
    """
    prefix = taxonPrefix(taxname)
    expected_rank = taxonRank(taxname)

    if prefix == "" or expected_rank is None:
        return None

    candidate_taxids = []
    prefix_with_space = f"{prefix} "

    for taxid, ncbi_name, name_class in zip(names_df[0], names_df[1], names_df[3]):
        if str(name_class) not in NCBI_ALIAS_NAME_CLASSES:
            continue

        normalised_name = normaliseTaxonName(ncbi_name)
        if normalised_name != prefix and not normalised_name.startswith(prefix_with_space):
            continue

        if getTaxonRank(taxid, nodes_df) != expected_rank:
            continue

        candidate_taxids.append(str(taxid))

    candidate_taxids = sorted(set(candidate_taxids), key=lambda x: int(x))

    if len(candidate_taxids) != 1:
        return None

    return candidate_taxids[0]


def searchTaxon(taxname, names_df, name_index=None):
    """ Takes a taxonomy id and searches for it in the names dataframe
    """
    taxname = normaliseTaxonName(taxname)

    if name_index is not None:
        return name_index.get(taxname)

    taxon_row = names_df[names_df[1] == taxname]

    if taxon_row.empty:
        return None
    else:
        return str(taxon_row.iloc[0, 0])


def getTaxidNames(taxname, mother_clade, names_df, nodes_df, name_index=None, taxadd_state=None):
    """ Get ncbi taxonomy ID for this taxon
    """

    taxname = normaliseTaxonName(taxname)

    taxid = searchTaxon(taxname, names_df, name_index=name_index)

    if taxid == None:
        taxid = searchTaxonAuthorityPrefix(taxname, names_df, nodes_df)
        if taxid is not None:
            scientific_name = getScientificName(taxid, names_df)
            print(
                f"Resolved {taxname} to existing NCBI taxID {taxid}"
                f" ({scientific_name}) using an authority/synonym name."
            )

    if taxid == None:
        # print(f"Cannot find {taxname} in ncbi taxonomy database.")

        ###  TAXADD BLOCK  ###
        ### IN DEVELOPMENT ###

        ## check to see if a mother clade is given
        if mother_clade != None and mother_clade.strip(" ") != "":
            mother_clade_unformatted = normaliseTaxonName(mother_clade)
            print(f"Attempting to find {mother_clade_unformatted} in ncbi taxonomy database...", end=" ")
            mother_taxid = searchTaxon(mother_clade_unformatted, names_df, name_index=name_index)

            if mother_taxid == None:
                genus = normaliseTaxonName(mother_clade).split(" ")[0]
                print(f"Attempting to find {genus} in ncbi taxonomy database...")
                genus_taxid = searchTaxon(genus, names_df, name_index=name_index)

                if genus_taxid == None:
                    print(f"{genus} taxon not found! Failed to add to ncbi taxonomy database...")

                else:
                    # print(f"Found {genus_taxid}.", end="\n")
                    mother_taxid, names_df, nodes_df = addTaxon(mother_clade_unformatted, genus_taxid, names_df, nodes_df, name_index=name_index, taxadd_state=taxadd_state)
                    taxid, names_df, nodes_df = addTaxon(taxname, mother_taxid, names_df, nodes_df, name_index=name_index, taxadd_state=taxadd_state)

            else:
                # print(f"Found {mother_taxid}.", end="\n")
                taxid, names_df, nodes_df = addTaxon(taxname, mother_taxid, names_df, nodes_df, name_index=name_index, taxadd_state=taxadd_state)

        else:
            genus = taxname.split(" ")[0]
            # print(f"Attempting to find genus {genus} in ncbi taxonomy database...", end=" ")
            genus_taxid = searchTaxon(genus, names_df, name_index=name_index)

            if genus_taxid == None:
                print(f"{genus} taxon not found! Failed to add to ncbi taxonomy database...")

            else:
                taxid, names_df, nodes_df = addTaxon(taxname, genus_taxid, names_df, nodes_df, name_index=name_index, taxadd_state=taxadd_state)


    return taxid, names_df, nodes_df


def writeDmp(names_df, nodes_df, outdir):
    """ Writes a new dmp files
    """

    with open(f"{outdir}/names.dmp", 'w') as fout:
        for index, row in names_df.iterrows():
            tmp_row = list([str(r) for r in row])[:-1]
            newrow = "\t|\t".join(tmp_row)
            newrow += "\t|"
            print(newrow, file=fout)

    with open(f"{outdir}/nodes.dmp", 'w') as fout:
        for index, row in nodes_df.iterrows():
            tmp_row = list([str(r) for r in row])[:-1]
            newrow = "\t|\t".join(tmp_row)
            newrow += "\t|"
            print(newrow, file=fout)


def readDmp(dmp_file):
    """ Read in nodes/names dmp files and return as a pandas dataframe
    """

    ## strip NCBI field padding before rewriting
    dmp_df = pd.read_csv(dmp_file, header=None, sep="|", dtype=str, keep_default_na=False)

    dmp_df = dmp_df.apply(lambda col: col.map(lambda value: value.strip() if isinstance(value, str) else value))
    dmp_df[0] = pd.to_numeric(dmp_df[0], errors="coerce").astype("Int64")

    return dmp_df


def taxadd_Main(fasta_WDir, fasta_db_path, names_path, nodes_path):

    scaffold_taxid_dict = defaultdict(list)

    fasta_dict = defaultdict(list)

    mapping_dict = defaultdict(str)

    names_df = readDmp(names_path)
    nodes_df = readDmp(nodes_path)
    name_index = buildNameIndex(names_df)
    taxadd_state = {
        "next_taxid": int(max(names_df[0])) + 1,
        "taxonomy_changed": False,
    }

    for dir, subdirs, fastas in walk(fasta_db_path):
        rel_dir = path.relpath(dir, fasta_db_path)

        if rel_dir == ".":
            taxon_rank = []
            taxname = path.basename(path.abspath(dir))
        else:
            taxon_rank = rel_dir.split(path.sep)
            taxname = taxon_rank[-1]

        if len(taxon_rank) > 1:
            mother_clade = taxon_rank[-2]
        else:
            mother_clade = None

        for fasta in fastas:
            if not fasta.endswith(FASTA_EXTS):
                continue

            infasta = path.join(dir, fasta)
            fasta_dict[taxname].append(infasta)

            taxid, names_df, nodes_df = getTaxidNames(taxname, mother_clade, names_df, nodes_df, name_index=name_index, taxadd_state=taxadd_state)

            if taxid == None:
                continue

            mapping_dict[taxid] = normaliseTaxonName(taxname)

            editFasta(infasta, fasta_WDir, taxid, normaliseTaxonName(taxname))

    dirpath = path.dirname(names_path)

    if taxadd_state["taxonomy_changed"]:
        rename(names_path, f"{dirpath}/names.dmp.backup")
        rename(nodes_path, f"{dirpath}/nodes.dmp.backup")

        writeDmp(names_df, nodes_df, dirpath)

    return fasta_dict, mapping_dict

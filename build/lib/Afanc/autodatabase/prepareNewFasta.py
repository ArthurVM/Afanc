from os import path, mkdir, chdir
from Afanc.utilities.runCommands import command


def getTaxonomy(args, ncbiDate):
    """
    Download the NCBI taxonomy, identified using args.ncbiDate

    INPUT:
        ncbiDate <str> : the date of the ncbi taxonomy to download

    OUTPUT:
        names <str> : names.dmp path
        nodes <str> : nodes.dmp path
    """

    taxonomy_ftp_path = f"ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_{ncbiDate}.zip"

    get_taxonomy_runline = f"wget {taxonomy_ftp_path} ; unzip -o taxdmp_{ncbiDate}.zip"
    stdout, stderr = command(get_taxonomy_runline, "GET_TAXONOMY").run_comm(1, args.stdout, args.stderr)

    names = path.abspath("./names.dmp")
    nodes = path.abspath("./nodes.dmp")

    ## add M. tomidae
    addTaxon(args, names, nodes, "Mycobacterium tomidae", 120793)

    return names, nodes


def addTaxon(args, names, nodes, taxon, taxa_number):
    """ Add a taxon to nodes.dmp and names.dmp

    INPUT:
        names <str> : names.dmp path
        nodes <str> : nodes.dmp path

    OUTPUT:
        None

    TODO:
        Generalise this function for any taxon, rather than just Mycobacterium.
        This with require automating the process of identifying the taxa_number
    """
    largestTax, stderr = command(f"sort -t't' -k1nr {names} | head -1 | cut -f1", "GET_TAXONOMY").run_comm_quiet(1, args.stdout, args.stderr)
    taxID = int(largestTax)+1

    with open(names, 'a') as names_fout, open(nodes, 'a') as nodes_fout:
        names_fout.write(f"{taxID}\t|\t{taxon}\t|\t\t|\tscientific name\t|")
        nodes_fout.write(f"{taxID}\t|\t{taxa_number}\t|\tspecies\t|\t\t|\t11\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|\t1\t|")

import json
import numpy as np

from collections import defaultdict
from os import path, listdir

from Afanc.utilities.runCommands import command
from Afanc.screen.report.parseK2report import readK2report, parseK2line


def fastANI(args, fastas_prefix):
    """ run fastANI to calculate parent-child distance.
    """
    outfile = f"{fastas_prefix}.fastANI.txt"
    runline = f"fastANI --rl {fastas_prefix}.fastas.txt --ql {fastas_prefix}.fastas.txt -t 4 -o {outfile}"
    stdout, stderr = command(runline, "VARIANT-INDEX").run_comm_quiet(1, args.stdout, args.stderr)

    return read_fastANI_out(outfile)


def read_fastANI_out(fastANI_out):
    """ Read the output file from fastANI
    """

    child_taxID = fastANI_out.split(".")[0]

    parent_similarity = []
    sibling_similarity = []

    for r, q, similarity, _, _ in read_csv(fastANI_out, delimiter="\t"):
        r = r.split("/")[-1]
        q = q.split("/")[-1]
        r_taxID = r.split("_")[0]
        q_taxID = q.split("_")[0]

        if child_taxID == r_taxID and r_taxID != q_taxID:
            parent_similarity.append(float(similarity))
        elif child_taxID == q_taxID and r_taxID == q_taxID:
            sibling_similarity.append(float(similarity))

    # print(fastANI_out, np.mean(parent_similarity), np.mean(sibling_similarity))

    return parent_similarity, sibling_similarity


def read_csv(csv, delimiter=","):
    """ simple generator function returning lines from a csv file.
    """
    with open(csv, 'r') as fin:
        for line in fin.readlines():
            yield line.strip("\n").split(delimiter)


def make_variant_index(args, base_nodes):
    """ Construct the variant index using fastANI results
    """

    fasta_dict = defaultdict(list)

    scores_dict = { "variant_index" : {} }

    ## make taxID : [ fastas ] dictionary
    for fasta in listdir(args.cleanFasta_WDir):
        taxID = int(fasta.split("_")[0])
        fasta_dict[taxID].append(path.join(args.cleanFasta_WDir, fasta))

    for taxID, fastabox in fasta_dict.items():

        ## capture warnings
        warning_box = []

        if base_nodes[taxID].parent.ncbi_taxID in fasta_dict:
            fastas_prefix = f"{taxID}.{base_nodes[taxID].name.replace(' ', '_')}"
            with open(fastas_prefix + ".fastas.txt", 'w') as fout:
                print("\n".join(fastabox), file=fout)
                print("\n".join(fasta_dict[base_nodes[taxID].parent.ncbi_taxID]), file=fout)

            parent_similarity, sibling_similarity = fastANI(args, fastas_prefix)

            # print(parent_similarity, sibling_similarity)

            parent_mean = np.mean(parent_similarity)
            sibling_mean = np.mean(sibling_similarity)

            if parent_mean >= sibling_mean:
                warning = f"Intra-taxon variation for {base_nodes[taxID].name} ({parent_mean}) exceeds mean distance from parent taxon {base_nodes[taxID].parent.name} ({sibling_mean}).\
                            This may lead to erroneous results within this clade."
                print(warning)
                warning_box.append(warning)

            scores_dict["variant_index"][taxID] = { "name" : str(base_nodes[taxID].name),
                                                    "parent" : str(base_nodes[taxID].parent.name),
                                                    "parent_index" : { "sim_array" : parent_similarity, "mean" : parent_mean, "median" : np.median(parent_similarity), "range" : np.max(parent_similarity) - np.min(parent_similarity) },
                                                    "sibling_index" : { "sim_array" : sibling_similarity, "mean" : sibling_mean, "median" : np.median(sibling_similarity), "range" : np.max(sibling_similarity) - np.min(sibling_similarity),
                                                    "warnings" : warning_box }
                                                  }

    with open(f"{args.autoDB_WDir}/variant_index.json", 'w') as fout:
        json.dump(scores_dict, fout, indent = 4)

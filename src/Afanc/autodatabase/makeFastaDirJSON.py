import json
from collections import defaultdict
from os import path, walk

def make_fasta_dir_JSON(fasta_dir):
    """ Takes the inpit FASTA directory structure and generates a JSON listing
    files and their hierarchical structure.
    """

    json_dict = defaultdict(list)

    for dir, subdirs, files in walk(fasta_dir):

        if len(files) == 0:
            continue

        taxon = path.basename(dir)

        ## remove base paths
        filenames = [path.basename(f) for f in files]
        json_dict[taxon] = files

    with open("assemblies.json", "w") as fout:
        json.dump(json_dict, fout, indent=4)

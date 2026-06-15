import csv
import json
import numpy as np

from collections import defaultdict
from os import listdir, makedirs, path

from Afanc.utilities.runCommands import command


def read_dmp(dmp_file):
    """Read an NCBI .dmp file."""
    dmp_rows = []

    with open(dmp_file, "r") as fin:
        for line in fin.readlines():
            fields = [field.strip() for field in line.split("|")]

            while len(fields) > 0 and fields[-1] == "":
                fields.pop()

            dmp_rows.append(fields)

    return dmp_rows


def read_csv(csv_file, delimiter="\t"):
    """Yield rows from a delimited file."""
    with open(csv_file, "r") as fin:
        for line in fin.readlines():
            yield line.strip("\n").split(delimiter)


def build_taxonomy(names_dmp, nodes_dmp):
    """Build taxonomy lookup dictionaries from NCBI dmp files."""
    names_dict = {}
    parent_dict = {}
    rank_dict = {}

    for row in read_dmp(names_dmp):
        if len(row) >= 4 and row[3] == "scientific name":
            names_dict[int(row[0])] = row[1]

    for row in read_dmp(nodes_dmp):
        if len(row) >= 3:
            tax_id = int(row[0])
            parent_dict[tax_id] = int(row[1])
            rank_dict[tax_id] = row[2]

    return names_dict, parent_dict, rank_dict


def get_taxon_name(names_dict, tax_id):
    """Get a display name for a taxID."""
    return str(names_dict.get(tax_id, tax_id))


def get_ancestor_by_rank(tax_id, parent_dict, rank_dict, target_rank):
    """Walk up a lineage until the requested rank is found."""
    current_tax_id = tax_id

    while current_tax_id in parent_dict:
        if rank_dict.get(current_tax_id) == target_rank:
            return current_tax_id

        parent_tax_id = parent_dict[current_tax_id]

        if parent_tax_id == current_tax_id:
            break

        current_tax_id = parent_tax_id

    return None


def collect_fasta_dict(fasta_dir):
    """Collect a taxID:fasta mapping from fasta_dir filenames."""
    fasta_dict = defaultdict(list)
    fasta_exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")

    for fasta_name in sorted(listdir(fasta_dir)):
        if not fasta_name.endswith(fasta_exts):
            continue

        tax_id_field = fasta_name.split("_", 1)[0]

        try:
            tax_id = int(tax_id_field)
        except ValueError:
            continue

        fasta_dict[tax_id].append(path.abspath(path.join(fasta_dir, fasta_name)))

    return fasta_dict


def write_fasta_list(fasta_box, fasta_list_file):
    """Write a Mash fasta list file."""
    with open(fasta_list_file, "w") as fout:
        print("\n".join(fasta_box), file=fout)


def run_mash_all_vs_all(
    mash_bin,
    fasta_box,
    out_dir,
    prefix="all_genomes",
    threads=4,
    kmer=16,
    sketch_size=10000,
):
    """Run mash sketch and mash dist for all genomes."""
    makedirs(out_dir, exist_ok=True)

    fasta_list_file = path.join(out_dir, f"{prefix}.fastas.txt")
    sketch_prefix = path.join(out_dir, prefix)
    mash_out = path.join(out_dir, f"{prefix}.mash.txt")

    write_fasta_list(fasta_box, fasta_list_file)

    sketch_cmd = [
        str(mash_bin), "sketch",
        "-o", sketch_prefix,
        "-k", str(kmer),
        "-s", str(sketch_size),
        "-p", str(threads),
        "-l", fasta_list_file,
    ]

    dist_cmd = [
        str(mash_bin), "dist",
        "-p", str(threads),
        f"{sketch_prefix}.msh",
        f"{sketch_prefix}.msh",
    ]

    command(sketch_cmd, "MASH").run_comm_quiet(0)
    dist_stdout, _ = command(dist_cmd, "MASH").run_comm_quiet(1)

    with open(mash_out, "wb") as fout:
        fout.write(dist_stdout)

    return mash_out


def read_mash_out(mash_out):
    """Read Mash distance output."""
    distance_dict = {}

    for ref_fasta, query_fasta, distance, _, _ in read_csv(mash_out, delimiter="\t"):
        distance_dict[(path.abspath(ref_fasta), path.abspath(query_fasta))] = float(distance)

    return distance_dict


def get_matrix_value(distance_dict, fasta_a, fasta_b):
    """Fetch a Mash distance from the distance dictionary."""
    if fasta_a == fasta_b:
        return 0.0

    if (fasta_a, fasta_b) in distance_dict:
        return distance_dict[(fasta_a, fasta_b)]

    if (fasta_b, fasta_a) in distance_dict:
        return distance_dict[(fasta_b, fasta_a)]

    return np.nan


def make_genome_box(fasta_dict, names_dict, parent_dict, rank_dict):
    """Construct genome records from the fasta dictionary and taxonomy."""
    genome_box = []

    for tax_id, fasta_box in fasta_dict.items():
        genus_tax_id = get_ancestor_by_rank(tax_id, parent_dict, rank_dict, "genus")
        species_tax_id = get_ancestor_by_rank(tax_id, parent_dict, rank_dict, "species")

        for fasta in fasta_box:
            genome_box.append({
                "tax_id": tax_id,
                "name": get_taxon_name(names_dict, tax_id),
                "fasta": path.abspath(fasta),
                "label": path.basename(fasta),
                "genus_tax_id": genus_tax_id,
                "species_tax_id": species_tax_id,
            })

    return genome_box


def write_matrix(matrix_file, labels, matrix, header):
    """Write a square matrix to TSV."""
    with open(matrix_file, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow([header] + labels)

        for label, row in zip(labels, matrix):
            out_row = [label]

            for value in row:
                if np.isnan(value):
                    out_row.append("NA")
                else:
                    out_row.append(f"{value:.6f}")

            writer.writerow(out_row)


def write_rank_heterogeneity(heterogeneity_file, heterogeneity_rows, header):
    """Write within-group heterogeneity summaries to TSV."""
    with open(heterogeneity_file, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow([header, "n_genomes", "within_group_mean_distance"])

        for row in heterogeneity_rows:
            if np.isnan(row["within_group_mean_distance"]):
                mean_distance = "NA"
            else:
                mean_distance = f"{row['within_group_mean_distance']:.6f}"

            writer.writerow([row["label"], row["n_genomes"], mean_distance])


def make_all_vs_all_matrix(genome_box, distance_dict):
    """Construct a genome all-vs-all Mash distance matrix."""
    labels = [genome["label"] for genome in genome_box]
    matrix = []

    for genome_i in genome_box:
        row = []

        for genome_j in genome_box:
            row.append(get_matrix_value(distance_dict, genome_i["fasta"], genome_j["fasta"]))

        matrix.append(row)

    return labels, matrix


def get_rank_members(genome_box, rank):
    """Group genomes by nearest genus or species ancestor."""
    members_dict = defaultdict(list)

    for genome in genome_box:
        rank_tax_id = genome[f"{rank}_tax_id"]

        if rank_tax_id is None:
            continue

        members_dict[rank_tax_id].append(genome)

    return members_dict


def get_within_group_mean_distance(members, distance_dict):
    """Compute mean within-group distance, excluding self comparisons."""
    pair_box = []

    for genome_i in members:
        for genome_j in members:
            if genome_i["fasta"] == genome_j["fasta"]:
                continue

            distance = get_matrix_value(distance_dict, genome_i["fasta"], genome_j["fasta"])

            if not np.isnan(distance):
                pair_box.append(distance)

    if pair_box:
        return float(np.mean(pair_box))

    return np.nan


def make_rank_matrix(genome_box, names_dict, distance_dict, rank):
    """Construct rank-level between-group matrices and within-group heterogeneity summaries."""
    members_dict = get_rank_members(genome_box, rank)
    rank_tax_ids = sorted(members_dict)
    labels = [f"{tax_id}|{get_taxon_name(names_dict, tax_id)}" for tax_id in rank_tax_ids]
    matrix = []
    heterogeneity_rows = []

    for tax_id, label in zip(rank_tax_ids, labels):
        heterogeneity_rows.append({
            "tax_id": tax_id,
            "label": label,
            "n_genomes": len(members_dict[tax_id]),
            "within_group_mean_distance": get_within_group_mean_distance(members_dict[tax_id], distance_dict),
        })

    for tax_id_i in rank_tax_ids:
        row = []

        for tax_id_j in rank_tax_ids:
            if tax_id_i == tax_id_j:
                row.append(0.0)
                continue

            pair_box = []

            for genome_i in members_dict[tax_id_i]:
                for genome_j in members_dict[tax_id_j]:
                    distance = get_matrix_value(distance_dict, genome_i["fasta"], genome_j["fasta"])

                    if not np.isnan(distance):
                        pair_box.append(distance)

            if pair_box:
                row.append(float(np.mean(pair_box)))
            else:
                row.append(np.nan)

        matrix.append(row)

    return labels, matrix, heterogeneity_rows


def summarise_similarity(similarity_box):
    """Summarise a similarity array for variant index output."""
    if len(similarity_box) == 0:
        return {"sim_array": [], "mean": None, "median": None, "range": None}

    sim_array = [float(value) for value in similarity_box]

    return {
        "sim_array": sim_array,
        "mean": float(np.mean(sim_array)),
        "median": float(np.median(sim_array)),
        "range": float(np.max(sim_array) - np.min(sim_array)),
    }


def make_variant_index(fasta_dict, names_dict, parent_dict, distance_dict):
    """Construct the Afanc-style parent-child similarity index."""
    variant_index = {"variant_index": {}}

    for tax_id, fasta_box in fasta_dict.items():
        if tax_id not in parent_dict:
            continue

        parent_tax_id = parent_dict[tax_id]

        if parent_tax_id == tax_id:
            continue

        if parent_tax_id not in fasta_dict:
            continue

        parent_similarity = []
        sibling_similarity = []
        warning_box = []

        for child_fasta in fasta_box:
            for parent_fasta in fasta_dict[parent_tax_id]:
                distance = get_matrix_value(distance_dict, child_fasta, parent_fasta)

                if not np.isnan(distance):
                    parent_similarity.append(100 - float(distance))

        for ref_fasta in fasta_box:
            for query_fasta in fasta_box:
                distance = get_matrix_value(distance_dict, ref_fasta, query_fasta)

                if not np.isnan(distance):
                    sibling_similarity.append(100 - float(distance))

        parent_index = summarise_similarity(parent_similarity)
        sibling_index = summarise_similarity(sibling_similarity)

        if parent_index["mean"] is not None and sibling_index["mean"] is not None:
            if parent_index["mean"] >= sibling_index["mean"]:
                warning = f"Mean parent-child similarity for {get_taxon_name(names_dict, tax_id)} ({parent_index['mean']}) exceeds mean intrataxon similarity ({sibling_index['mean']}). This clade may behave poorly during read commuting."
                print(warning)
                warning_box.append(warning)

        sibling_index["warnings"] = warning_box

        variant_index["variant_index"][str(tax_id)] = {
            "name": get_taxon_name(names_dict, tax_id),
            "parent": get_taxon_name(names_dict, parent_tax_id),
            "parent_index": parent_index,
            "sibling_index": sibling_index,
        }

    return variant_index


def write_variant_index(variant_index, out_json):
    """Write the parent-child similarity index to JSON."""
    with open(out_json, "w") as fout:
        json.dump(variant_index, fout, indent=4)


def make_mash_taxonomic_matrices(
    fasta_dir,
    names_dmp,
    nodes_dmp,
    mash_bin,
    out_dir,
    threads=4,
    kmer=16,
    sketch_size=10000,
    prefix="all_genomes",
):
    """Run the Mash workflow from fasta_dir and NCBI dmp taxonomy."""
    fasta_dict = collect_fasta_dict(fasta_dir)
    names_dict, parent_dict, rank_dict = build_taxonomy(names_dmp, nodes_dmp)
    genome_box = make_genome_box(fasta_dict, names_dict, parent_dict, rank_dict)
    fasta_box = [genome["fasta"] for genome in genome_box]

    mash_out = run_mash_all_vs_all(
        mash_bin,
        fasta_box,
        out_dir,
        prefix=prefix,
        threads=threads,
        kmer=kmer,
        sketch_size=sketch_size,
    )

    distance_dict = read_mash_out(mash_out)

    genome_labels, genome_matrix = make_all_vs_all_matrix(genome_box, distance_dict)
    genus_labels, genus_matrix, genus_heterogeneity = make_rank_matrix(genome_box, names_dict, distance_dict, "genus")
    species_labels, species_matrix, species_heterogeneity = make_rank_matrix(genome_box, names_dict, distance_dict, "species")
    variant_index = make_variant_index(fasta_dict, names_dict, parent_dict, distance_dict)

    write_matrix(path.join(out_dir, "mash_all_vs_all.tsv"), genome_labels, genome_matrix, header="sample")
    write_matrix(path.join(out_dir, "mash_genus_matrix.tsv"), genus_labels, genus_matrix, header="genus")
    write_matrix(path.join(out_dir, "mash_species_matrix.tsv"), species_labels, species_matrix, header="species")
    write_rank_heterogeneity(path.join(out_dir, "mash_genus_heterogeneity.tsv"), genus_heterogeneity, header="genus")
    write_rank_heterogeneity(path.join(out_dir, "mash_species_heterogeneity.tsv"), species_heterogeneity, header="species")
    write_variant_index(variant_index, path.join(out_dir, "variant_index.json"))

    return {
        "mash_out": mash_out,
        "distance_dict": distance_dict,
        "genome_box": genome_box,
        "mash_all_vs_all": {"labels": genome_labels, "matrix": genome_matrix},
        "mash_genus_matrix": {"labels": genus_labels, "matrix": genus_matrix},
        "mash_genus_heterogeneity": genus_heterogeneity,
        "mash_species_matrix": {"labels": species_labels, "matrix": species_matrix},
        "mash_species_heterogeneity": species_heterogeneity,
        "variant_index": variant_index,
    }

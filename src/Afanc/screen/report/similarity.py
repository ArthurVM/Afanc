import csv
import json
from pathlib import Path


def read_mvi_variant_index(mvi):
    """Read variant_index.json from an MVI directory or mapping."""
    if isinstance(mvi, dict):
        if "variant_index" in mvi:
            return mvi["variant_index"]
        return mvi

    mvi_path = Path(mvi)
    variant_index_path = mvi_path / "variant_index.json" if mvi_path.is_dir() else None
    if variant_index_path is None or not variant_index_path.exists():
        raise ValueError("variant_index was not supplied and variant_index.json was not found in the MVI.")

    with open(variant_index_path, "r") as fin:
        data = json.load(fin)

    return data["variant_index"]


def read_mvi_mash_distances(mvi):
    """Read Mash distance matrices from an MVI path into a taxID-pair lookup.

    MVI input may be a TSV matrix, a directory containing Mash matrix TSVs, or a
    prebuilt metric lookup. Matrix values are native Mash distances and are
    converted to weighting similarities as ``1 - distance``.
    """
    if isinstance(mvi, dict):
        return mvi

    matrix_path = Path(mvi)
    if matrix_path.is_dir():
        matrix_files = [
            path for path in sorted(matrix_path.glob("*.tsv"))
            if "matrix" in path.name or "all_vs_all" in path.name
        ]
    else:
        matrix_files = [matrix_path]

    lookup = {}
    for matrix_file in matrix_files:
        read_mash_matrix_file(matrix_file, lookup)

    return lookup


def read_mash_matrix_file(matrix_file, lookup):
    with open(matrix_file, "r", newline="") as fin:
        reader = csv.reader(fin, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return

        col_labels = header[1:]
        col_taxids = [taxid_from_matrix_label(label) for label in col_labels]

        for row in reader:
            if not row:
                continue

            row_taxid = taxid_from_matrix_label(row[0])
            if row_taxid is None:
                continue

            for col_taxid, value in zip(col_taxids, row[1:]):
                if col_taxid is None or value in {None, ""}:
                    continue

                raw_value = float(value)
                metric = metric_from_mash_distance(raw_value, source=str(matrix_file))
                store_best_metric(lookup, row_taxid, col_taxid, metric)


def taxid_from_matrix_label(label):
    label = str(label).strip()
    if not label:
        return None

    if "|" in label:
        label = label.split("|", 1)[0]
    else:
        label = Path(label).name
        label = label.split("_", 1)[0]

    try:
        return int(label)
    except ValueError:
        return None


def metric_from_mash_distance(value, source):
    distance = max(float(value), 0.0)
    similarity = max(1.0 - distance, 0.0)

    return {
        "similarity": similarity,
        "distance": distance,
        "source": source,
    }


def store_best_metric(lookup, taxid_a, taxid_b, metric):
    key = (int(taxid_a), int(taxid_b))
    reverse_key = (int(taxid_b), int(taxid_a))

    existing = lookup.get(key)
    if existing is None or metric["similarity"] > existing["similarity"]:
        lookup[key] = metric
        lookup[reverse_key] = metric


def similarity_between_subtrees(donor, recipient, similarity_lookup):
    """Return direct or best descendant-pair similarity metric for two nodes."""
    donor_taxid = int(donor.ncbi_taxID)
    recipient_taxid = int(recipient.ncbi_taxID)

    if donor_taxid == recipient_taxid:
        return {"similarity": 1.0, "distance": 0.0, "source": "self"}

    direct = similarity_lookup.get((donor_taxid, recipient_taxid))
    if direct is not None:
        return direct

    metrics = []
    for donor_leaf in donor.leaf_descendants():
        for recipient_leaf in recipient.leaf_descendants():
            metric = similarity_lookup.get((int(donor_leaf.ncbi_taxID), int(recipient_leaf.ncbi_taxID)))
            if metric is not None:
                metrics.append(metric)

    if metrics:
        return max(metrics, key=lambda metric: metric["similarity"])

    return None

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict


def load_allele_json(input_json: Path) -> Dict[str, Any]:
    with Path(input_json).open("r") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("Allele JSON must be an object.")

    if "alleles" in payload:
        alleles = payload["alleles"]
    elif "allele" in payload:
        alleles = payload["allele"]
    else:
        raise ValueError("Allele JSON must contain 'allele' or 'alleles'.")

    missing = payload.get("missing", [])
    if not isinstance(alleles, list):
        raise ValueError("Allele JSON 'allele'/'alleles' field must be a list.")
    if not isinstance(missing, list):
        raise ValueError("Allele JSON 'missing' field must be a list.")

    normalised_missing = []
    for item in missing:
        if isinstance(item, str):
            parts = item.rsplit(".", 1)
            if len(parts) != 2:
                raise ValueError("String missing positions must have form 'chrom.pos'.")
            normalised_missing.append([parts[0], parts[1]])
        elif isinstance(item, (list, tuple)) and len(item) == 2:
            normalised_missing.append([str(item[0]), str(item[1])])
        else:
            raise ValueError("Missing positions must be [chrom, pos] pairs or 'chrom.pos' strings.")

    return {
        "alleles": [str(allele_id) for allele_id in alleles],
        "missing": normalised_missing,
    }

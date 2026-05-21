from os import path

from Afanc.utilities.runCommands import command


def kraken_report_to_krona_text(kraken_report, krona_text):
    """Convert a Kraken-like report to ktImportText input.

    The input report is expected to contain the standard Kraken columns:
    percentage, clade reads, direct taxon reads, rank, taxID, and indented name.
    Krona receives only direct taxon reads to avoid double-counting ancestors.
    """
    rows_written = 0
    lineage_by_level = {}

    with open(kraken_report, "r") as fin, open(krona_text, "w") as fout:
        for line in fin:
            parsed = _parse_kraken_report_line(line)
            if parsed is None:
                continue

            level, taxon_reads, name = parsed
            if level == 0 and name == "root":
                lineage_by_level = {}
                continue

            for stale_level in list(lineage_by_level):
                if stale_level >= level:
                    del lineage_by_level[stale_level]
            lineage_by_level[level] = name
            if taxon_reads <= 0:
                continue

            lineage = [lineage_by_level[key] for key in sorted(lineage_by_level)]
            print("\t".join([str(taxon_reads), *lineage]), file=fout)
            rows_written += 1

    return rows_written


def run_krona_from_kraken_report(
    kraken_report,
    output_html,
    stdout=None,
    stderr=None,
    subprocess_id="KRONA",
):
    """Generate a Krona chart from a Kraken-like report using ktImportText."""
    krona_text = path.splitext(output_html)[0] + ".krona.txt"
    rows_written = kraken_report_to_krona_text(kraken_report, krona_text)
    if rows_written == 0:
        return {
            "status": "not_run",
            "reason": "no_nonzero_taxon_reads",
            "input": path.abspath(kraken_report),
            "krona_text": path.abspath(krona_text),
            "html": path.abspath(output_html),
            "rows": rows_written,
        }

    runline = ["ktImportText", krona_text, "-o", output_html]
    command(runline, subprocess_id).run_comm(0, stdout, stderr, raise_on_error=False)

    return {
        "status": "created" if path.exists(output_html) else "failed",
        "input": path.abspath(kraken_report),
        "krona_text": path.abspath(krona_text),
        "html": path.abspath(output_html),
        "rows": rows_written,
    }


def _parse_kraken_report_line(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 6:
        return None

    try:
        taxon_reads = int(round(float(fields[2])))
    except ValueError:
        return None

    raw_name = fields[-1]
    spaces = 0
    name = raw_name
    for char in raw_name:
        if char == " ":
            name = name[1:]
            spaces += 1
        else:
            break

    return int(spaces / 2), taxon_reads, name

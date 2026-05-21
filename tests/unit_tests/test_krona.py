from Afanc.utilities.krona import kraken_report_to_krona_text


def test_kraken_report_to_krona_text_uses_direct_reads_and_lineage_paths(tmp_path):
    report = tmp_path / "sample.k2.report.txt"
    report.write_text(
        "\n".join(
            [
                "100.00\t10\t0\tR\t1\troot",
                "100.00\t10\t2\tD\t2\t  Bacteria",
                "80.00\t8\t0\tG\t10\t    Genus",
                "80.00\t8\t8\tS\t20\t      Species",
            ]
        )
        + "\n"
    )
    krona_text = tmp_path / "sample.krona.txt"

    rows = kraken_report_to_krona_text(report, krona_text)

    assert rows == 2
    assert krona_text.read_text().splitlines() == [
        "2\tBacteria",
        "8\tBacteria\tGenus\tSpecies",
    ]


def test_kraken_report_to_krona_text_resets_same_level_sibling_branches(tmp_path):
    report = tmp_path / "mixed_myco.k2.report.txt"
    report.write_text(
        "\n".join(
            [
                "100.00\t10\t0\tR\t1\troot",
                "100.00\t10\t0\tD\t2\t  Bacteria",
                "100.00\t10\t1\tG\t1763\t    Mycobacterium",
                "60.00\t6\t0\tG1\t120793\t      Mycobacterium avium complex (MAC)",
                "30.00\t3\t3\tS\t1764\t        Mycobacterium avium",
                "30.00\t3\t3\tS\t1767\t        Mycobacterium intracellulare",
                "20.00\t2\t0\tG1\t2249310\t      Mycobacterium simiae complex",
                "20.00\t2\t2\tS\t33895\t        Mycobacterium interjectum",
                "20.00\t2\t0\tG1\t77643\t      Mycobacterium tuberculosis complex",
                "20.00\t2\t2\tS\t1773\t        Mycobacterium tuberculosis",
                "10.00\t1\t1\tS\t1768\t      Mycobacterium kansasii",
            ]
        )
        + "\n"
    )
    krona_text = tmp_path / "mixed_myco.krona.txt"

    rows = kraken_report_to_krona_text(report, krona_text)

    assert rows == 6
    assert krona_text.read_text().splitlines() == [
        "1\tBacteria\tMycobacterium",
        "3\tBacteria\tMycobacterium\tMycobacterium avium complex (MAC)\tMycobacterium avium",
        "3\tBacteria\tMycobacterium\tMycobacterium avium complex (MAC)\tMycobacterium intracellulare",
        "2\tBacteria\tMycobacterium\tMycobacterium simiae complex\tMycobacterium interjectum",
        "2\tBacteria\tMycobacterium\tMycobacterium tuberculosis complex\tMycobacterium tuberculosis",
        "1\tBacteria\tMycobacterium\tMycobacterium kansasii",
    ]

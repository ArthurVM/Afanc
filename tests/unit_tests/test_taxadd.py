import pandas as pd

from Afanc.autodatabase.taxadd import (
    getTaxidNames,
    searchTaxonAuthorityPrefix,
    taxonPrefix,
    taxonRank,
)


def _names_df():
    return pd.DataFrame(
        [
            [1763, "Mycobacterium", "", "scientific name"],
            [1774, "Mycobacterium chelonae Bergey et al. 1923 emend. Nouioui et al. 2018", "", "authority"],
            [1774, "Mycobacteroides chelonae", "", "scientific name"],
            [2038731, "Mycobacterium chelonae subsp. bovis Kim et al. 2017", "", "authority"],
            [2038731, "[Mycobacterium] chelonae subsp. bovis", "", "scientific name"],
        ]
    )


def _nodes_df():
    return pd.DataFrame(
        [
            [1763, 1762, "genus"],
            [1774, 670516, "species"],
            [2038731, 1774, "subspecies"],
        ]
    )


def test_taxon_prefix_and_rank_are_conservative():
    assert taxonPrefix("Mycobacterium chelonae") == "Mycobacterium chelonae"
    assert taxonRank("Mycobacterium chelonae") == "species"
    assert taxonPrefix("Mycobacterium chelonae subsp. bovis") == ""
    assert taxonRank("Mycobacterium chelonae subsp. bovis") is None


def test_authority_prefix_resolves_old_binomial_to_existing_species_taxid():
    taxid = searchTaxonAuthorityPrefix("Mycobacterium chelonae", _names_df(), _nodes_df())

    assert taxid == "1774"


def test_authority_prefix_does_not_resolve_below_species_labels():
    taxid = searchTaxonAuthorityPrefix("Mycobacterium tuberculosis Lineage 4.9", _names_df(), _nodes_df())

    assert taxid is None


def test_get_taxid_names_uses_authority_prefix_before_adding_taxon():
    taxadd_state = {"next_taxid": 3000000, "taxonomy_changed": False}

    taxid, names_df, nodes_df = getTaxidNames(
        "Mycobacterium chelonae",
        None,
        _names_df(),
        _nodes_df(),
        name_index={},
        taxadd_state=taxadd_state,
    )

    assert taxid == "1774"
    assert len(names_df) == len(_names_df())
    assert len(nodes_df) == len(_nodes_df())
    assert taxadd_state["next_taxid"] == 3000000
    assert taxadd_state["taxonomy_changed"] is False

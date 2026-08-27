"""The enzyme table and the digest report.

These are user-facing strings, and one of them was actively misleading before:
the site count was labelled "cuts", which invited the reader to conclude an
enzyme had worked when the number says nothing of the sort. Several tests here
exist to keep that from coming back.
"""
import pytest

from porec_tools.digest import DigestStats, find_cut_points
from porec_tools.enzymes import iter_usable, resolve_enzymes
from porec_tools.report import (
    describe_cut,
    describe_overhang,
    enzyme_table,
    expected_site_spacing,
    list_enzymes,
)


def one(name):
    return resolve_enzymes(name)[0]


@pytest.mark.parametrize("name,expected", [
    ("DpnII", "N^GATC_N"),
    ("NlaIII", "_CATG^"),
    ("HindIII", "A^AGCT_T"),
    ("HinfI", "G^ANT_C"),
    ("EcoRV", "GAT^_ATC"),
])
def test_cut_is_shown_on_both_strands(name, expected):
    assert describe_cut(one(name)) == expected


def test_cut_never_shows_fst5_to_the_user():
    """'fst5=0' is meaningless at the bench; the cut must be drawn."""
    for name in ("DpnII", "NlaIII", "HindIII", "HinfI", "MseI", "BsrI", "AcuI"):
        assert "fst5" not in describe_cut(one(name))
        assert "^" in describe_cut(one(name))


@pytest.mark.parametrize("name,expected", [
    ("DpnII", "5' GATC"),      # 5' overhang
    ("NlaIII", "3' CATG"),     # 3' overhang
    ("HindIII", "5' AGCT"),
    ("EcoRV", "blunt"),
    ("SmaI", "blunt"),
])
def test_overhang_reported_with_polarity(name, expected):
    assert describe_overhang(one(name)) == expected


def test_incompatible_ends_are_flagged():
    """DpnII leaves 5' GATC, NlaIII 3' CATG: they cannot ligate to each other.

    That is why junction motifs separate cleanly by enzyme, so it is worth
    saying out loud in the table.
    """
    text = "\n".join(enzyme_table(resolve_enzymes("DpnII,NlaIII")))
    assert "do not ligate to each other" in text


def test_no_compatibility_claim_for_a_single_enzyme():
    text = "\n".join(enzyme_table(resolve_enzymes("DpnII")))
    assert "ligate" not in text


def test_table_columns_stay_aligned():
    lines = enzyme_table(resolve_enzymes("DpnII,NlaIII,HindIII"))
    header, rows = lines[0], lines[1:4]
    start = header.index("site")
    for row in rows:
        assert row[start - 1] == " ", "column must not run into the previous"


def test_expected_spacing_counts_ambiguity():
    assert expected_site_spacing("GATC") == 256
    assert expected_site_spacing("AAGCTT") == 4096
    assert expected_site_spacing("GANTC") == 256     # N contributes nothing
    assert expected_site_spacing("CCWGG") == 512     # W allows 2 of 4
    assert expected_site_spacing("GGYRCC") == 1024   # two 2-fold codes


def _report(spec, seq):
    stats = DigestStats()
    stats.n_concatemers, stats.n_bases = 1, len(seq)
    enzymes = resolve_enzymes(spec)
    find_cut_points(seq, enzymes, stats)
    return "\n".join(stats.summary_lines(enzymes))


def test_report_does_not_say_cuts():
    """Regression guard on the flaw this module was written to fix."""
    text = _report("DpnII,HindIII", "GATC" + "TTTTAAGCTTTT" * 20)
    assert "sites found" in text
    assert "cuts (" not in text
    assert "not proof an enzyme cut" in text


def test_report_points_at_the_test_that_works():
    text = _report("DpnII,NlaIII", "GATCCATG" * 20)
    assert "ALIGNED" in text
    assert "INTEGRATION.md" in text


def test_report_warns_only_when_a_site_is_truly_absent():
    absent = _report("DpnII,HindIII", "GATC" * 40)          # no AAGCTT at all
    assert "WARNING" in absent and "HindIII" in absent
    present = _report("DpnII,HindIII", ("GATC" + "AAGCTT") * 20)
    assert "WARNING" not in present


def test_report_survives_an_empty_digest():
    stats = DigestStats()
    lines = stats.summary_lines(resolve_enzymes("DpnII"))
    assert lines and "0 concatemers" in lines[0]


# --------------------------------------------------------------------------
# the enzyme catalogue
# --------------------------------------------------------------------------
def test_catalogue_shows_the_3c_workhorses():
    text = "\n".join(list_enzymes())
    for name in ("DpnII", "NlaIII", "HinfI", "MseI", "HindIII"):
        assert name in text


def test_catalogue_says_it_is_not_a_restriction():
    """A short list must not read as "these are the only ones allowed"."""
    text = "\n".join(list_enzymes())
    assert "not a restriction" in text
    assert "--all" in text


def test_catalogue_lists_every_usable_enzyme():
    assert len(list(iter_usable())) > 700


def test_catalogue_excludes_what_resolve_would_reject():
    """Listing an enzyme the digest refuses would only invite a bad error."""
    names = {e.name for e in iter_usable()}
    assert "BcgI" not in names, "cuts twice"
    assert "ScoDS2II" not in names, "no defined cut position"
    for name in names:
        resolve_enzymes(name)          # must not raise


def test_search_by_name_fragment():
    assert {e.name for e in iter_usable("Dpn")} == {"DpnI", "DpnII"}


def test_search_by_recognition_site():
    found = {e.name for e in iter_usable("CATG")}
    assert "NlaIII" in found and "FatI" in found
    assert all(e.site == "CATG" for e in iter_usable("CATG"))


def test_search_is_case_insensitive():
    assert {e.name for e in iter_usable("dpn")} == \
           {e.name for e in iter_usable("DPN")}


def test_search_with_no_match_explains_what_to_try():
    text = "\n".join(list_enzymes("ZZZZ"))
    assert "No usable enzyme matches" in text
    assert "GATC" in text and "--all" in text


def test_no_ligation_note_in_a_catalogue():
    """DpnI and DpnII are alternatives here, not a combination.

    The compatibility warning belongs to enzymes the user chose to digest
    with; in a listing it is a non-sequitur.
    """
    assert "do not ligate" not in "\n".join(list_enzymes("Dpn"))


def test_ligation_note_still_appears_for_a_chosen_set():
    text = "\n".join(enzyme_table(resolve_enzymes("DpnII,NlaIII")))
    assert "do not ligate" in text

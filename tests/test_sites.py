"""Site location: overlaps, IUPAC codes, boundaries, version independence."""
import random

import pytest
from Bio import Restriction
from Bio.Seq import Seq

from pore_c_aqb.enzymes import _resolve_one
from pore_c_aqb.sites import IUPAC_CLASSES, find_cuts_for_enzyme, site_regex

from conftest import random_seq

#: every enzyme in this list is palindromic and used in real 3C/Hi-C protocols
PALINDROMIC = [
    "DpnII", "NlaIII", "HinfI", "MseI", "HindIII", "EcoRI", "BglII", "NcoI",
    "Csp6I", "ApoI", "DdeI", "Sau3AI", "MboI", "BamHI", "XhoI", "AluI",
    "HaeIII", "TaqI", "NotI", "SmaI",
]


def _site_touches_edge(pos, enzyme, length):
    """True when the site producing `pos` runs into either end of the read."""
    start = pos - enzyme.fst5
    return start <= 0 or start + len(enzyme.site) >= length


@pytest.mark.parametrize("name", PALINDROMIC)
def test_agrees_with_biopython_away_from_edges(name):
    """The contract: identical to Biopython except at sequence boundaries.

    Boundary sites are the one documented divergence (Biopython changed its
    own answer there between 1.82 and 1.88); everywhere else we must match
    exactly, whatever Biopython version is installed.
    """
    enzyme = _resolve_one(name)
    rng = random.Random(hash(name) % 10_000)
    for _ in range(60):
        length = rng.randint(30, 2000)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        reference = {x - 1 for x in getattr(Restriction, name).search(Seq(seq))}
        got = find_cuts_for_enzyme(seq, enzyme)
        for pos in (reference ^ got):
            assert _site_touches_edge(pos, enzyme, length), (
                f"{name}: disagreement at {pos} is not a boundary case"
            )


def test_overlapping_sites_are_all_reported():
    """A lookahead scan must not skip past an overlapping match."""
    # TTAA overlaps itself in TTAATTAA only at distinct offsets, but AluI's
    # AGCT overlaps in AGCTAGCT at 0 and 4; use a self-overlapping case:
    seq = "AAA" + "CATGCATG" + "AAA"
    enzyme = _resolve_one("NlaIII")
    cuts = find_cuts_for_enzyme(seq, enzyme)
    assert len(cuts) == 2


def test_site_at_position_zero_is_reported():
    """Biopython 1.88 drops this; we keep it, matching 1.82 and the pipeline."""
    enzyme = _resolve_one("NlaIII")          # CATG^, fst5 = 4
    cuts = find_cuts_for_enzyme("CATG" + "TTTT" * 5, enzyme)
    assert 4 in cuts


def test_site_at_the_very_end_is_reported():
    enzyme = _resolve_one("DpnII")           # ^GATC, fst5 = 0
    seq = "TTTT" * 5 + "GATC"
    cuts = find_cuts_for_enzyme(seq, enzyme)
    assert (len(seq) - 4) in cuts


def test_cuts_never_fall_outside_the_read():
    enzyme = _resolve_one("NlaIII")
    for seq in ("CATG", "C", "", "CATGCATG"):
        for pos in find_cuts_for_enzyme(seq, enzyme):
            assert 0 <= pos <= len(seq)


def test_lowercase_sequence_is_handled():
    enzyme = _resolve_one("DpnII")
    upper = find_cuts_for_enzyme("TTTTGATCTTTT", enzyme)
    lower = find_cuts_for_enzyme("ttttgatctttt", enzyme)
    assert upper == lower and len(upper) == 1


def test_non_acgt_characters_never_match():
    enzyme = _resolve_one("DpnII")
    assert find_cuts_for_enzyme("NNNNNNNNNNNN", enzyme) == set()


def test_ambiguous_site_matches_every_variant():
    """HinfI is GANTC: all four middle bases must be recognised."""
    enzyme = _resolve_one("HinfI")
    for middle in "ACGT":
        seq = "TTTT" + f"GA{middle}TC" + "TTTT"
        assert len(find_cuts_for_enzyme(seq, enzyme)) == 1, middle


def test_regex_is_cached():
    assert site_regex("GATC") is site_regex("GATC")


def test_iupac_table_covers_every_usable_enzyme():
    """Every enzyme we accept must have a translatable recognition site.

    Enzymes we reject (no defined cut position, two alternative sites, cuts
    twice) are excluded: they never reach the scanner. HpyUM037X for instance
    has site 'TNGGNAG|GTGGNAG' and fst5 = None.
    """
    from pore_c_aqb.enzymes import EnzymeSpecError, _resolve_one
    used, n_ok = set(), 0
    for name in dir(Restriction):
        if name.startswith("_"):
            continue
        try:
            enzyme = _resolve_one(name)
        except (EnzymeSpecError, TypeError, AttributeError):
            continue
        used |= set(enzyme.site.upper())
        n_ok += 1
    assert n_ok > 500, "expected Biopython to expose many usable enzymes"
    assert used <= set(IUPAC_CLASSES), \
        f"missing codes: {used - set(IUPAC_CLASSES)}"


def test_enzyme_with_alternative_sites_is_rejected():
    """HpyUM037X has two recognition sequences and no single cut position."""
    from pore_c_aqb.enzymes import EnzymeSpecError, _resolve_one
    with pytest.raises(EnzymeSpecError):
        _resolve_one("HpyUM037X")


def test_non_palindromic_delegates_to_biopython():
    """Type IIS enzymes are handed back to Biopython, and must still work."""
    enzyme = _resolve_one("BsrI")
    seq = random_seq(1500, seed=3)
    expected = {x - 1 for x in Restriction.BsrI.search(Seq(seq))}
    assert find_cuts_for_enzyme(seq, enzyme) == expected

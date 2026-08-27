"""Cut-point union, interval splitting and monomer emission."""
import random

import pysam
import pytest
from Bio import Restriction
from Bio.Seq import Seq

from porec_tools.digest import DigestStats, digest_sequence, find_cut_points
from porec_tools.enzymes import resolve_enzymes

from conftest import make_read, random_seq


def _biopython_cuts(seq, names):
    cuts = set()
    for n in names:
        cuts |= {x - 1 for x in getattr(Restriction, n).search(Seq(seq))}
    return sorted(cuts)


def test_single_enzyme_matches_biopython(spaced_sites):
    got = find_cut_points(spaced_sites, resolve_enzymes("DpnII"))
    assert got == _biopython_cuts(spaced_sites, ["DpnII"])


def test_union_is_the_union(spaced_sites):
    enzymes = resolve_enzymes("DpnII,NlaIII")
    got = find_cut_points(spaced_sites, enzymes)
    assert got == _biopython_cuts(spaced_sites, ["DpnII", "NlaIII"])


def test_union_is_superset_of_each_part(spaced_sites):
    dpn = set(find_cut_points(spaced_sites, resolve_enzymes("DpnII")))
    nla = set(find_cut_points(spaced_sites, resolve_enzymes("NlaIII")))
    both = set(find_cut_points(spaced_sites, resolve_enzymes("DpnII,NlaIII")))
    assert both == dpn | nla
    assert both >= dpn and both >= nla


def test_cut_points_sorted_and_unique(spaced_sites):
    cuts = find_cut_points(spaced_sites, resolve_enzymes("DpnII,NlaIII,HinfI"))
    assert cuts == sorted(set(cuts))


def test_overlapping_sites_are_both_found():
    """CATGATC contains CATG at 0 and GATC at 3: a naive scan misses one.

    This is exactly the failure mode of a single regex alternation pass,
    which is why the implementation searches per enzyme.
    """
    seq = "TTTT" + "CATGATC" + "TTTT"
    cuts = find_cut_points(seq, resolve_enzymes("DpnII,NlaIII"))
    assert cuts == _biopython_cuts(seq, ["DpnII", "NlaIII"])
    assert len(cuts) == 2, "both the CATG and the GATC cut must be present"


def test_no_sites_gives_one_monomer():
    seq = "AAAACCCCTTTT" * 5
    reads = list(digest_sequence(
        make_read("r", seq), resolve_enzymes("DpnII")))
    assert len(reads) == 1
    assert reads[0].query_sequence == seq


def test_monomers_tile_the_read_without_gap_or_overlap(spaced_sites):
    reads = list(digest_sequence(
        make_read("r", spaced_sites), resolve_enzymes("DpnII,NlaIII")))
    spans = [tuple(r.get_tag("Xc")[:2]) for r in reads]
    assert spans[0][0] == 0
    assert spans[-1][1] == len(spaced_sites)
    for (_, end), (start, _) in zip(spans, spans[1:]):
        assert end == start, "monomers must tile the read exactly"
    assert "".join(r.query_sequence for r in reads) == spaced_sites


def test_more_enzymes_never_fewer_monomers(spaced_sites):
    n1 = len(list(digest_sequence(
        make_read("r", spaced_sites), resolve_enzymes("DpnII"))))
    n2 = len(list(digest_sequence(
        make_read("r", spaced_sites), resolve_enzymes("DpnII,NlaIII"))))
    assert n2 >= n1


def test_monomer_names_are_lexicographically_sortable(spaced_sites):
    reads = list(digest_sequence(
        make_read("read1", spaced_sites), resolve_enzymes("DpnII,NlaIII")))
    names = [r.query_name for r in reads]
    assert names == sorted(names)
    for n in names:
        assert n.startswith("read1:")


def test_concatemer_tag_preserved(spaced_sites):
    reads = list(digest_sequence(
        make_read("read1", spaced_sites), resolve_enzymes("DpnII,NlaIII")))
    assert all(r.get_tag("MI") == "read1" for r in reads)


def test_qualities_are_trimmed_with_the_sequence(spaced_sites):
    qual = "".join(chr(33 + (i % 40)) for i in range(len(spaced_sites)))
    reads = list(digest_sequence(
        make_read("r", spaced_sites, qual), resolve_enzymes("DpnII,NlaIII")))
    for r in reads:
        assert len(r.query_qualities) == len(r.query_sequence)
    rebuilt = "".join(
        pysam.qualities_to_qualitystring(r.query_qualities) for r in reads)
    assert rebuilt == qual


def test_stats_count_each_enzyme(spaced_sites):
    stats = DigestStats()
    enzymes = resolve_enzymes("DpnII,NlaIII")
    find_cut_points(spaced_sites, enzymes, stats)
    assert stats.sites_per_enzyme["DpnII"] >= 1
    assert stats.sites_per_enzyme["NlaIII"] >= 1
    assert stats.n_cut_points >= 2


def test_stats_flag_an_enzyme_with_no_site_at_all():
    """Zero sites means the enzyme name is wrong or the data is not this run."""
    seq = "GATC".join(["TTTT"] * 20)      # DpnII sites only, no AAGCTT
    stats = DigestStats()
    stats.n_concatemers = 1
    stats.n_bases = len(seq)
    enzymes = resolve_enzymes("DpnII,HindIII")
    find_cut_points(seq, enzymes, stats)
    assert stats.sites_per_enzyme["DpnII"] > 0
    assert stats.sites_per_enzyme.get("HindIII", 0) == 0
    report = "\n".join(stats.summary_lines(enzymes))
    assert "WARNING" in report and "HindIII" in report


def test_report_never_claims_an_enzyme_cut():
    """The count is sites present in the reads, not proof of cutting.

    Overclaiming here would be worse than saying nothing: a user reading
    "1,497 cuts" for an enzyme that failed draws exactly the wrong conclusion.
    """
    seq = "GATC".join(["TTTTAAGCTTTTTT"] * 20)   # both motifs present
    stats = DigestStats()
    stats.n_concatemers = 1
    stats.n_bases = len(seq)
    enzymes = resolve_enzymes("DpnII,HindIII")
    find_cut_points(seq, enzymes, stats)
    report = "\n".join(stats.summary_lines(enzymes))
    assert "sites found" in report
    assert "cuts (" not in report, "must not label site counts as cuts"
    assert "not proof an enzyme cut" in report
    assert "means nothing on its own" in report
    assert "ALIGNED" in report, "must point at the test that actually works"


def test_shared_cut_point_counted_once():
    """Two enzymes cutting the same base is one cut in the tube."""
    stats = DigestStats()
    # DpnII and MboI are isoschizomers; resolve_enzymes collapses them, so
    # build the degenerate case by hand to exercise the accounting.
    from porec_tools.enzymes import _resolve_one
    enzymes = [_resolve_one("DpnII"), _resolve_one("MboI")]
    seq = "TTTT" + "GATC" + "TTTT"
    cuts = find_cut_points(seq, enzymes, stats)
    assert len(cuts) == 1
    assert stats.n_shared_cut_points == 1


@pytest.mark.parametrize("seed", range(12))
def test_random_sequences_match_biopython(seed):
    """Property test: the union must equal Biopython's, always."""
    seq = random_seq(random.Random(seed).randint(60, 4000), seed=seed)
    names = ["DpnII", "NlaIII", "HinfI", "MseI"]
    got = find_cut_points(seq, resolve_enzymes(",".join(names)))
    assert got == _biopython_cuts(seq, names)


@pytest.mark.parametrize("seed", range(8))
def test_random_sequences_tile_exactly(seed):
    seq = random_seq(random.Random(seed).randint(60, 3000), seed=seed + 100)
    reads = list(digest_sequence(
        make_read("r", seq), resolve_enzymes("DpnII,NlaIII")))
    assert "".join(r.query_sequence for r in reads) == seq


def test_expected_spacing_matches_site_complexity():
    from porec_tools.report import expected_site_spacing
    assert expected_site_spacing("GATC") == 256        # 4^4
    assert expected_site_spacing("AAGCTT") == 4096     # 4^6
    assert expected_site_spacing("GANTC") == 256       # N is free
    assert expected_site_spacing("CCWGG") == 512       # W allows 2 of 4


def test_observed_spacing_is_reported(spaced_sites):
    stats = DigestStats()
    stats.n_concatemers = 1
    stats.n_bases = len(spaced_sites)
    enzymes = resolve_enzymes("DpnII,NlaIII")
    find_cut_points(spaced_sites, enzymes, stats)
    report = "\n".join(stats.summary_lines(enzymes))
    assert "1 per (observed)" in report and "1 per (chance)" in report

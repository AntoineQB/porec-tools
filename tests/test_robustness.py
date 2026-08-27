"""Awkward inputs: empty files, wrong formats, odd sequences, missing paths.

A tool that other people will run gets fed things its author never imagined.
Every case here was found by deliberately trying to break the commands, and
several of them did break something before these tests existed.

The rule: fail with a sentence the user can act on, or handle it correctly.
Never a bare traceback.
"""
import array
import io
import contextlib

import pysam
import pytest

from pore_c_aqb.cli import main as cli
from pore_c_aqb.digest import digest_sequence, find_cut_points
from pore_c_aqb.enzymes import resolve_enzymes
from pore_c_aqb.junctions import main as junctions_main
from pore_c_aqb.merge import main as merge_main
from pore_c_aqb.reads import iter_concatemers, read_span, UNKNOWN_SPAN

from conftest import make_read, write_bam

HEADER = {"HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
          "SQ": [{"SN": "chr1", "LN": 300000}, {"SN": "chr2", "LN": 300000}]}


def aligned(name, chrom, start, length=200, mapq=60, rev=False,
            mi=None, xc=None, unmapped=False):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = "ACGT" * (length // 4)
    a.query_qualities = pysam.qualitystring_to_array("I" * (length // 4 * 4))
    if unmapped:
        a.flag = 4
    else:
        a.flag = 16 if rev else 0
        a.reference_id = 0 if chrom == "chr1" else 1
        a.reference_start = start
        a.cigar = [(0, length // 4 * 4)]
        a.mapping_quality = mapq
    if mi:
        a.set_tag("MI", mi)
    if xc is not None:
        a.set_tag("Xc", array.array("I", xc))
    return a


def aligned_bam(path, records):
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
        for r in records:
            out.write(r)
    return str(path)


def quiet(fn):
    """Run something, swallowing its chatter, and give back the return code."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(), buf.getvalue()


# --------------------------------------------------------------------------
# odd sequences
# --------------------------------------------------------------------------
@pytest.mark.parametrize("seq,expected", [
    ("ttttgatcttttcatgtttt", [4, 16]),      # lowercase input
    ("TTttGaTcTTTTcAtGTTTT", [4, 16]),      # mixed case
    ("TTTTNNNNGATCNNNNTTTT", [8]),          # N must not match a literal base
    ("", []),                               # empty
    ("GA", []),                             # shorter than any site
    ("GATCTTTT", [0]),                      # site flush against the start
    ("TTTTGATC", [4]),                      # site flush against the end
    ("GATCGATCGATC", [0, 4, 8]),            # nothing but sites
])
def test_site_search_survives_odd_sequences(seq, expected):
    assert find_cut_points(seq, resolve_enzymes("DpnII,NlaIII")) == expected


@pytest.mark.parametrize("seq", [
    "ttttgatctttt", "A", "NNNNGATCNNNN", "AAAACCCCTTTT", "GATC",
])
def test_digest_reconstructs_any_read(seq):
    """Monomers must always tile the read exactly, whatever it contains.

    Compare against ``read.query_sequence``, not the input string: a BAM
    encodes DNA on 4 bits, so pysam hands sequences back uppercase and case is
    simply not a property a BAM can carry.
    """
    read = make_read("r", seq)
    stored = read.query_sequence
    monomers = list(digest_sequence(read, resolve_enzymes("DpnII,NlaIII")))
    assert "".join(m.query_sequence for m in monomers) == stored


def test_digest_end_to_end_on_awkward_reads(tmp_path):
    reads = [make_read("lower", "ttttgatcttttcatgtttt"),
             make_read("withN", "TTTTNNNNGATCTTTT"),
             make_read("tiny", "A"),
             make_read("nosite", "AAAACCCC")]
    originals = {r.query_name: r.query_sequence for r in reads}
    src = write_bam(tmp_path / "edge.bam", reads)
    out = tmp_path / "out.bam"
    assert cli(["digest", "DpnII,NlaIII", src,
                "--output", str(out), "--quiet"]) == 0
    got = {}
    with pysam.AlignmentFile(str(out), check_sq=False) as fh:
        for m in fh:
            got.setdefault(m.get_tag("MI"), []).append(m.query_sequence)
    assert set(got) == set(originals)
    for name, pieces in got.items():
        assert "".join(pieces) == originals[name]


# --------------------------------------------------------------------------
# missing, empty and malformed files
# --------------------------------------------------------------------------
def test_merge_missing_bam(tmp_path):
    with pytest.raises(SystemExit, match="Input file not found"):
        merge_main([str(tmp_path / "nope.bam"), "--output", str(tmp_path / "o")])


def test_merge_file_that_is_not_a_bam(tmp_path):
    junk = tmp_path / "junk.bam"
    junk.write_text("this is not a BAM")
    with pytest.raises(SystemExit, match="Cannot read"):
        merge_main([str(junk), "--output", str(tmp_path / "o")])


def test_merge_missing_sizes_file(tmp_path):
    bam = aligned_bam(tmp_path / "a.bam",
                      [aligned("r:0000:0200", "chr1", 1000, mi="r",
                               xc=[0, 200, 0, 0, 1])])
    with pytest.raises(SystemExit, match="sizes file not found"):
        merge_main([bam, "--output", str(tmp_path / "o"),
                    "--pairs", str(tmp_path / "p"),
                    "--sizes", str(tmp_path / "nope.sizes")])


def test_merge_empty_bam(tmp_path):
    bam = aligned_bam(tmp_path / "empty.bam", [])
    out = tmp_path / "o.tsv"
    rc, _ = quiet(lambda: merge_main([bam, "--output", str(out)]))
    assert rc == 0
    assert out.read_text().strip().count("\n") == 0, "header only"


def test_merge_only_unmapped_reads(tmp_path):
    bam = aligned_bam(tmp_path / "u.bam",
                      [aligned("r:0000:0200", None, 0, mi="r", unmapped=True)])
    out = tmp_path / "o.tsv"
    rc, _ = quiet(lambda: merge_main([bam, "--output", str(out)]))
    assert rc == 0


def test_merge_single_monomer_molecule(tmp_path):
    bam = aligned_bam(tmp_path / "s.bam",
                      [aligned("r:0000:0200", "chr1", 1000, mi="r",
                               xc=[0, 200, 0, 0, 1])])
    out = tmp_path / "o.tsv"
    rc, _ = quiet(lambda: merge_main([bam, "--output", str(out)]))
    assert rc == 0
    assert len(out.read_text().strip().split("\n")) == 2, "header + one row"


def test_digest_missing_input(tmp_path):
    with pytest.raises(SystemExit, match="Input file not found"):
        cli(["digest", "DpnII", str(tmp_path / "nope.bam"),
             "--output", str(tmp_path / "o.bam")])


def test_junctions_missing_reference(tmp_path):
    bam = aligned_bam(tmp_path / "j.bam",
                      [aligned("r:0000:0200", "chr1", 1000, mi="r",
                               xc=[0, 200, 0, 0, 1])])
    rc, text = quiet(lambda: junctions_main(
        [bam, str(tmp_path / "nope.fa"), "--enzymes", "DpnII"]))
    assert rc == 2
    assert "faidx" in text or "cannot open" in text.lower()


# --------------------------------------------------------------------------
# BAMs that did not come from this tool
# --------------------------------------------------------------------------
def test_read_span_unknown_when_nothing_says_where(tmp_path):
    a = aligned("noPositionInfo", "chr1", 1000)
    assert read_span(a) == UNKNOWN_SPAN


def test_read_span_survives_a_non_numeric_suffix():
    a = aligned("some:weird:name", "chr1", 1000)
    assert read_span(a) == UNKNOWN_SPAN


def test_merge_falls_back_to_file_order(tmp_path):
    """No Xc and unparseable names: keep the order records came in, and warn."""
    bam = aligned_bam(tmp_path / "w.bam",
                      [aligned("weirdA", "chr1", 1000),
                       aligned("weirdB", "chr2", 5000)])
    out = tmp_path / "o.tsv"
    rc, _ = quiet(lambda: merge_main([bam, "--output", str(out)]))
    assert rc == 0
    rows = out.read_text().strip().split("\n")[1:]
    assert len(rows) == 2
    assert [r.split("\t")[3] for r in rows] == ["chr1", "chr2"]
    assert all(r.split("\t")[10] == "-1" for r in rows), \
        "read positions are unknown and must be reported as such, not invented"


def test_unknown_read_positions_are_warned_about(tmp_path, caplog):
    bam = aligned_bam(tmp_path / "w.bam", [aligned("weirdA", "chr1", 1000)])
    with caplog.at_level("WARNING"):
        list(iter_concatemers(bam))
    assert any("unknown" in r.message.lower() or "Xc" in r.message
               for r in caplog.records)


# --------------------------------------------------------------------------
# parameter extremes
# --------------------------------------------------------------------------
@pytest.mark.parametrize("gap", ["-5", "0", "999999999"])
def test_merge_gap_extremes_do_not_crash(tmp_path, gap):
    bam = aligned_bam(tmp_path / f"g{gap}.bam", [
        aligned("r:0000:0200", "chr1", 1000, mi="r", xc=[0, 200, 0, 0, 2]),
        aligned("r:0200:0400", "chr1", 1200, mi="r", xc=[200, 400, 0, 1, 2])])
    out = tmp_path / f"o{gap}.tsv"
    rc, _ = quiet(lambda: merge_main(
        [bam, "--output", str(out), "--merge-gap", gap]))
    assert rc == 0


def test_a_huge_gap_never_merges_across_chromosomes(tmp_path):
    """The gap is a distance, and distance is meaningless between chromosomes."""
    bam = aligned_bam(tmp_path / "x.bam", [
        aligned("r:0000:0200", "chr1", 1000, mi="r", xc=[0, 200, 0, 0, 2]),
        aligned("r:0200:0400", "chr2", 1000, mi="r", xc=[200, 400, 0, 1, 2])])
    out = tmp_path / "o.tsv"
    quiet(lambda: merge_main(
        [bam, "--output", str(out), "--merge-gap", "999999999"]))
    assert len(out.read_text().strip().split("\n")) == 3, "header + two rows"


@pytest.mark.parametrize("mapq", ["0", "255"])
def test_mapq_extremes(tmp_path, mapq):
    bam = aligned_bam(tmp_path / f"q{mapq}.bam",
                      [aligned("r:0000:0200", "chr1", 1000, mi="r",
                               xc=[0, 200, 0, 0, 1])])
    out = tmp_path / f"o{mapq}.tsv"
    rc, _ = quiet(lambda: merge_main(
        [bam, "--output", str(out), "--mapq", mapq]))
    assert rc == 0


def test_pairs_without_sizes_still_produces_a_valid_file(tmp_path):
    bam = aligned_bam(tmp_path / "p.bam", [
        aligned("r:0000:0200", "chr1", 1000, mi="r", xc=[0, 200, 0, 0, 2]),
        aligned("r:0200:0400", "chr2", 9000, mi="r", xc=[200, 400, 0, 1, 2])])
    pairs = tmp_path / "c.pairs"
    quiet(lambda: merge_main([bam, "--output", str(tmp_path / "o.tsv"),
                              "--pairs", str(pairs)]))
    lines = pairs.read_text().splitlines()
    assert lines[0] == "## pairs format v1.0"
    assert len([ln for ln in lines if not ln.startswith("#")]) == 1

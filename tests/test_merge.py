"""Merging fragments the in silico digest split but the enzyme never cut.

The bug this guards against is silent and expensive: a Hi-C map whose diagonal
is inflated by contacts that never happened. Several tests below encode
mistakes actually made while developing this, so they cannot come back.
"""
import array
import gzip
import json

import pysam
import pytest

from pore_c_aqb.merge import Fragment, main, merge_adjacent
from pore_c_aqb.reads import iter_concatemers

CHROMS = {"chr1": 300000, "chr2": 300000}


def frag(read_start, chrom, ref_start, length=200, strand="+", mapq=60):
    return Fragment(read_start, read_start + length, chrom,
                    ref_start, ref_start + length, strand, mapq)


# --------------------------------------------------------------------------
# the merging rule itself
# --------------------------------------------------------------------------
def test_contiguous_monomers_are_merged():
    """Two pieces that touch on the genome were one uncut fragment."""
    got = merge_adjacent([frag(0, "chr1", 1000), frag(200, "chr1", 1200)], 100)
    assert len(got) == 1
    assert (got[0].ref_start, got[0].ref_end) == (1000, 1400)
    assert got[0].n_monomers == 2


def test_distant_monomers_are_a_real_contact():
    got = merge_adjacent([frag(0, "chr1", 1000), frag(200, "chr1", 900000)], 100)
    assert len(got) == 2


def test_different_chromosomes_never_merge():
    got = merge_adjacent([frag(0, "chr1", 1000), frag(200, "chr2", 1200)], 100)
    assert len(got) == 2


def test_different_strands_never_merge():
    """Opposite orientation means a ligation, not an uncut site."""
    got = merge_adjacent(
        [frag(0, "chr1", 1000), frag(200, "chr1", 1200, strand="-")], 100)
    assert len(got) == 2


def test_gap_threshold_is_respected():
    close = merge_adjacent([frag(0, "chr1", 1000), frag(200, "chr1", 1250)], 100)
    far = merge_adjacent([frag(0, "chr1", 1000), frag(200, "chr1", 1500)], 100)
    assert len(close) == 1 and len(far) == 2


def test_merging_works_backwards_on_the_minus_strand():
    """A minus-strand fragment runs backwards along the reference.

    The next monomer's start then *precedes* the previous one's start, so the
    gap has to be measured in both directions or the merge silently misses
    every reverse-aligned pair.
    """
    got = merge_adjacent(
        [frag(0, "chr1", 1200, strand="-"), frag(200, "chr1", 1000, strand="-")],
        100)
    assert len(got) == 1
    assert (got[0].ref_start, got[0].ref_end) == (1000, 1400)


def test_a_whole_chain_collapses_to_one():
    chain = [frag(200 * i, "chr1", 1000 + 200 * i) for i in range(6)]
    got = merge_adjacent(chain, 100)
    assert len(got) == 1 and got[0].n_monomers == 6


def test_merged_block_keeps_the_best_mapq():
    """A block containing one confident monomer must survive the MAPQ filter."""
    got = merge_adjacent(
        [frag(0, "chr1", 1000, mapq=0), frag(200, "chr1", 1200, mapq=60)], 100)
    assert got[0].mapq == 60


def test_read_order_is_what_counts():
    """Monomers are merged along the read, whatever order they arrive in."""
    a, b = frag(0, "chr1", 1000), frag(200, "chr1", 1200)
    assert len(merge_adjacent([b, a], 100)) == 1


def test_empty_input():
    assert merge_adjacent([], 100) == []


def test_merging_never_increases_the_count():
    chain = [frag(200 * i, "chr1", 1000 + 500 * i) for i in range(10)]
    for gap in (0, 10, 100, 1000, 10000):
        assert len(merge_adjacent(chain, gap)) <= len(chain)


def test_larger_gap_never_yields_more_fragments():
    """Monotonicity: loosening the threshold can only merge more."""
    chain = [frag(200 * i, "chr1", 1000 + 300 * i) for i in range(8)]
    counts = [len(merge_adjacent(chain, g)) for g in (0, 50, 100, 500, 1000)]
    assert counts == sorted(counts, reverse=True)


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------
def write_bam(path, records):
    """records: [(read_id, [(chrom, start, end, rev, mapq), ...]), ...]"""
    header = {"HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
              "SQ": [{"SN": n, "LN": ln} for n, ln in CHROMS.items()]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for read_id, monomers in records:
            offset = 0
            for i, (chrom, start, end, rev, mapq) in enumerate(monomers):
                n = end - start
                a = pysam.AlignedSegment(out.header)
                a.query_name = f"{read_id}:{offset:04d}:{offset + n:04d}"
                a.query_sequence = "A" * n
                a.query_qualities = pysam.qualitystring_to_array("I" * n)
                a.reference_id = out.header.get_tid(chrom)
                a.reference_start = start
                a.cigar = [(0, n)]
                a.mapping_quality = mapq
                a.flag = 16 if rev else 0
                a.set_tag("MI", read_id)
                a.set_tag("Xc", array.array(
                    "I", [offset, offset + n, 0, i, len(monomers)]))
                out.write(a)
                offset += n
    return str(path)


@pytest.fixture
def split_locus_bam(tmp_path):
    """A molecule whose first locus was cut into 4 pieces in silico only.

    Real structure: one 800 bp fragment on chr1 in contact with chr2. The
    virtual digest split the chr1 fragment into 4 monomers, so the naive
    reading sees 5 loci and manufactures 10 pairs instead of 1.
    """
    monomers = [("chr1", 1000 + 200 * i, 1200 + 200 * i, False, 60)
                for i in range(4)]
    monomers.append(("chr2", 50000, 50200, False, 60))
    return write_bam(tmp_path / "split.bam", [("readA", monomers)])


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]


def test_end_to_end_collapses_the_split_locus(tmp_path, split_locus_bam):
    out = tmp_path / "frags.tsv"
    assert main([split_locus_bam, "--output", str(out)]) == 0
    rows = read_tsv(out)
    assert len(rows) == 2, "4 chr1 pieces + 1 chr2 piece must become 2 loci"
    chr1 = next(r for r in rows if r["chrom"] == "chr1")
    assert chr1["n_monomers_merged"] == "4"
    assert (int(chr1["start"]), int(chr1["end"])) == (1000, 1800)


def test_pairs_output_has_one_contact_not_ten(tmp_path, split_locus_bam):
    """The whole point: 5 monomers give 10 pairs, 2 fragments give 1."""
    out, pairs = tmp_path / "f.tsv", tmp_path / "c.pairs"
    main([split_locus_bam, "--output", str(out), "--pairs", str(pairs)])
    body = [ln for ln in pairs.read_text().splitlines()
            if not ln.startswith("#")]
    assert len(body) == 1
    fields = body[0].split("\t")
    assert fields[1] == "chr1" and fields[3] == "chr2"


def test_pairs_header_is_valid_4dn(tmp_path, split_locus_bam):
    out, pairs = tmp_path / "f.tsv", tmp_path / "c.pairs"
    sizes = tmp_path / "sizes.genome"
    sizes.write_text("".join(f"{c}\t{ln}\n" for c, ln in CHROMS.items()))
    main([split_locus_bam, "--output", str(out), "--pairs", str(pairs),
          "--sizes", str(sizes)])
    lines = pairs.read_text().splitlines()
    assert lines[0] == "## pairs format v1.0"
    assert "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2" in lines
    assert "#chromsize: chr1 300000" in lines


def test_pairs_positions_are_one_based(tmp_path, split_locus_bam):
    out, pairs = tmp_path / "f.tsv", tmp_path / "c.pairs"
    main([split_locus_bam, "--output", str(out), "--pairs", str(pairs)])
    row = [ln for ln in pairs.read_text().splitlines()
           if not ln.startswith("#")][0].split("\t")
    rows = read_tsv(out)
    mid = {r["chrom"]: int(r["midpoint"]) for r in rows}
    assert int(row[2]) == mid["chr1"] + 1
    assert int(row[4]) == mid["chr2"] + 1


# --------------------------------------------------------------------------
# the ordering trap
# --------------------------------------------------------------------------
def test_mapq_filter_runs_after_merging(tmp_path):
    """Filtering first breaks the chain and invents a contact.

    Four contiguous chr1 monomers where the second maps poorly. Filter first
    and the chain splits into two blocks that look like two distinct loci in
    contact - a contact that does not exist. Merge first and it stays one
    fragment, because the block carries the best MAPQ of its parts.
    """
    monomers = [("chr1", 1000, 1200, False, 60),
                ("chr1", 1200, 1400, False, 0),     # unreliable, in the middle
                ("chr1", 1400, 1600, False, 60),
                ("chr2", 50000, 50200, False, 60)]
    bam = write_bam(tmp_path / "chain.bam", [("readB", monomers)])
    out = tmp_path / "f.tsv"
    main([bam, "--output", str(out), "--mapq", "30"])
    rows = read_tsv(out)
    assert len(rows) == 2, "must be one chr1 fragment and one chr2 fragment"
    chr1 = next(r for r in rows if r["chrom"] == "chr1")
    assert (int(chr1["start"]), int(chr1["end"])) == (1000, 1600)
    assert chr1["n_monomers_merged"] == "3"


def test_stricter_mapq_never_creates_fragments(tmp_path):
    """A stricter filter cannot yield *more* pieces. It once did.

    That impossibility is what revealed the merge/filter order was wrong: the
    molecule count went from 61,677 to 70,735 as the threshold was tightened.
    Filtering first removes a middle monomer, the chain around it splits, and
    one fragment becomes two.
    """
    records = []
    for r in range(30):
        # every other monomer is unreliable and sits inside a contiguous run
        mons = [("chr1", 1000 + 200 * i, 1200 + 200 * i, False,
                 0 if i % 2 else 60) for i in range(7)]
        mons.append(("chr2", 50000 + r, 50200 + r, False, 60))
        records.append((f"read{r:03d}", mons))
    bam = write_bam(tmp_path / "many.bam", records)
    counts = []
    for q in (0, 1, 30, 60):
        out = tmp_path / f"f{q}.tsv"
        main([bam, "--output", str(out), "--mapq", str(q)])
        counts.append(len(read_tsv(out)))
    assert counts == sorted(counts, reverse=True), \
        f"fragment count rose with a stricter MAPQ: {counts}"


# --------------------------------------------------------------------------
# options and guards
# --------------------------------------------------------------------------
def test_min_fragments_filters_molecules(tmp_path, split_locus_bam):
    out = tmp_path / "f.tsv"
    main([split_locus_bam, "--output", str(out), "--min-fragments", "3"])
    assert read_tsv(out) == [], "the molecule has 2 fragments after merging"


def test_min_length_drops_short_fragments(tmp_path, split_locus_bam):
    out = tmp_path / "f.tsv"
    main([split_locus_bam, "--output", str(out), "--min-length", "500"])
    rows = read_tsv(out)
    assert [r["chrom"] for r in rows] == ["chr1"], "chr2 piece is 200 bp"


def test_min_sep_drops_close_cis_pairs(tmp_path):
    monomers = [("chr1", 1000, 1200, False, 60),
                ("chr1", 6000, 6200, False, 60)]
    bam = write_bam(tmp_path / "cis.bam", [("readC", monomers)])
    out, pairs = tmp_path / "f.tsv", tmp_path / "c.pairs"
    main([bam, "--output", str(out), "--pairs", str(pairs),
          "--min-sep", "10000"])
    assert [ln for ln in pairs.read_text().splitlines()
            if not ln.startswith("#")] == []


def test_stats_json_reports_the_collapse(tmp_path, split_locus_bam):
    out, stats = tmp_path / "f.tsv", tmp_path / "s.json"
    main([split_locus_bam, "--output", str(out), "--stats", str(stats)])
    data = json.loads(stats.read_text())
    assert data["aligned_monomers"] == 5
    assert data["fragments_after_merge"] == 2
    assert data["collapse_ratio"] == 2.5


def test_gzip_output_is_written(tmp_path, split_locus_bam):
    out = tmp_path / "f.tsv.gz"
    main([split_locus_bam, "--output", str(out)])
    assert out.read_bytes()[:2] == b"\x1f\x8b"
    assert len(read_tsv(out)) == 2


def test_coordinate_sorted_bam_is_rejected(tmp_path, split_locus_bam):
    cs = tmp_path / "cs.bam"
    pysam.sort("-o", str(cs), split_locus_bam)
    with pytest.raises(SystemExit, match="coordinate-sorted"):
        list(iter_concatemers(str(cs)))


def test_unmapped_and_secondary_records_are_ignored(tmp_path):
    bam = write_bam(tmp_path / "x.bam",
                    [("readD", [("chr1", 1000, 1200, False, 60),
                                ("chr2", 50000, 50200, False, 60)])])
    with pysam.AlignmentFile(bam, check_sq=False) as f:
        header, reads = f.header, list(f)
    with pysam.AlignmentFile(str(tmp_path / "y.bam"), "wb",
                             header=header) as out:
        for r in reads:
            out.write(r)
        junk = pysam.AlignedSegment(header)
        junk.query_name = "readD:9999:9999"
        junk.query_sequence = "AAAA"
        junk.flag = 4
        junk.set_tag("MI", "readD")
        out.write(junk)
    got = list(iter_concatemers(str(tmp_path / "y.bam")))
    assert len(got) == 1 and len(got[0][1]) == 2

"""The junction-motif analysis behind `porec-junctions`.

This is the command the digest report points at when it says it cannot tell
whether an enzyme cut. It therefore has to give the right answer on data where
the answer is known, which is what these tests build: a synthetic genome cut at
chosen DpnII sites, joined into concatemers, and aligned.
"""
import array
import random

import pysam
import pytest

from porec_tools.enzymes import resolve_enzymes
from porec_tools.junctions import (
    junction_boundaries,
    main,
    motif_at,
    read_order_key,
)

CHROMS = {"chr1": 60000, "chr2": 60000}
FRAG = 400          # distance between the cut sites we plant


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """A reference with DpnII sites every FRAG bases and no AAGCTT anywhere."""
    rng = random.Random(11)
    path = tmp_path_factory.mktemp("ref") / "ref.fa"
    seqs = {}
    for name, length in CHROMS.items():
        bases = [rng.choice("ACGT") for _ in range(length)]
        for i in range(0, length - FRAG, FRAG):
            bases[i:i + 4] = list("GATC")          # plant the cut sites
        seq = "".join(bases).replace("AAGCTT", "AAGCTA")   # no HindIII site
        seqs[name] = seq
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(path))
    return str(path), seqs


def write_bam(path, records):
    """records: list of concatemers, each a list of (chrom, start, end, rev)."""
    header = {"HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
              "SQ": [{"SN": n, "LN": ln} for n, ln in CHROMS.items()]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, monomers in enumerate(records):
            offset = 0
            for j, (chrom, start, end, rev) in enumerate(monomers):
                a = pysam.AlignedSegment(out.header)
                a.query_name = f"read{i:04d}:{offset:04d}:{offset + end - start:04d}"
                a.query_sequence = "A" * (end - start)
                a.query_qualities = pysam.qualitystring_to_array(
                    "I" * (end - start))
                a.reference_id = out.header.get_tid(chrom)
                a.reference_start = start
                a.cigar = [(0, end - start)]
                a.mapping_quality = 60
                a.flag = 16 if rev else 0
                a.set_tag("MI", f"read{i:04d}")
                a.set_tag("Xc", array.array(
                    "H", [offset, offset + end - start, 0, j]))
                out.write(a)
                offset += end - start
    return str(path)


@pytest.fixture(scope="module")
def cut_bam(tmp_path_factory, genome):
    """Concatemers whose monomers start and end exactly on planted GATC sites."""
    rng = random.Random(3)
    records = []
    for _ in range(400):
        monomers = []
        for _ in range(3):
            chrom = rng.choice(list(CHROMS))
            k = rng.randrange(1, CHROMS[chrom] // FRAG - 2)
            start = k * FRAG                        # a GATC cut point
            monomers.append((chrom, start, start + FRAG, rng.random() < 0.5))
        records.append(monomers)
    return write_bam(tmp_path_factory.mktemp("cut") / "cut.bam", records)


@pytest.fixture(scope="module")
def random_bam(tmp_path_factory, genome):
    """Same shape, but boundaries fall wherever - no enzyme made them."""
    rng = random.Random(4)
    records = []
    for _ in range(400):
        monomers = []
        for _ in range(3):
            chrom = rng.choice(list(CHROMS))
            start = rng.randrange(1000, CHROMS[chrom] - 2000)
            monomers.append((chrom, start, start + FRAG, rng.random() < 0.5))
        records.append(monomers)
    return write_bam(tmp_path_factory.mktemp("rnd") / "rnd.bam", records)


def run(capsys, bam, ref, *extra):
    assert main([bam, ref, "--enzymes", "DpnII,HindIII", *extra]) == 0
    return capsys.readouterr().out


def parse(out):
    rows = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) >= 6 and parts[0] in ("DpnII", "HindIII", "NlaIII"):
            rows[parts[0]] = {
                "at_junctions": float(parts[2].rstrip("%")),
                "enrichment": float(parts[4].rstrip("x")),
                "verdict": " ".join(parts[5:]),
            }
    return rows


def test_detects_the_enzyme_that_made_the_junctions(capsys, genome, cut_bam):
    rows = parse(run(capsys, cut_bam, genome[0]))
    assert rows["DpnII"]["at_junctions"] > 90
    assert rows["DpnII"]["enrichment"] > 5
    assert rows["DpnII"]["verdict"] == "yes, main enzyme"


def test_absent_enzyme_is_called_absent(capsys, genome, cut_bam):
    """HindIII has no site in this genome at all: it must not be credited."""
    rows = parse(run(capsys, cut_bam, genome[0]))
    assert rows["HindIII"]["at_junctions"] == 0.0
    assert rows["HindIII"]["verdict"] == "no"


def test_random_boundaries_give_no_enrichment(capsys, genome, random_bam):
    """The negative control: junctions no enzyme made must not look enzymatic.

    Without this, a test that only ever sees real cut sites would pass just as
    happily on a script that always says yes.
    """
    rows = parse(run(capsys, random_bam, genome[0]))
    assert rows["DpnII"]["enrichment"] < 2.5
    assert rows["DpnII"]["verdict"] in ("no", "unclear")


def test_reverse_strand_monomers_are_handled(capsys, genome, tmp_path):
    """A reverse alignment's read-5' end is its reference_end, not its start."""
    recs = [[("chr1", k * FRAG, k * FRAG + FRAG, True),
             ("chr2", (k + 7) * FRAG, (k + 7) * FRAG + FRAG, True)]
            for k in range(5, 100)]
    bam = write_bam(tmp_path / "rev.bam", recs)
    rows = parse(run(capsys, bam, genome[0]))
    assert rows["DpnII"]["at_junctions"] > 90


def test_read_ends_are_excluded(genome, cut_bam):
    """n monomers give n-1 junctions, so 2(n-1) boundaries - never 2n."""
    got = list(junction_boundaries(cut_bam, mapq=20, min_jump=1000))
    bam = pysam.AlignmentFile(cut_bam, "rb")
    n_monomers = sum(1 for _ in bam)
    bam.close()
    assert len(got) <= 2 * (n_monomers - n_monomers // 3)


def test_close_boundaries_are_not_junctions(genome, tmp_path):
    """Adjacent fragments of one locus are uncut sites, not ligation events."""
    recs = [[("chr1", 10000, 10400, False), ("chr1", 10400, 10800, False)]]
    bam = write_bam(tmp_path / "near.bam", recs)
    assert list(junction_boundaries(bam, mapq=20, min_jump=1000)) == []


def test_mapq_filter_applies(genome, cut_bam):
    strict = list(junction_boundaries(cut_bam, mapq=61, min_jump=1000))
    assert strict == []


def test_read_order_from_xc_tag(genome, cut_bam):
    bam = pysam.AlignmentFile(cut_bam, "rb")
    first = next(iter(bam))
    assert read_order_key(first) == first.get_tag("Xc")[0]
    bam.close()


def test_read_order_falls_back_to_the_name():
    a = pysam.AlignedSegment()
    a.query_name = "somereadname:0348:0547"
    assert read_order_key(a) == 348


def test_motif_at_respects_tolerance(genome):
    ref = pysam.FastaFile(genome[0])
    dpn = resolve_enzymes("DpnII")[0]
    assert motif_at(ref, "chr1", 2 * FRAG, dpn, tol=0)
    assert not motif_at(ref, "chr1", 2 * FRAG + 40, dpn, tol=0)
    assert motif_at(ref, "chr1", 2 * FRAG + 3, dpn, tol=4)
    ref.close()


def test_bad_enzyme_name_exits_cleanly(capsys, genome, cut_bam):
    assert main([cut_bam, genome[0], "--enzymes", "Nonsense9"]) == 2
    assert "Unknown restriction enzyme" in capsys.readouterr().err


def test_missing_reference_exits_cleanly(capsys, cut_bam, tmp_path):
    assert main([cut_bam, str(tmp_path / "nope.fa"), "--enzymes", "DpnII"]) == 2
    assert "faidx" in capsys.readouterr().err


def test_no_junction_found_is_reported(capsys, genome, tmp_path):
    bam = write_bam(tmp_path / "single.bam", [[("chr1", 4000, 4400, False)]])
    assert main([bam, genome[0], "--enzymes", "DpnII"]) == 1
    assert "No ligation junction" in capsys.readouterr().err


def test_result_is_reproducible(capsys, genome, cut_bam):
    a = parse(run(capsys, cut_bam, genome[0], "--seed", "1"))
    b = parse(run(capsys, cut_bam, genome[0], "--seed", "1"))
    assert a == b

"""Equivalence with upstream ``pore-c-py`` in single-enzyme mode.

This is the load-bearing test of the whole project. The claim being made is:

    with one enzyme, this tool produces exactly what pore-c-py produces.

If that holds, then a user switching to ``pore-c-aqb`` risks nothing, and the
multi-enzyme behaviour is a strict extension rather than a reimplementation.

Two levels of checking:

* :func:`test_matches_reference_implementation` re-implements upstream's
  algorithm inline (it is short) and compares record by record. Always runs.
* :func:`test_matches_upstream_container` runs the *real* pore-c-py inside the
  published wf-pore-c image and diffs the BAMs. Skipped when Docker or the
  image is unavailable, so CI stays green on machines without it.

The image is deliberately ``wf-pore-c``, which carries pore-c-py **2.0.6** —
the version the workflow actually runs, and the one ``_vendored.py`` was taken
from. The standalone ``ontresearch/pore-c-py`` image is a later 2.1.x that
changed the mod-base tags (``ML`` became a uint8 array and ``MN`` was added);
comparing against it reports differences that are upstream's version bump, not
a defect here.
"""
from __future__ import annotations

import array
import copy
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest
from Bio import Restriction
from Bio.Seq import Seq

from pore_c_aqb import _vendored
from pore_c_aqb.digest import digest_sequence
from pore_c_aqb.enzymes import resolve_enzymes

from conftest import make_read, random_seq, write_bam

UPSTREAM_IMAGE = (
    "ontresearch/wf-pore-c:sha3787c234c0cacf66a67fb77da223cc2e1cb0baf0"
)


# --------------------------------------------------------------------------
# upstream algorithm, transcribed from pore_c_py/digest.py 2.0.6
# --------------------------------------------------------------------------
def upstream_digest(align, enzyme_name):
    """Exactly pore-c-py's digest_sequence, for one enzyme."""
    enzyme = getattr(Restriction, enzyme_name)
    tags_remove = {'mv'}
    concatemer_id = align.query_name
    cut_points = [x - 1 for x in enzyme.search(Seq(align.query_sequence))]
    read_length = len(align.query_sequence)
    num_digits = len(str(read_length))
    intervals = _vendored.splits_to_intervals(cut_points, read_length)
    num_intervals = len(intervals)

    for idx, (start, end) in enumerate(intervals):
        read = copy.copy(align)
        seq = align.query_sequence[start:end]
        qual = align.query_qualities[start:end] if align.query_qualities else None
        read.query_sequence = seq
        read.query_qualities = qual
        if ('Mm' in tags_remove) or ('MM' in tags_remove):
            for tag in ('Mm', 'Ml', 'MM', 'ML'):
                read.set_tag(tag, None)
        else:
            mm, ml = _vendored.get_subread_modified_bases(align, start, end)
            for tag in ('Mm', 'Ml'):
                read.set_tag(tag, None)
            read.set_tag("MM", mm)
            read.set_tag("ML", ml)
        read.query_name = \
            f"{concatemer_id}:{start:0{num_digits}d}:{end:0{num_digits}d}"
        read.set_tag(_vendored.MONOMER_DATA_TAG,
                     [start, end, read_length, idx, num_intervals])
        read.set_tag(_vendored.CONCATEMER_ID_TAG, concatemer_id, "Z")
        yield read


def _comparable(read):
    """The fields that must match, in a form pytest can diff readably."""
    return (
        read.query_name,
        read.query_sequence,
        None if read.query_qualities is None
        else pysam.qualities_to_qualitystring(read.query_qualities),
        read.get_tag("MI"),
        tuple(read.get_tag("Xc")),
        read.get_tag("MM") if read.has_tag("MM") else None,
        read.get_tag("ML") if read.has_tag("ML") else None,
    )


def _every_field(read):
    """Everything, including tag value types.

    ``_comparable`` names the fields we care about; this one takes no view and
    so catches a divergence in a tag nobody thought to list - including the
    *type* a tag is written with, which is how an ML written as a string
    instead of a uint8 array would show up.
    """
    return (
        read.query_name,
        read.query_sequence,
        None if read.query_qualities is None
        else pysam.qualities_to_qualitystring(read.query_qualities),
        read.flag,
        sorted((t[0], repr(t[1]), t[2])
               for t in read.get_tags(with_value_type=True)),
    )


def methylated_read(name, seq):
    """A read carrying MM/ML tags, so mod-base handling is exercised.

    Without this the container comparison runs on tag-free synthetic reads and
    would pass even if the MM/ML recomputation were completely broken.
    """
    positions = [i for i, b in enumerate(seq) if b == "C"]
    if not positions:
        raise ValueError("sequence has no C to modify")
    deltas = ",".join("0" for _ in positions)
    read = make_read(name, seq)
    read.set_tag("MM", f"C+m,{deltas};", "Z")
    read.set_tag("ML", array.array("B", [(7 * i) % 256 for i in
                                         range(len(positions))]))
    return read


@pytest.mark.parametrize("enzyme", ["DpnII", "NlaIII", "HinfI", "HindIII"])
@pytest.mark.parametrize("seed", range(6))
def test_matches_reference_implementation(enzyme, seed):
    seq = random_seq(400 + seed * 137, seed=seed)
    read = make_read(f"read{seed}", seq)

    ours = [_comparable(r) for r in
            digest_sequence(copy.copy(read), resolve_enzymes(enzyme))]
    theirs = [_comparable(r) for r in upstream_digest(copy.copy(read), enzyme)]
    assert ours == theirs


def test_matches_reference_with_modified_bases():
    """Methylation tags must survive trimming identically.

    The MM/ML recomputation is the subtlest part of the digest: an off-by-one
    there would silently corrupt every downstream methylation analysis.
    """
    seq = random_seq(600, seed=42)
    read = make_read("modread", seq)
    # mark a handful of Cs as 5mC, in the interim tag form the code upgrades
    c_positions = [i for i, b in enumerate(seq) if b == "C"][:20]
    deltas, prev = [], 0
    all_c = [i for i, b in enumerate(seq) if b == "C"]
    for pos in c_positions:
        rank = all_c.index(pos)
        deltas.append(rank - prev)
        prev = rank + 1
    read.set_tag("MM", "C+m," + ",".join(map(str, deltas)) + ";", "Z")
    read.set_tag("ML", array.array("B", [200] * len(c_positions)))

    ours = [_comparable(r) for r in
            digest_sequence(copy.copy(read), resolve_enzymes("DpnII"))]
    theirs = [_comparable(r) for r in upstream_digest(copy.copy(read), "DpnII")]
    assert ours == theirs
    assert any(t[5] for t in ours), "expected at least one MM tag to be set"


def test_multi_enzyme_is_a_strict_refinement():
    """Adding an enzyme only ever splits monomers further, never merges."""
    seq = random_seq(2000, seed=11)
    read = make_read("r", seq)
    one = [tuple(r.get_tag("Xc")[:2]) for r in
           digest_sequence(copy.copy(read), resolve_enzymes("DpnII"))]
    two = [tuple(r.get_tag("Xc")[:2]) for r in
           digest_sequence(copy.copy(read), resolve_enzymes("DpnII,NlaIII"))]
    # every boundary of the single-enzyme digest survives in the two-enzyme one
    one_bounds = {b for span in one for b in span}
    two_bounds = {b for span in two for b in span}
    assert one_bounds <= two_bounds


# --------------------------------------------------------------------------
# against the real upstream tool
# --------------------------------------------------------------------------
def _docker_available() -> bool:
    if shutil.which("docker") is None:
        return False
    try:
        subprocess.run(
            ["docker", "image", "inspect", UPSTREAM_IMAGE],
            check=True, capture_output=True, timeout=120)
        return True
    except Exception:
        return False


def _run_upstream(tmp_path: Path, in_bam: str, enzyme: str) -> Path:
    theirs = tmp_path / "theirs.bam"
    subprocess.run(
        ["docker", "run", "--rm", "-v", f"{tmp_path}:/data", UPSTREAM_IMAGE,
         "bash", "-lc",
         f"pore-c-py digest /data/{Path(in_bam).name} {enzyme} "
         f"--output /data/theirs.bam"],
        check=True, capture_output=True, timeout=900)
    return theirs


@pytest.mark.skipif(not _docker_available(),
                    reason="Docker or the upstream wf-pore-c image unavailable")
@pytest.mark.parametrize("enzyme", ["DpnII", "NlaIII"])
def test_matches_upstream_container_with_modified_bases(tmp_path: Path, enzyme):
    """The comparison that matters: real MM/ML tags, every field, every tag.

    MM/ML recomputation is the subtlest part of the digest and the part most
    likely to drift from upstream, so it is checked against the real tool
    rather than against our own transcription of it.
    """
    reads = [methylated_read(f"read{i}", random_seq(400 + 91 * i, seed=i))
             for i in range(20)]
    in_bam = write_bam(tmp_path / "mod.bam", reads)

    ours = tmp_path / "ours.bam"
    subprocess.run(
        [sys.executable, "-m", "pore_c_aqb.cli", "digest", enzyme, in_bam,
         "--output", str(ours), "--quiet"],
        check=True, capture_output=True)
    theirs = _run_upstream(tmp_path, in_bam, enzyme)

    with pysam.AlignmentFile(str(ours), "rb", check_sq=False) as a, \
         pysam.AlignmentFile(str(theirs), "rb", check_sq=False) as b:
        got = [_every_field(r) for r in a]
        want = [_every_field(r) for r in b]

    assert len(got) == len(want), "different number of monomers"
    assert got == want
    assert any(t[0] == "ML" for row in want for t in row[4]), \
        "the fixture must actually carry ML tags for this test to mean anything"


@pytest.mark.skipif(not _docker_available(),
                    reason="Docker or the upstream wf-pore-c image unavailable")
@pytest.mark.parametrize("enzyme", ["DpnII", "NlaIII"])
def test_matches_upstream_container(tmp_path: Path, enzyme):
    """Diff our BAM against one produced by the published pore-c-py."""
    reads = [make_read(f"read{i}", random_seq(300 + 89 * i, seed=i))
             for i in range(25)]
    in_bam = write_bam(tmp_path / "in.bam", reads)

    ours = tmp_path / "ours.bam"
    subprocess.run(
        [sys.executable, "-m", "pore_c_aqb.cli", "digest", enzyme, in_bam,
         "--output", str(ours), "--quiet"],
        check=True, capture_output=True)

    theirs = tmp_path / "theirs.bam"
    subprocess.run(
        ["docker", "run", "--rm", "-v", f"{tmp_path}:/data", UPSTREAM_IMAGE,
         "bash", "-lc",
         f"pore-c-py digest /data/in.bam {enzyme} "
         f"--output /data/theirs.bam"],
        check=True, capture_output=True, timeout=900)

    with pysam.AlignmentFile(str(ours), "rb", check_sq=False) as a, \
         pysam.AlignmentFile(str(theirs), "rb", check_sq=False) as b:
        got = [_comparable(r) for r in a]
        want = [_comparable(r) for r in b]

    assert len(got) == len(want), "different number of monomers"
    assert got == want

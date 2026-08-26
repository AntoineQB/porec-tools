"""Shared fixtures: synthetic concatemers with known cut sites."""
from __future__ import annotations

import random

import pysam
import pytest

HEADER = {"HD": {"VN": "1.6", "SO": "unknown"}, "SQ": []}


def make_read(name: str, seq: str, qual: str | None = None,
              tags: list | None = None) -> pysam.AlignedSegment:
    """An unaligned BAM record carrying `seq`."""
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = 4
    a.query_qualities = pysam.qualitystring_to_array(qual or "I" * len(seq))
    if tags:
        for tag, value, vtype in tags:
            a.set_tag(tag, value, vtype)
    return a


def write_bam(path, reads) -> str:
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as fh:
        for r in reads:
            fh.write(r)
    return str(path)


def random_seq(n: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def spaced_sites():
    """A sequence with GATC and CATG sites at positions we control.

    Layout (0-based):
        0    filler 20
        20   GATC        <- DpnII cuts at 20
        60   CATG        <- NlaIII cuts at 64 (CATG^, fst5=4)
        100  GATC        <- DpnII cuts at 100
    """
    filler = random_seq(200, seed=7)
    # scrub any accidental sites from the filler so positions are exact
    scrubbed = (filler.replace("GATC", "GTTC").replace("CATG", "CTTG")
                      .replace("GATC", "GTTC").replace("CATG", "CTTG"))
    s = list(scrubbed[:140])
    s[20:24] = list("GATC")
    s[60:64] = list("CATG")
    s[100:104] = list("GATC")
    seq = "".join(s)
    assert "GATC" in seq and "CATG" in seq
    return seq

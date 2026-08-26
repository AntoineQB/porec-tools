# This Source Code Form is subject to the terms of the Oxford Nanopore
# Technologies PLC. Public License Version 1.0. See the LICENSE file.
"""Multi-enzyme digestion of unaligned concatemers.

Derived from ``pore_c_py.digest`` (Oxford Nanopore Technologies PLC).

The only behavioural change is that :func:`find_cut_points` takes a *list* of
enzymes and returns the sorted union of their cut positions. Everything
downstream — interval splitting, sequence and quality trimming, MM/ML tag
recomputation, monomer naming and tagging — is unchanged, which is what lets
:mod:`tests.test_equivalence` assert byte-identical output against the original
in the single-enzyme case.

Why the union, and not one pass per enzyme
------------------------------------------
A concatemer is digested once, by all the enzymes present in the reaction at
the same time. The resulting fragments are therefore delimited by the union of
all recognition sites. Digesting separately and merging afterwards would be
wrong: it would produce overlapping monomer sets rather than one partition.
"""
from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterator, Sequence

from Bio.Seq import Seq

from pore_c_aqb import _vendored
from pore_c_aqb.enzymes import ResolvedEnzyme
from pore_c_aqb.sites import find_cuts_for_enzyme

__all__ = ["DigestStats", "find_cut_points", "digest_sequence",
           "get_concatemer_seqs"]


@dataclass
class DigestStats:
    """Counters gathered while digesting.

    ``cuts_per_enzyme`` answers a question that is otherwise very hard to ask:
    *did this enzyme actually cut?* A reaction where one enzyme failed produces
    almost no cuts attributable to it, and that shows up here immediately
    rather than as a silent loss of contacts downstream.
    """

    n_concatemers: int = 0
    n_monomers: int = 0
    n_cut_points: int = 0
    cuts_per_enzyme: Counter = field(default_factory=Counter)
    #: cut positions found by more than one enzyme, counted once in the digest
    n_shared_cut_points: int = 0

    def summary_lines(self, enzymes: Sequence[ResolvedEnzyme]) -> list[str]:
        """Human-readable report, one line per enzyme."""
        out = [
            f"Digested {self.n_concatemers:,} concatemers into "
            f"{self.n_monomers:,} monomers "
            f"({self.n_cut_points:,} cut points)."
        ]
        if not self.n_cut_points:
            return out
        for enz in enzymes:
            n = self.cuts_per_enzyme.get(enz.name, 0)
            pct = 100.0 * n / self.n_cut_points
            flag = "   <-- no sites found, check the protocol" if n == 0 else ""
            out.append(
                f"  {enz.name:<12} {enz.site:<10} {n:>12,} cuts "
                f"({pct:5.1f}%){flag}"
            )
        if self.n_shared_cut_points:
            out.append(
                f"  {'(shared)':<12} {'':<10} {self.n_shared_cut_points:>12,} "
                f"positions cut by more than one enzyme, counted once"
            )
        return out


def find_cut_points(
    sequence: str | Seq,
    enzymes: Sequence[ResolvedEnzyme],
    stats: DigestStats | None = None,
) -> list[int]:
    """Sorted, de-duplicated 0-based cut positions for all enzymes.

    Positions follow ``pore-c-py``'s convention: 0-based, the cut falling
    immediately before the returned index. Site location is delegated to
    :func:`pore_c_aqb.sites.find_cuts_for_enzyme`, which is reproducible
    across Biopython versions for palindromic enzymes.

    Two enzymes can cut at the same position. That is one cut in the tube, so
    it must be one cut here — hence the set.
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence)

    cuts: set[int] = set()
    per_enzyme: dict[str, set[int]] = {}
    for enz in enzymes:
        found = find_cuts_for_enzyme(sequence, enz)
        per_enzyme[enz.name] = found
        cuts |= found

    if stats is not None:
        for name, found in per_enzyme.items():
            stats.cuts_per_enzyme[name] += len(found)
        # a position claimed by k enzymes is counted k times above but once
        # in the digest; record the excess so the report stays honest
        total_claims = sum(len(v) for v in per_enzyme.values())
        stats.n_shared_cut_points += total_claims - len(cuts)
        stats.n_cut_points += len(cuts)

    return sorted(cuts)


def digest_sequence(align, enzymes: Sequence[ResolvedEnzyme],
                    tags_remove=None, stats: DigestStats | None = None):
    """Split one concatemer into monomers. Mirrors pore-c-py exactly."""
    # the move tag massively bloats files, and we don't care for
    # it or handle it in trimming, so force its removal by default.
    if tags_remove is None:
        tags_remove = {'mv'}

    concatemer_id = align.query_name
    cut_points = find_cut_points(align.query_sequence, enzymes, stats)
    read_length = len(align.query_sequence)
    num_digits = len(str(read_length))
    intervals = _vendored.splits_to_intervals(cut_points, read_length)
    num_intervals = len(intervals)

    for idx, (start, end) in enumerate(intervals):
        read = copy.copy(align)
        # trim the sequence and quality
        seq = align.query_sequence[start:end]
        qual = None
        if align.query_qualities:
            qual = align.query_qualities[start:end]
        read.query_sequence = seq
        read.query_qualities = qual
        # deal with mods, upgrading tag from interim to approved spec
        if ('Mm' in tags_remove) or ('MM' in tags_remove):
            for tag in ('Mm', 'Ml', 'MM', 'ML'):
                read.set_tag(tag, None)
        else:
            mm, ml = _vendored.get_subread_modified_bases(align, start, end)
            for tag in ('Mm', 'Ml'):
                read.set_tag(tag, None)
            read.set_tag("MM", mm)
            read.set_tag("ML", ml)
        # lexographically sortable monomer ID
        read.query_name = \
            f"{concatemer_id}:{start:0{num_digits}d}:{end:0{num_digits}d}"
        read.set_tag(
            _vendored.MONOMER_DATA_TAG,
            [start, end, read_length, idx, num_intervals])
        _vendored.set_monomer_data(
            read, start, end, read_length, idx, num_intervals)
        read.set_tag(_vendored.CONCATEMER_ID_TAG, concatemer_id, "Z")
        yield read


def get_concatemer_seqs(
    input_file,
    enzymes: Sequence[ResolvedEnzyme],
    remove_tags=None,
    stats: DigestStats | None = None,
) -> Iterator:
    """Digest concatemers into unaligned monomers.

    :param input_file: pysam.AlignmentFile input
    :param enzymes: enzymes resolved by :func:`pore_c_aqb.enzymes.resolve_enzymes`
    :param remove_tags: additional SAM tags to strip
    :param stats: optional :class:`DigestStats` to accumulate counters into
    """
    if stats is None:
        stats = DigestStats()
    tags_remove = {"mv"}
    if remove_tags:
        tags_remove.update(set(remove_tags))
    for align in input_file.fetch(until_eof=True):
        stats.n_concatemers += 1
        for read in digest_sequence(align, enzymes, tags_remove, stats):
            stats.n_monomers += 1
            yield read
